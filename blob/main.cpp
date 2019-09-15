#include <cstdio>
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "gdal/gdal_priv.h"
#include "gdal/ogrsf_frmts.h"
#include "gtest/gtest.h"

#include "admin_patch.h"
#include "component_data.h"
#include "gdal_raster.h"
#include "gdal_vector.h"
#include "on_demand_raster.h"
#include "projection.h"
#include "sparse_settlements.h"


using namespace dd_harp;
using namespace std;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

/* An Alpha_shape_2 is a Delaunay_triangulation_2 which is a Triangulation_2.
 * .classify(const Point&) returns the Classification_type.
 * Classification types are EXTERIOR (not in alpha-complex).
 * SINGULAR - boundary, not in a 2d face, REGULAR, edge of 2D face, INTERIOR of face.
 * The segment() accessor is in CGAL/Triangulation_2.h.
 * A Segment has .source() and .target() returning a Point_2().
 * It also has .vertex(0) and .vertex(1) returning source and target vertices.
 * A Point_2 has .x() and .y()
 */


po::options_description parser(const map<string,fs::path>& path_argument)
{
    po::options_description options("blob create");
    options.add_options()
            ("help", "write help message")
            ("test", "run all tests")
            ("tile-subset", po::value<int>(), "how many tiles to use")
            ("population-cutoff", po::value<double>(), "minimum people per pixel")
            ("population-per-patch", po::value<double>(), "number of people in each patch")
            ("admin-limit", po::value<int>(), "an upper limit on how many admin units to process")
            ;
    for (auto const& [name, default_path] : path_argument) {
        options.add_options()(name.c_str(), po::value<fs::path>(), name.c_str());
    }
    return options;
}


/*! Given command line args, assign file paths, using defaults.
 *
 * @param vm Variables initialized from the command line.
 * @param input_path A map from a keyword to a default filename for each file type.
 * @return Whether all the input files specified as defaults or on command line do exist.
 */
bool read_paths_from_command_line_args(po::variables_map& vm, map<string,fs::path>& input_path) {
    for (auto path_iter=input_path.begin(); path_iter != input_path.end(); ++path_iter) {
        if (vm.count(path_iter->first)) {
            path_iter->second = vm[path_iter->first].as<fs::path>();
        }
        if (!fs::exists(path_iter->second)) {
            cout << "path for " << path_iter->first << " not found: " << path_iter->second << endl;
            return false;
        }
    }
    return true;
}


/*! A sanity check for longitude and latitude.
 *  Uganda is near longitude=33, latitude=2
 */
void assert_x_is_longitude(const OGRMultiPolygon * geometry) {
    OGREnvelope bounds;
    geometry->getEnvelope(&bounds);
    if (bounds.MinX < 10 || bounds.MinY > 10) {
        throw std::runtime_error("Polygon X is not longitude for Uganda.");
    }
}


// Reads a list of points and returns a list of segments
// corresponding to the Alpha shape.

int entry(int argc, char* argv[])
{
    ios_base::sync_with_stdio(false);
    const char* user = getenv("USER");
    string home{"/home/"};
    if (user != nullptr) {
        home += user;
    } else {
        home += "adolgert";
    }
    map<string,fs::path> input_path = {
            {"settlement", home + "/dev/spacepop/data/hrsl/hrsl_uga_pop.tif"},
            {"pfpr", home + "/dev/spacepop/data/PfPR/Raster Data/PfPR_rmean/2019_Global_PfPR_2017.tif"},
            {"admin", home + "/dev/spacepop/data/uga_admbnda/uga_admbnda_adm3_UBOS_v5.shp"},
    };
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, parser(input_path)), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << parser(input_path) << endl;
        return 0;
    }
    if (!read_paths_from_command_line_args(vm, input_path)) {
        return 3;
    }

    int use_subset_of_tiles = 0;  // Use only a few tiles in the corner.
    if (vm.count("tile-subset")) {
        use_subset_of_tiles = vm["tile-subset"].as<int>();
    }
    double population_cutoff = 0.1;
    if (vm.count("population-cutoff")) {
        population_cutoff = vm["population-cutoff"].as<double>();
    }
    double population_per_patch = 500;
    if (vm.count("population-per-patch")) {
        population_per_patch = vm["population-per-patch"].as<double>();
    }
    int admin_limit{std::numeric_limits<int>::max()};
    if (vm.count("admin-limit")) {
        admin_limit = vm["admin-limit"].as<int>();
    }
    if (vm.count("test")) {
        ::testing::InitGoogleTest(&argc, argv);
        return RUN_ALL_TESTS();
    }

    // This initializes GDAL's list of drivers to read and write files.
    GDALAllRegister();
    cout << "hrsl " << input_path.at("settlement") << endl;
    auto settlement_dataset = OpenGeoTiff(input_path.at("settlement"));
    GDALRasterBand* settlement_pop_band = settlement_dataset->GetRasterBand(1);
    vector<double> settlement_geo_transform(6);
    if (settlement_dataset->GetGeoTransform(&settlement_geo_transform[0]) != CE_None) {
        cout << "Could not get settlement transform" << endl;
        return 9;
    } else {
        cout << "settlement geotransform ";
        for (auto sgt: settlement_geo_transform) {
            cout << sgt << " ";
        }
        cout << endl;
    }

    auto pfpr_dataset = OpenGeoTiff(input_path.at("pfpr"));
    GDALRasterBand* pfpr_band = pfpr_dataset->GetRasterBand(1);
    vector<double> pfpr_geo_transform(6);
    if (pfpr_dataset->GetGeoTransform(&pfpr_geo_transform[0]) != CE_None) {
        cout << "Could not get pfpr transform" << endl;
        return 10;
    }

    auto projection_ref = settlement_dataset->GetProjectionRef();
    OGRSpatialReference lat_long_srs(projection_ref);
    auto to_projected = reproject(&lat_long_srs);

    auto admin_dataset = GDALDatasetUniquePtr(static_cast<GDALDataset *>(GDALOpenEx(
            input_path.at("admin").c_str(),
            GDAL_OF_VECTOR | GDAL_OF_READONLY | GDAL_OF_VERBOSE_ERROR,
            nullptr,
            nullptr,
            nullptr
    )), GDALDatasetUniquePtrDeleter());
    if (admin_dataset == nullptr) {
        cout << "Could not load dataset from " << input_path.at("admin") << endl;
        return 7;
    }
    if (admin_dataset->GetLayerCount() == 0) {
        cout << "Admin dataset has no layers" << endl;
        return 8;
    }

    // Loop over the admin layers, making patches within each layer.
    auto settlement_arr = OnDemandRaster(settlement_pop_band, settlement_geo_transform);
    auto pfpr_arr = OnDemandRaster(pfpr_band, pfpr_geo_transform);

    OGRLayer* first_admin_layer = *admin_dataset->GetLayers().begin();
    OGRMultiPolygon mp_buffer;
    int poly_idx{0};
    vector<ComponentData> components;
    first_admin_layer->ResetReading();
    OGRFeature* geom_feature = first_admin_layer->GetNextFeature();
    for (
            int feature_idx=0;
            (geom_feature != nullptr) && (feature_idx != admin_limit);
            ++feature_idx, geom_feature = first_admin_layer->GetNextFeature()
                    ) {
        cout << "polygon " << poly_idx << endl;
        auto geometry = geom_feature->GetGeometryRef();  // reference, not owned
        if (geometry != nullptr) {
            auto geometry_type = geometry->getGeometryType();
            if (geometry_type != wkbPolygon && geometry_type != wkbMultiPolygon) {
                throw std::runtime_error("Geometry wasn't a polygon.");
            }
            OGRMultiPolygon* multi_polygon{nullptr};
            if (geometry_type == wkbPolygon) {
                // Directly means it doesn't clone the geometry, but adds it.
                mp_buffer.addGeometryDirectly(geometry);
                multi_polygon = &mp_buffer;
            } else {
                // Converts type but does not promote polygon to multipolygon.
                multi_polygon = geometry->toMultiPolygon();
            }
            assert_x_is_longitude(multi_polygon);
            vector<PixelData> settlement_pfpr = sparse_settlements(
                    settlement_arr, pfpr_arr, multi_polygon, settlement_geo_transform, population_cutoff
            );
            auto patch_components = CreatePatches(
                    multi_polygon, settlement_pfpr, settlement_geo_transform, population_per_patch
                    );
            std::copy(patch_components.begin(), patch_components.end(),
                    back_inserter(components)
                    );
            if (geometry_type == wkbPolygon) {
                bool do_not_delete{false};  // because it belongs to the layer.
                mp_buffer.removeGeometry(0, do_not_delete);
            }  // else no cleanup if geometry was a multipolygon.
        } else {
            cout << "geometry was null?" << endl;
        }
        ++poly_idx;
    }
    WriteVector("ug_patches.shp", components);

    return 0;
//    std::vector<Point> points;
//    gdal_raster_points(
//            std::back_inserter(points), settlement_pop_band, population_cutoff, use_subset_of_tiles
//            );
//
//    // The width of a line.
//    auto scan_length = settlement_dataset->GetRasterXSize();
//
//    std::cout << "settlements " << points.size() << std::endl;
//
//    double mosquito_meters = 200;
//    double pixel_side_meters = 30;
//    double alpha = std::pow(mosquito_meters / pixel_side_meters, 2);
//    auto pixel_sets = PixelSets(points, alpha, scan_length);
//    std::cout << "patches " << pixel_sets.size() << std::endl;
//
//    double cost_distance = 1.0;
//    PolylineComponents(pixel_sets, alpha, cost_distance, scan_length);
//    return 0;
}

int main(int argc, char* argv[])
{
    int retval{0};
    try {
        retval = entry(argc, argv);
    } catch (const std::exception& exc) {
        cout << exc.what() << endl;
        retval = 1;
    } catch (...) {
        cout << "unknown exception" << endl;
    }
    return retval;
}
