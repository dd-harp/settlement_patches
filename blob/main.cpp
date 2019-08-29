
#include <filesystem>
#include <iostream>
#include <boost/program_options.hpp>

#include "gdal/gdal_priv.h"
#include "gdal/ogrsf_frmts.h"
#include "gtest/gtest.h"

#include "gdal_raster.h"
#include "projection.h"

using namespace spacepop;
using namespace std;
namespace fs = std::filesystem;
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
            ;
    for (auto const& [name, default_path] : path_argument) {
        options.add_options()(name.c_str(), po::value<fs::path>(), name.c_str());
    }
    return options;
}


/*! Given command line args, assign file paths, using defaults.
 *
 * @param vm Variables initialized from teh command line.
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


// Reads a list of points and returns a list of segments
// corresponding to the Alpha shape.
int main(int argc, char* argv[]) {
    ios_base::sync_with_stdio(false);
    map<string,fs::path> input_path = {
            {"settlement", "/home/adolgert/dev/spacepop/data/hrsl/hrsl_uga_pop.tif"},
            {"pfpr", "/home/adolgert/dev/spacepop/data/PfPR/Raster Data/PfPR_rmean/2019_Global_PfPR_2017.tif"},
            {"admin", "/home/adolgert/dev/spacepop/data/uga_admbnda/uga_admbnda_adm3_UBOS_v5.shp"},
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
    if (vm.count("test")) {
        ::testing::InitGoogleTest(&argc, argv);
        return RUN_ALL_TESTS();
    }

    // This initializes GDAL's list of drivers to read and write files.
    GDALAllRegister();
    cout << "hrsl " << input_path.at("settlement") << endl;
    auto settlement_dataset = OpenGeoTiff(input_path.at("settlement"));
    GDALRasterBand* settlement_pop_band = settlement_dataset->GetRasterBand(1);
    double settlement_geo_transform[6];
    if (settlement_dataset->GetGeoTransform(settlement_geo_transform) != CE_None) {
        cout << "Could not get settlement transform" << endl;
        return 9;
    }

    auto pfpr_dataset = OpenGeoTiff(input_path.at("pfpr"));
    GDALRasterBand* pfpr_band = pfpr_dataset->GetRasterBand(1);
    double pfpr_geo_transform[6];
    if (settlement_dataset->GetGeoTransform(pfpr_geo_transform) != CE_None) {
        cout << "Could not get pfpr transform" << endl;
        return 10;
    }

    auto projection_ref = settlement_dataset->GetProjectionRef();
    OGRSpatialReference lat_long_srs(projection_ref);
    auto to_projected = reproject(&lat_long_srs);

    auto admin_dataset = static_cast<GDALDataset *>(GDALOpenEx(
            input_path.at("admin").c_str(),
            GDAL_OF_VECTOR | GDAL_OF_READONLY | GDAL_OF_VERBOSE_ERROR,
            nullptr,
            nullptr,
            nullptr
    ));
    if (admin_dataset == nullptr) {
        cout << "Could not load dataset from " << input_path.at("admin") << endl;
        return 7;
    }
    if (admin_dataset->GetLayerCount() == 0) {
        cout << "Admin dataset has no layers" << endl;
        return 8;
    }

    // Loop over the admin layers, making patches within each layer.
    OGRLayer* first_admin_layer = *admin_dataset->GetLayers().begin();
    for (auto& admin_geometry: first_admin_layer) {
        auto geometry = admin_geometry->GetGeometryRef();
        if (geometry != nullptr) {
            auto geometry_type = geometry->getGeometryType();
            if (geometry_type == wkbPolygon || geometry_type == wkbMultiPolygon) {
                OGRMultiPolygon* multi_polygon = geometry->toMultiPolygon();

            } else {
                cout << "geometry wasn't a polygon!" << endl;
            }
        } else {
            cout << "geometry was null?" << endl;
        }
    }

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
