
#include <filesystem>
#include <iostream>
#include <list>
#include <vector>
#include <cmath>
#include <boost/program_options.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include "gdal/gdal_priv.h"
#include "gdal_raster.h"
#include "connected_settlements.h"
#include "simple_patches.h"
#include "gtest/gtest.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point;
typedef CGAL::Alpha_shape_vertex_base_2<Kernel> Vb;
typedef CGAL::Alpha_shape_face_base_2<Kernel> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<Kernel, Tds> Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2> Alpha_shape_2;

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
    map<string,fs::path> input_path = {
            {"settlement", "/home/adolgert/dev/spacepop/data/hrsl/hrsl_uga_pop.tif"},
            {"pfpr", "/home/adolgert/dev/spacepop/data/PfPR/Raster Data/PfPR_rmean/2019_Global_PfPR_2017.tif"},
            {"admin", "/home/adolgert/dev/spacepop/data/maplibrary/UGA_boundaries.shp"},
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
    auto dataset = OpenGeoTiff(input_path.at("settlement"));
    GDALRasterBand* band = dataset->GetRasterBand(1);
    std::vector<Point> points;
    gdal_raster_points(std::back_inserter(points), band, population_cutoff, use_subset_of_tiles);

    // The width of a line.
    auto scan_length = dataset->GetRasterXSize();

    std::cout << "settlements " << points.size() << std::endl;

    double mosquito_meters = 200;
    double pixel_side_meters = 30;
    double alpha = std::pow(mosquito_meters / pixel_side_meters, 2);
    auto pixel_sets = PixelSets(points, alpha, scan_length);
    std::cout << "patches " << pixel_sets.size() << std::endl;

    double cost_distance = 1.0;
    PolylineComponents(pixel_sets, alpha, cost_distance, scan_length);
    return 0;
}
