#include <fstream>
#include <filesystem>
#include <iostream>
#include <list>
#include <algorithm>
#include <vector>
#include <map>
#include <cmath>
#include <cassert>
#include <array>
#include <string.h>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/program_options.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/assertions.h>
#include "gdal/gdal_priv.h"
#include "gdal/cpl_conv.h"
#include "geotiff/xtiffio.h"
#include "geotiff/geotiffio.h"
#include "read_tiff.h"
#include "gdal_raster.h"
#include "pixel.h"
#include "connected_settlements.h"
#include "simple_patches.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_2 Point;
typedef Kernel::Segment_2 Segment;
typedef CGAL::Alpha_shape_vertex_base_2<Kernel> Vb;
typedef CGAL::Alpha_shape_face_base_2<Kernel> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<Kernel, Tds> Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2> Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;
typedef Alpha_shape_2::Alpha_shape_vertices_iterator Alpha_shape_vertices_iterator;
typedef Alpha_shape_2::Vertex_handle Vertex_handle;

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

/*!
 *
 * Call as:
 *     std::vector <Segment> segments;
 *     write_apha_edges(alpha_shape, std::back_inserter(segments));
 *
 * @tparam OutputIterator
 * @param A
 * @param out
 */
template<class OutputIterator>
void write_apha_edges(const Alpha_shape_2 &A, OutputIterator out) {
    /*! Write each edge of an alpha shape into the output. */
    auto it = A.alpha_shape_edges_begin();
    auto end = A.alpha_shape_edges_end();
    for (; it != end; ++it, ++out) {
        *out = A.segment(*it);
    }
}



template<class OutputIterator>
bool file_input(OutputIterator out) {
    std::ifstream is("./data/fin", std::ios::in);
    if (is.fail()) {
        std::cerr << "unable to open file for input" << std::endl;
        return false;
    }
    int n;
    is >> n;
    std::cout << "Reading " << n << " points from file" << std::endl;
    CGAL::cpp11::copy_n(std::istream_iterator<Point>(is), n, out);
    return true;
}


void print_gdal_drivers() {
    for (int driver_idx=0; driver_idx < GDALGetDriverCount(); ++driver_idx) {
        auto driver = GDALGetDriver(driver_idx);
        cout << "driver " << driver_idx << " " << GDALGetDriverShortName(driver) << " "
             << GDALGetDriverLongName(driver) << endl;
    }
    if (!GDALGetDriverByName("GTiff")) {
        cout << "Cannot load GTiff driver" << endl;
    }
}


void read_geotiff(const fs::path& geotiff_filename) {
    auto parent_directory = geotiff_filename.parent_path();
    auto file_stem = geotiff_filename.stem().string();
    auto without_extension = parent_directory / file_stem;
    cout << "Opening " << geotiff_filename << endl;
    const char *allowed_drivers[3];
    allowed_drivers[0] = "GTiff";
    allowed_drivers[1] = "GeoTIFF";
    allowed_drivers[2] = nullptr;

    auto parent_iter = fs::directory_iterator(parent_directory);
    file_stem.append(".");
    cout << "stem is " << file_stem << endl;
    auto compare_len = file_stem.size();
    const auto sibling_cnt = count_if(begin(parent_iter), end(parent_iter), [&] (auto dirent) {
        return dirent.path().filename().string().compare(0, compare_len, file_stem) == 0;
    }) - 1;  // for the file itself.
    cout << "sibling count " << sibling_cnt << endl;
    char* siblings[sibling_cnt + 1];
    size_t sibling_idx = 0;
    for (auto& sibling : fs::directory_iterator(parent_directory)) {
        bool matches = sibling.path().filename().string().compare(0, compare_len, file_stem) == 0;
        bool same_file = sibling.path() == geotiff_filename;
        if (matches && !same_file) {
            auto name = sibling.path().string();
            siblings[sibling_idx] = new char[name.size() + 1];
            strcpy(siblings[sibling_idx++], &name[0]);
            cout << "sibling " << name << endl;
        }
    }
    assert(sibling_idx == sibling_cnt);
    siblings[sibling_cnt] = nullptr;

    auto dataset = static_cast<GDALDataset*>(GDALOpenEx(
            geotiff_filename.c_str(),
            GDAL_OF_RASTER | GDAL_OF_READONLY | GDAL_OF_VERBOSE_ERROR,
            allowed_drivers,
            nullptr,
            siblings
            ));

    for (size_t sib_del=0; sib_del < sibling_cnt; ++sib_del) {
        delete[] siblings[sib_del];
    }

    if (!dataset) {
        cout << "GDAL could not open " << geotiff_filename << endl;
        return;
    }

    const char* projection = dataset->GetProjectionRef();
    cout << "projection: " << projection << endl;

    const int xsize = dataset->GetRasterXSize();
    const int ysize = dataset->GetRasterYSize();
    cout << "GDAL says x,y (" << xsize << ", " << ysize << ")" << endl;

    vector<double> geoTransform(6);
    if (dataset->GetGeoTransform(&geoTransform[0]) == CE_None) {
        cout << "GeoTransform top left x " << geoTransform[0] << endl;
        cout << "GeoTransform west-east pixel resolution " << geoTransform[1] << endl;
        cout << "GeoTransform top left y " << geoTransform[3] << endl;
        cout << "GeoTransform north-south pixel resolution  " << geoTransform[5] << endl;
    }

    cout << dataset->GetRasterCount() << " raster bands" << endl;
    GDALRasterBand* band = dataset->GetRasterBand(1);
    array<int,2> block_size;
    band->GetBlockSize(&block_size[0], &block_size[1]);
    cout << "block size (" << block_size[0] << ", " << block_size[1] << ")" << endl;
    cout << "raster type " << GDALGetDataTypeName(band->GetRasterDataType()) << endl;

    array<int,2> block_cnt = {
            (band->GetXSize() + block_size[0] - 1) / block_size[0],
            (band->GetYSize() + block_size[1] - 1) / block_size[1],
    };
    cout << "block count (" << block_cnt[0] << ", " << block_cnt[1] << ")" << endl;

    GDALClose(dataset);
}


po::options_description parser()
{
    po::options_description options("blob create");
    options.add_options()
            ("help", "write help message")
            ("settlement", po::value<fs::path>(), "settlement layer GeoTIFF")
            ;
    return options;
}


// Reads a list of points and returns a list of segments
// corresponding to the Alpha shape.
int main(int argc, char* argv[]) {
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, parser()), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << parser() << endl;
        return 0;
    }
    fs::path filename{"/home/adolgert/dev/spacepop/data/hrsl_uga_pop.tif"};
    if (vm.count("settlement")) {
        filename = vm["settlement"].as<fs::path>();
    } else {
        cout << "Using default file " << filename << endl;
    }

    if (!fs::exists(filename)) {
        cout << "Could not find file " << filename << "." << endl;
        return 3;
    }
    GDALAllRegister();

    int subset = 40;  // Use only a few tiles in the corner.
    auto dataset = OpenGeoTiff(filename);
    GDALRasterBand* band = dataset->GetRasterBand(1);
    std::vector<Point> points;
    gdal_raster_points(std::back_inserter(points), band, 0.1, subset);

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
