#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/assertions.h>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/property_map/property_map.hpp>
#include <fstream>
#include <iostream>
#include <list>
#include <vector>
#include <map>
#include <cmath>
#include "tiffio.h"
#include "read_tiff.h"
#include "pixel.h"
#include "connected_settlements.h"

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


// Reads a list of points and returns a list of segments
// corresponding to the Alpha shape.
int main() {
    std::string filename{"hrsl_uga_pop.tif"};
    int subset = 40;  // Use only a few tiles in the corner.

    TIFF* tif = TIFFOpen(filename.c_str(), "r");
    if (nullptr == tif) {
        return 7;
    }
    // The width of a line.
    auto scan_length = ImageSizes(tif)[0];
    std::vector<Point> points;
    tiff_input(std::back_inserter(points), tif, 0.1, subset);
    TIFFClose(tif);
    std::cout << "settlements " << points.size() << std::endl;

    double mosquito_meters = 200;
    double pixel_side_meters = 30;
    double alpha = std::pow(mosquito_meters / pixel_side_meters, 2);
    auto pixel_sets = PixelSets(points, alpha, scan_length);
    std::cout << "patches " << pixel_sets.size() << std::endl;
    return 0;
}
