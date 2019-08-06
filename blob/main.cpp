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

//! Use a single integer to identify an (x, y) pixel.
using PixelKey = size_t;

//! Convert an (x, y) into a single integer.
template<typename Point>
PixelKey PointToPixel(const Point& p, uint32 scan_length) {
    auto x = p.x();
    auto y = p.y();
    return std::lround(y) * scan_length + std::lround(x);
}


/*! Is this edge or point on the boundary of the complex?
 * or is it internal? */
template<typename ComplexElement, typename Complex>
bool EdgePoint(const ComplexElement& point, const Complex& complex) {
    auto classification = complex.classify(point);
    return classification == Alpha_shape_2::SINGULAR || classification == Alpha_shape_2::REGULAR;
}


/*! Make sets of pixels that are connected.
 *
 * @param complex
 * @param scan_length
 * @return
 */
std::map<PixelKey,std::shared_ptr<std::vector<size_t>>>
PixelSets(const Alpha_shape_2& complex, uint32 scan_length) {
    using RankMap = std::map<PixelKey, int>;
    using RankPMap = boost::associative_property_map<RankMap>;
    using ParentMap = std::map<PixelKey,PixelKey>;
    using ParentPMap = boost::associative_property_map<ParentMap>;

    RankMap rank;
    RankPMap rank_pmap(rank);
    ParentMap parent;
    ParentPMap parent_pmap(parent);

    // This finds connected components in the graph using the Union-Find algorithm.
    boost::disjoint_sets<RankPMap,ParentPMap> nearby_sets{rank_pmap, parent_pmap};
    std::vector<PixelKey> pixels;

    // First, each point is its own set.
    auto vertex = complex.alpha_shape_vertices_begin();
    auto vertex_end = complex.alpha_shape_vertices_end();
    for (; vertex != vertex_end; ++vertex) {
        // Exclude interior points from disjoint sets. Maybe a little less data, but no more correct.
        Vertex_handle vh = *vertex;
        if (EdgePoint(*vertex, complex)) {
            auto new_pixel = PointToPixel(vh->point(), scan_length);
            nearby_sets.make_set(new_pixel);
            pixels.push_back(new_pixel);
        }
    }

    // Then we union those sets if they share an edge.
    auto edge_iter = complex.alpha_shape_edges_begin();
    auto edge_end = complex.alpha_shape_edges_end();
    for (; edge_iter != edge_end; ++edge_iter) {
        if (EdgePoint(*edge_iter, complex)) {
            auto segment = complex.segment(*edge_iter);
            nearby_sets.union_set(
                    PointToPixel(segment.source(), scan_length),
                    PointToPixel(segment.target(), scan_length)
                    );
        }
    }

    // Clean up the data structure so that each parent is the
    // representative of the set, so all pixels in the same group
    // have the same parent.
    nearby_sets.compress_sets(pixels.begin(), pixels.end());

    // Create a data structure to hold all of the sets.
    std::map<PixelKey,std::shared_ptr<std::vector<size_t>>> sets;
    for (auto set_size: rank) {
        sets[set_size.first] = std::make_shared<std::vector<size_t>>(set_size.second);
    }

    // Fill in which pixels belong to which sets.
    for (auto element: parent) {
        sets[element.second]->push_back(element.first);
    }
    return sets;
}


std::vector<uint32> ImageSizes(TIFF* tif) {
    uint32 imageWidth=0;
    uint32 imageLength=0;
    uint32 tileWidth=0;
    uint32 tileLength=0;

    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &imageWidth);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imageLength);
    TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileWidth);
    TIFFGetField(tif, TIFFTAG_TILELENGTH, &tileLength);
    return {imageWidth, imageLength, tileWidth, tileLength};
}


template<class PointOutputIterator>
void tiff_input(PointOutputIterator out, TIFF* tif, double cutoff) {
    /*! Loop through a TIFF file, tile by tile, and insert values into output. */
    auto sizes = ImageSizes(tif);
    auto buffer = static_cast<double*>(_TIFFmalloc(TIFFTileSize(tif)));
    for (uint32 y=0; y<sizes[1]; y+= sizes[3]) {
        for (uint32 x=0; x<sizes[0]; x+= sizes[2]) {
            // Pull one tile at a time.
            uint32 slice_index=0;  // Only have one slice.
            tsample_t sample_plane=0;  // Data aren't in separate planes.
            TIFFReadTile(tif, buffer, x, y, slice_index, sample_plane);

            // Insert points from that tile into our point set.
            for (uint32 j=0; j<sizes[3]; j++) {
                for (uint32 i=0; i<sizes[2]; i++) {
                    if (buffer[j*sizes[2] + i] > cutoff) {
                        *out++ = Point(i + x, j+y);
                    }
                }
            }
        }
    }
    _TIFFfree(buffer);
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
    std::list<Point> points;
    std::string filename{"hrsl_uga_pop.tif"};

    TIFF* tif = TIFFOpen(filename.c_str(), "r");
    if (nullptr == tif) {
        return 7;
    }
    // The width of a line.
    auto scan_length = ImageSizes(tif)[0];
    tiff_input(std::back_inserter(points), tif, 1.0);
    TIFFClose(tif);

    Alpha_shape_2 alpha_shape(
            points.begin(), points.end(),
            FT(10000),
            Alpha_shape_2::GENERAL
    );
    auto pixel_sets = PixelSets(alpha_shape, scan_length);
    return 0;
}
