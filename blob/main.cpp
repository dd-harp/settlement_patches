#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/assertions.h>
#include <fstream>
#include <iostream>
#include <list>
#include <vector>
#include "tiffio.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment;
typedef CGAL::Alpha_shape_vertex_base_2<K> Vb;
typedef CGAL::Alpha_shape_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds> Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2> Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;


template<class OutputIterator>
void alpha_edges(const Alpha_shape_2 &A, OutputIterator out) {
    Alpha_shape_edges_iterator it = A.alpha_shape_edges_begin(),
            end = A.alpha_shape_edges_end();
    for (; it != end; ++it)
        *out++ = A.segment(*it);
}


template<class OutputIterator>
bool tiff_input(OutputIterator out, const std::string& filename, double cutoff) {
    TIFF* tif = TIFFOpen(filename.c_str(), "r");
    if (nullptr == tif) {
        return false;
    }
    uint32 imageWidth=0;
    uint32 imageLength=0;
    uint32 tileWidth=0;
    uint32 tileLength=0;

    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &imageWidth);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imageLength);
    TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileWidth);
    TIFFGetField(tif, TIFFTAG_TILELENGTH, &tileLength);
    auto buffer = static_cast<double*>(_TIFFmalloc(TIFFTileSize(tif)));
    for (uint32 y=0; y<imageLength; y+= tileLength) {
        for (uint32 x=0; x<imageWidth; x+= tileWidth) {
            // Pull one tile at a time.
            uint32 slice_index=0;  // Only have one slice.
            tsample_t sample_plane=0;  // Data aren't in separate planes.
            TIFFReadTile(tif, buffer, x, y, slice_index, sample_plane);

            // Insert points from that tile into our point set.
            for (uint32 j=0; j<tileLength; j++) {
                for (uint32 i=0; i<tileWidth; i++) {
                    if (buffer[j*tileWidth + i] > cutoff) {
                        out << Point(i + x, j+y);
                    }
                }
            }
        }
    }
    _TIFFfree(buffer);
    TIFFClose(tif);
    return true;
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
    if (!tiff_input(std::back_inserter(points), "hrsl_uga_pop.tif", 1.0))
        return -1;
    Alpha_shape_2 A(points.begin(), points.end(),
                    FT(10000),
                    Alpha_shape_2::GENERAL);
    std::vector<Segment> segments;
    alpha_edges(A, std::back_inserter(segments));
    std::cout << "Alpha Shape computed" << std::endl;
    std::cout << segments.size() << " alpha shape edges" << std::endl;
    std::cout << "Optimal alpha: " << *A.find_optimal_alpha(1) << std::endl;
    return 0;
}
