#ifndef BLOB_CONNECTED_SETTLEMENTS_H
#define BLOB_CONNECTED_SETTLEMENTS_H
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

namespace dd_harp {

/*! Is this edge or point on the boundary of the complex?
 * or is it internal? */
    template<typename ComplexElement, typename Complex>
    bool EdgePoint(const ComplexElement &point, const Complex &complex) {
        auto classification = complex.classify(point);
        return classification == Complex::SINGULAR || classification == Complex::REGULAR;
    }


/*! Make sets of pixels that are connected.
 *
 * @param complex
 * @param scan_length
 * @return
 */
    template<typename PointList>
    std::map<PixelKey, std::shared_ptr<std::vector<size_t>>>
    PixelSets(const PointList &points, double alpha, uint32 scan_length) {
        typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
        typedef Kernel::FT FT;
        typedef CGAL::Alpha_shape_vertex_base_2<Kernel> Vb;
        typedef CGAL::Alpha_shape_face_base_2<Kernel> Fb;
        typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
        typedef CGAL::Delaunay_triangulation_2<Kernel, Tds> Triangulation_2;
        typedef CGAL::Alpha_shape_2<Triangulation_2> Alpha_shape_2;
        typedef Alpha_shape_2::Vertex_handle Vertex_handle;

        Alpha_shape_2 complex(
                points.begin(), points.end(),
                FT(alpha),
                Alpha_shape_2::GENERAL
        );

        using RankMap = std::map<PixelKey, int>;
        using RankPMap = boost::associative_property_map<RankMap>;
        using ParentMap = std::map<PixelKey, PixelKey>;
        using ParentPMap = boost::associative_property_map<ParentMap>;

        RankMap rank;
        RankPMap rank_pmap(rank);
        ParentMap parent;
        ParentPMap parent_pmap(parent);

        // This finds connected components in the graph using the Union-Find algorithm.
        boost::disjoint_sets<RankPMap, ParentPMap> nearby_sets{rank_pmap, parent_pmap};
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
        std::map<PixelKey, std::shared_ptr<std::vector<size_t>>> sets;
        for (auto set_size: rank) {
            sets[set_size.first] = std::make_shared<std::vector<size_t>>(set_size.second);
        }

        // Fill in which pixels belong to which sets.
        for (auto element: parent) {
            sets[element.second]->push_back(element.first);
        }
        return sets;
    }
}

#endif //BLOB_CONNECTED_SETTLEMENTS_H
