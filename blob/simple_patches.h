#ifndef BLOB_SIMPLE_PATCHES_H
#define BLOB_SIMPLE_PATCHES_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <vector>

#include <boost/pending/disjoint_sets.hpp>
#include <boost/property_map/property_map.hpp>

#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "CGAL/Alpha_shape_2.h"
#include "CGAL/Alpha_shape_vertex_base_2.h"
#include "CGAL/Alpha_shape_face_base_2.h"
#include "CGAL/Delaunay_triangulation_2.h"
#include "CGAL/Constrained_Delaunay_triangulation_2.h"
#include "CGAL/Constrained_triangulation_plus_2.h"
#include "CGAL/Polyline_simplification_2/simplify.h"
#include "CGAL/Polyline_simplification_2/Squared_distance_cost.h"
#include "CGAL/algorithm.h"
#include "CGAL/assertions.h"

#include "tiffio.h"
#include "read_tiff.h"
#include "pixel.h"

namespace dd_harp {
    typedef CGAL::Exact_predicates_inexact_constructions_kernel PolygonKernel;
    typedef CGAL::Polygon_2<PolygonKernel> Polygon_2;

    template<typename PIXELS>
    Polygon_2 pixels_to_polygon(PIXELS pixels, int scan_length, double alpha) {
        typedef PolygonKernel::Point_2 Point;
        typedef PolygonKernel::Segment_2 Segment;
        // Set up the alpha shapes
        typedef CGAL::Alpha_shape_vertex_base_2<PolygonKernel> Vb;
        typedef CGAL::Alpha_shape_face_base_2<PolygonKernel> Fb;
        typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
        typedef CGAL::Delaunay_triangulation_2<PolygonKernel, Tds> Triangulation_2;
        typedef CGAL::Alpha_shape_2<Triangulation_2> Alpha_shape_2;
        typedef Alpha_shape_2::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;
        typedef Alpha_shape_2::Alpha_shape_vertices_iterator Alpha_shape_vertices_iterator;
        typedef Alpha_shape_2::Vertex_handle Vertex_handle;

        std::vector<Point> pixel_bounds;
        for (auto &pixel: pixels) {
            auto [x, y] = PixelToPoint(pixel, scan_length);
            // It's a uint32 coming in. We need to use real coordinates near this line.
            std::cout << "pixel (" << x << ", " << y << ")" << std::endl;
            pixel_bounds.emplace_back(Point(x, y));
            pixel_bounds.emplace_back(Point(x + 1, y));
            pixel_bounds.emplace_back(Point(x, y + 1));
            pixel_bounds.emplace_back(Point(x + 1, y + 1));
        }
        // Make an alpha-shape that's all of those pixel corners.
        Alpha_shape_2 complex(
                pixel_bounds.begin(), pixel_bounds.end(),
                PolygonKernel::FT(alpha),
                Alpha_shape_2::GENERAL
        );
        std::set<Point> points;
        auto all_vertex = complex.alpha_shape_vertices_begin();
        auto all_vertex_end = complex.alpha_shape_vertices_end();
        for (; all_vertex != all_vertex_end; ++all_vertex) {
            points.insert((*all_vertex)->point());
        }
        std::cout << "found " << points.size() << " points" << std::endl;

        std::multimap<Point,Point> point_to_segment;
        auto all_edge = complex.alpha_shape_edges_begin();
        auto all_edge_end = complex.alpha_shape_edges_end();
        for (; all_edge != all_edge_end; ++all_edge) {
            if (complex.classify(*all_edge) == Alpha_shape_2::REGULAR) {
                auto segment = complex.segment(*all_edge);
                point_to_segment.insert(std::make_pair(segment.source(), segment.target()));
                point_to_segment.insert(std::make_pair(segment.target(), segment.source()));
            }
        }
        std::cout << "found " << point_to_segment.size() << " segments" << std::endl;

        std::vector<Polygon_2> polygon;
        while (!points.empty()) {
            Point seed_point;
            bool finding_seed{true};
            while (finding_seed && !points.empty()) {
                auto seed_iter = points.begin();
                seed_point = *seed_iter;
                points.erase(seed_iter);
                auto search_point = point_to_segment.equal_range(seed_point);
                auto neighboring_points = std::distance(search_point.first, search_point.second);
                if (neighboring_points == 2) {
                    finding_seed = false;
                } else {
                    std::cout << "search point " << seed_point << " distance " << neighboring_points << std::endl;
                }
            }

            bool looking_for_seed{true};
            auto current_point = seed_point;
            std::vector<Point> polypoints;
            while (!points.empty() && looking_for_seed) {
                polypoints.emplace_back(current_point);
                auto search_segment = point_to_segment.equal_range(current_point);
                auto search_point = search_segment.first;
                auto next_point{current_point};
                for (; search_point != search_segment.second; ++search_point) {
                    auto found_point = search_point->second;
                    if (polypoints.size() > 1) {
                        if (polypoints[polypoints.size() - 2] != found_point) {
                            next_point = found_point;
                        }
                    } else {
                        next_point = found_point;
                    }
                }
                if (next_point != current_point) {
                    looking_for_seed = false;
                    std::cout << "didn't find a next point " << current_point << std::endl;
                } else if (current_point == seed_point) {
                    std::cout << "found seed " << seed_point << std::endl;
                    looking_for_seed = false;
                }
            }
            polygon.emplace_back(Polygon_2{polypoints.begin(), polypoints.end()});
        }
        if (polygon.empty()) {
            polygon.emplace_back(Polygon_2{});
        } else if (polygon.size() > 1) {
            std::cout << "Polygon has " << polygon.size() << " parts" << std::endl;
        }
        return polygon[0];
    }

    template<typename PIXELSETMAP>
    void PolylineComponents(PIXELSETMAP pixel_set, double alpha, double cost_distance, uint32 scan_length) {

        namespace PS = CGAL::Polyline_simplification_2;

        typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
        typedef Kernel::FT FT;
        typedef Kernel::Point_2 Point;
        typedef Kernel::Segment_2 Segment;

        // Set up the polyline simplification.
        typedef PS::Vertex_base_2<Kernel> PSVb;
        typedef CGAL::Constrained_triangulation_face_base_2<Kernel> PSFb;
        typedef CGAL::Triangulation_data_structure_2<PSVb, PSFb> PSTDS;
        typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, PSTDS, CGAL::Exact_predicates_tag> CDT;
        typedef CGAL::Constrained_triangulation_plus_2<CDT> CT;
        typedef PS::Stop_below_count_ratio_threshold Stop;
        typedef PS::Squared_distance_cost Cost;

        CT constrained_triangulation;

        for (auto &pixel_iter: pixel_set) {
            // Create a list that has every corner of every pixel in a pixel set.

//            auto vertex = *complex.alpha_shape_vertices_begin();
//            auto point = vertex->point();

            // Put the alpha-complex of that shape into the list of shapes.
//            auto edge_iter = complex.alpha_shape_edges_begin();
//            auto edge_end = complex.alpha_shape_edges_end();
//            for (; edge_iter != edge_end; ++edge_iter) {
//                auto segment = complex.segment(*edge_iter);
//                constrained_triangulation.insert_constraint(segment.source(), segment.target());
//            }
        }
        std::cout << "smoothing" << std::endl;
        PS::simplify(constrained_triangulation, Cost(), Stop(std::pow(cost_distance, 2)));
    }
}

#endif //BLOB_SIMPLE_PATCHES_H
