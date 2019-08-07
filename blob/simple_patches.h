//
// Created by adolgert on 8/7/19.
//

#ifndef BLOB_SIMPLE_PATCHES_H
#define BLOB_SIMPLE_PATCHES_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/Polyline_simplification_2/Squared_distance_cost.h>
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

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_2 Point;
typedef Kernel::Segment_2 Segment;
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef CGAL::Alpha_shape_vertex_base_2<Kernel> Vb;
typedef CGAL::Alpha_shape_face_base_2<Kernel> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<Kernel, Tds> Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2> Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;
typedef Alpha_shape_2::Alpha_shape_vertices_iterator Alpha_shape_vertices_iterator;
typedef Alpha_shape_2::Vertex_handle Vertex_handle;

typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, Tds, CGAL::Exact_predicates_tag> CDT;
typedef CGAL::Constrained_triangulation_plus_2<CDT>     CT;


template<typename PIXELSETMAP>
void PolylineComponents(PIXELSETMAP pixel_set, double alpha, uint32 scan_length) {
    CT constrained_triangulation;

    for (auto& pixel_iter: pixel_set) {
        std::vector<Point> pixel_bounds;
        for (auto& pixel: *pixel_iter.second) {
            auto pt = PixelToPoint(pixel, scan_length);
            // It's a uint32 coming in. We need to use real coordinates near this line.
            double x = std::get<0>(pt);
            double y = std::get<1>(pt);
            pixel_bounds.emplace_back(Point(x, y));
            pixel_bounds.emplace_back(Point(x + 1, y));
            pixel_bounds.emplace_back(Point(x, y + 1));
            pixel_bounds.emplace_back(Point(x + 1, y + 1));
        }

        Alpha_shape_2 complex(
                pixel_bounds.begin(), pixel_bounds.end(),
                FT(alpha),
                Alpha_shape_2::GENERAL
        );

        auto edge_iter = complex.alpha_shape_edges_begin();
        auto edge_end = complex.alpha_shape_edges_end();
        for (; edge_iter != edge_end; ++edge_iter) {
            auto segment = complex.segment(*edge_iter);
            constrained_triangulation.insert_constraint(segment.source(), segment.target());
        }
    }
}

#endif //BLOB_SIMPLE_PATCHES_H
