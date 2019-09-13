/*! Demonstrate CGAL use b/c it's difficult to understand.
 *  CGAL documentation describes Concepts but doesn't tell you what
 *  functions you can call and what types they return, so this works
 *  through from the basic to the more complicated, working with
 *  the core types.
 *
 *  This tutorial seems on point.
 *  http://graphics.stanford.edu/courses/cs368-00-spring/manuals/CGAL_Tutorial.pdf
 */
#include <algorithm>
#include <iostream>
#include <iterator>

#include "CGAL/Homogeneous.h"
#include "CGAL/Point_2.h"

#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "CGAL/Kernel/global_functions.h"

#include "CGAL/Polygon_2.h"

#include <CGAL/Triangulation_2.h>

#include "gtest/gtest.h"


using namespace std;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// The Point class itself.
// https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Point__2.html
TEST(CGALKernel, PointCreation) {
    using Point_2 = Kernel::Point_2;
    Point_2 p{5, 9};
    Point_2 q{5.0, 9.0};
    Point_2 r{1, 1};
    ASSERT_EQ(p, q);
    using Segment_2 = Kernel::Segment_2;
    Segment_2 segment(q, r);
    cout << "Segment is " << segment << endl;
    using Vector_2 = Kernel::Vector_2;
    Vector_2 v = p - q;
    cout << "Vector " << v << endl;
}


TEST(CGALKernel, PointMath) {
    using Point_2 = Kernel::Point_2;
    using Vector_2 = Kernel::Vector_2;
    Point_2 p{5, 9};
    Point_2 q{5.0, 9.0};
    Point_2 r{1, 1};
    ASSERT_TRUE(r < q);
    Vector_2 v = p - q;
    Point_2 s = q + v;
    ASSERT_EQ(s, p);
    ASSERT_EQ(p.x(), 5.0);
    ASSERT_EQ(p.y(), 9.0);
    ASSERT_EQ(p.cartesian(0), 5.0);
    ASSERT_EQ(p.cartesian(1), 9.0);
    ASSERT_EQ(p.dimension(), 2);
    // Bbox is not part of the Cartesian kernel
    // /usr/include/CGAL/Cartesian/Cartesian_base.h
    auto bbox = p.bbox();
    cout << "bbox " << bbox << endl;
    ASSERT_EQ(bbox.xmin(), p.x());
    ASSERT_EQ(bbox.ymax(), p.y());
}


// Segment described here:
// https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Segment__2.html
TEST(CGALKernel, SegmentParts) {
    using Point_2 = Kernel::Point_2;
    using Segment_2 = Kernel::Segment_2;

    // Don't have to declare the points.
    Segment_2 seg({1, 0}, {3, 5});
    ASSERT_EQ(seg.source().x(), 1.0);
    ASSERT_EQ(seg.vertex(0).y(), 0.0);
    ASSERT_EQ(seg.point(1).x(), 3.0);

    Segment_2 flat({3, 7}, {10, 7});
    ASSERT_TRUE(flat.is_horizontal());
    ASSERT_TRUE(flat.has_on({5, 7}));
}


// There are functions in the kernel that act on these objects.
// https://doc.cgal.org/latest/Kernel_23/group__kernel__global__function.html
TEST(CGALKernel, GlobalFunctions) {
    using Point_2 = Kernel::Point_2;
    auto midpoint = CGAL::midpoint(Point_2{0, 0}, Point_2{3, 0});
    ASSERT_EQ(midpoint.x(), 1.5);

    ASSERT_TRUE(CGAL::collinear(Point_2{0, 1}, Point_2(0, 3), Point_2{0, 7}));
}


// bounded on circle... how can you figure out what this does?
// This page doesn't tell you what the return value is.
// https://doc.cgal.org/latest/Kernel_23/group__side__of__bounded__circle__grp.html#ga1a28c423cc75775dd0799f80700506c0
// Some luck with this find command:
// find . -type f -exec grep -Hi side_of_bounded_circle {} \; | grep -v _3 | grep -v typedef
// Which leads me to a return type of BoundedSide.
// find /usr/include/CGAL -type f -exec grep -H Bounded_side | grep enum
TEST(CGALKernel,BoundedOnCircle) {
    using Point = Kernel::Point_2;
    Point p1(0, -5), p2(3, -4), p3(4, 3), in(-1, 4), out(5, -1), on(0, 5);
    CGAL::Bounded_side inside, onside, outside;
    inside = CGAL::side_of_bounded_circle(p1, p2, p3, in);
    outside = CGAL::side_of_bounded_circle(p1, p2, p3, out);
    onside = CGAL::side_of_bounded_circle(p1, p2, p3, on);
    ASSERT_EQ(inside,  CGAL::ON_BOUNDED_SIDE);
    ASSERT_EQ(outside, CGAL::ON_UNBOUNDED_SIDE);
    ASSERT_EQ(onside, CGAL::ON_BOUNDARY);
}


// Homogeneous coordinates.
TEST(CGALHomogeneous,PointInLine) {
    using Point = CGAL::Point_2<CGAL::Homogeneous<long>>;
    Point p1{0, 0}, p2{3, 17, 10}, p3{9, 51, 10};
    ASSERT_EQ(CGAL::orientation(p1, p2, p3), CGAL::COLLINEAR);
}


// Polygons
// I can't find out what the iterators point to from this:
// https://doc.cgal.org/latest/Polygon/classCGAL_1_1Polygon__2.html
// It won't tell me the value type or anything. Kinda dislike y'all right now.
TEST(CGALPolygon,Construct) {
    using Point_2 = Kernel::Point_2;
    using Polygon_2 = CGAL::Polygon_2<Kernel>;

    Point_2 points[] = { {0,0}, {5.1,0}, {1,1}, {0.5,6}};
    Polygon_2 pgn(std::begin(points), std::end(points));
    // ASSERT_TRUE(pgn.is_simple());
    // ASSERT_FALSE(pgn.is_convex());

    double area{0};
    // The last argumen is a Polygon_2_traits Concept. That turns out to be
    // fulfilled by a Kernel object. How am I to know that?
    CGAL::area_2(std::begin(points), std::end(points), area, Kernel());
    cout << "area of polygon " << area << endl;

    // Returns a Vertex_const_iterator. What's that?
    auto vertex = *pgn.vertices_begin();
    // Turns out a vertex here is the Point type.
    ASSERT_EQ(vertex.x(), 0);
    ASSERT_EQ(vertex.y(), 0);
    // And the edge is a Segment type.
    auto edge = *pgn.edges_begin();
    ASSERT_EQ(edge.target().x(), 5.1);

    // circulators aren't too fancy. The end is the beginning.
    auto edge_iter = pgn.edges_circulator();
    auto edge_end = pgn.edges_circulator();
    for (; edge_iter != edge_end; ++edge_iter) {
        double x = edge_iter->source().x();
        double y = edge_iter->source().y();
    }
}


TEST(CGALTriangulation,CreateTriangulation) {
    using Point = Kernel::Point_2;
    using Triangulation = CGAL::Triangulation_2<Kernel>;
    Point points[] = {{0.4, 1}, {1, 0.3}, {0.0, -0.9}, {-1, 0},
                      {0, 0}, {-1, 1}, {1, 0.9}, {1.4, -0.3}, {0.6, 0}};
    Triangulation tr;
    tr.insert(std::begin(points), std::end(points));

    // The type of vertex is Triangulation_vertex_base_2, and it contains a
    // Triangulation_ds_vertex_base_2 in the templates. A find command gets this:
    // /usr/include/CGAL/Triangulation_vertex_base_2.h
    auto vertex = *tr.vertices_begin();
    cout << "triangulation point " << vertex.point() << endl;
    cout << "triangulation size " << std::distance(tr.vertices_begin(), tr.vertices_end()) << endl;
    auto vertex_handle = vertex.handle();

    // These are what it claims are the neighboring points.
    std::set<Point> by_point;
    auto vertex_circulator = vertex.incident_vertices();
    auto vertex_end{vertex_circulator};
    for (; vertex_circulator!=vertex_end; ++vertex_circulator) {
        by_point.emplace(vertex_circulator->point());
    }
    auto face_circulator = vertex.incident_faces();
    auto edge_circulator = vertex.incident_edges();
    auto face_handle_has_this_vertex = vertex.face();
    auto edge_end{edge_circulator};
    // find /usr/include/CGAL -type f -name "*.h" -exec grep -H Triangulation_ds_edge_circulator_2 {} \;
    for (; edge_circulator != edge_end; ++edge_circulator) {
        // How do I know that this is a pair? Some tutorial I finally found.
        // Again, finding things is hard with this library.
        auto edge_face_handle = edge_circulator->first;
        auto edge_index = edge_circulator->second;
        ; //if (edge_circulator->segment())
    }

    // We traverse faces to find the neighboring points.
    std::set<Point> neighbor_points;
    auto face_end{face_circulator};
    for (; face_circulator != face_end; ++face_circulator) {
        // Ask the face what our index is, and then ask it for the vertex
        // that comes after our index.
        auto next_vertex = face_circulator->vertex(face_circulator->ccw(face_circulator->index(vertex_handle)));
        neighbor_points.emplace(next_vertex->point());
    }

    ASSERT_EQ(neighbor_points, by_point);
}
