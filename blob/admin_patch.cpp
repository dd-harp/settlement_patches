#include <cmath>
#include <deque>
#include <iostream>

#include "boost/geometry.hpp"
#include "boost/geometry/geometries/point_xy.hpp"
#include "boost/geometry/geometries/polygon.hpp"
#include "boost/geometry/geometries/segment.hpp"
#include "boost/polygon/voronoi.hpp"
#include "boost/polygon/polygon.hpp"
#include "gdal/ogr_geometry.h"

#include "admin_patch.h"
#include "gdal_raster.h"
#include "on_demand_raster.h"
#include "projection.h"
#include "split_patches.h"
#include "sparse_settlements.h"

using namespace boost::geometry;
using namespace std;

namespace bg = boost::geometry;
using dpoint = model::d2::point_xy<double>;
using dpolygon = model::polygon<dpoint>;
using dmpolygon = model::multi_polygon<dpolygon>;

using ipoint = model::d2::point_xy<int>;
using isegment = model::segment<ipoint>;

using boost::polygon::voronoi_builder;
using boost::polygon::voronoi_diagram;

namespace boost::polygon {

    template <>
    struct geometry_concept<ipoint> { typedef point_concept type; };

    template <>
    struct point_traits<ipoint> {
        typedef int coordinate_type;

        static inline coordinate_type get(const ipoint& point, orientation_2d orient) {
            return (orient == HORIZONTAL) ? bg::get<0>(point) : bg::get<1>(point);
        }
    };

    template <>
    struct geometry_concept<isegment> { typedef segment_concept type; };

    template <>
    struct segment_traits<isegment> {
        typedef int coordinate_type;
        typedef ipoint point_type;

        static inline point_type get(const isegment& segment, direction_1d dir) {
            if (dir.to_int()) {
                return {bg::get<1, 0>(segment), bg::get<1, 1>(segment)};
            } else {
                return {bg::get<0,0>(segment), bg::get<0,1>(segment)};
            }
        }
    };
}

namespace dd_harp {
// These templates convince boost::geometry to treat
// an array of two doubles as a Point type.
//namespace boost::geometry::traits {
//    template<> struct tag<std::array<double, 2>>
//    { typedef point_tag type; };
//    template<> struct coordinate_type<std::array<double, 2>>
//    { typedef double type; };
//    template<> struct coordinate_system<std::array<double, 2>>
//    { typedef cs::cartesian type; };
//    template<> struct dimension<std::array<double, 2>> : boost::mpl::int_<2> {};
//    template<>
//    struct access<std::array<double, 2>, 0> {
//        static double get(std::array<double, 2> const &p) {
//            return p[0];
//        }
//
//        static void set(std::array<double, 2> &p, double value) {
//            p[0] = value;
//        }
//    };
//
//    template<>
//    struct access<std::array<double, 2>, 1> {
//        static double get(std::array<double, 2> const &p) {
//            return p[1];
//        }
//
//        static void set(std::array<double, 2> &p, double value) {
//            p[1] = value;
//        }
//    };
//}

const int X{0}, Y{1};


void create_neighbor_graph(map<array<int, 2>,PixelData>& settlement_pfpr) {
    array<double, 2> bmin{numeric_limits<double>::max(), numeric_limits<double>::max()};
    array<double, 2> bmax{numeric_limits<double>::min(), numeric_limits<double>::min()};
    for (const auto& [bounds_p, data_p]: settlement_pfpr) {
        if (data_p.overlap == Overlap::in || data_p.overlap == Overlap::on) {
            auto cpx = data_p.centroid_in.get<0>();
            bmin[0] = min(cpx, bmin[0]);
            bmax[0] = min(cpx, bmax[0]);
            auto cpy = data_p.centroid_in.get<1>();
            bmin[1] = min(cpx, bmin[1]);
            bmax[1] = min(cpx, bmax[1]);
        }
        if (data_p.overlap == Overlap::out || data_p.overlap == Overlap::on) {
            auto cpx = data_p.centroid_out.get<0>();
            bmin[0] = min(cpx, bmin[0]);
            bmax[0] = min(cpx, bmax[0]);
            auto cpy = data_p.centroid_out.get<1>();
            bmin[1] = min(cpx, bmin[1]);
            bmax[1] = min(cpx, bmax[1]);
        }
    }
    array<double, 2> scale = { (1<<30) / (bmax[0] - bmin[0]), (1<<30) / (bmax[1] - bmin[1])};

    std::vector<ipoint> points;
    for (const auto& [pix_key, sd]: settlement_pfpr) {
        if (sd.overlap == Overlap::in || sd.overlap == Overlap::on) {
            auto cpx = sd.centroid_in.get<0>();
            auto cpy = sd.centroid_in.get<1>();
            points.emplace_back(ipoint{
                    static_cast<int>(lround((cpx - bmin[0]) * scale[0])),
                    static_cast<int>(lround((cpy - bmin[1]) * scale[1]))
            });
        }
        if (sd.overlap == Overlap::out || sd.overlap == Overlap::on) {
            auto cpx = sd.centroid_out.get<0>();
            auto cpy = sd.centroid_out.get<1>();
            points.emplace_back(ipoint{
                    static_cast<int>(lround((cpx - bmin[0]) * scale[0])),
                    static_cast<int>(lround((cpy - bmin[1]) * scale[1]))
            });
        }
    }
    std::vector<isegment> segments;
    voronoi_diagram<double> vd;
    construct_voronoi(points.begin(), points.end(), segments.begin(), segments.end(), &vd);

    // Make a graph

    // Walk edges to turn them into graph links.
    for (auto edge_iter = vd.edges().begin(); edge_iter != vd.edges().end(); ++edge_iter) {
        if (edge_iter->is_primary()) {
            ;
        }
    }
}


void CreatePatches(
        OGRMultiPolygon* admin, map<array<int, 2>,PixelData>& settlement_pfpr,
        const std::vector<double>& settlement_geo_transform
        )
{
    // Work in projection where units are meters.
    OGREnvelope polygon_bounding_box;
    admin->getEnvelope(&polygon_bounding_box);
    auto [project, unproject] = projection_for_lat_long(polygon_bounding_box.MinY, polygon_bounding_box.MinX);

    // Transform the incoming multipolygon in place, not a copy.
    admin->transform(project.get());
    // Use Boost Polygon because it lets us create things and intersect and delaunay them.
    dmpolygon admin_bg = convert(admin);

    split_patches_retaining_pfpr(settlement_pfpr, settlement_geo_transform, admin_bg, project);

    // Should return a graph with an index into settlement_pfpr.
    create_neighbor_graph(settlement_pfpr);

    // Cluster on the graph, excluding nodes that are outside the polygon.
}

}
