#include <cmath>
#include <deque>
#include <iostream>
#include <limits>
#include <memory>

#include "boost/geometry.hpp"
#include "boost/geometry/geometries/point_xy.hpp"
#include "boost/geometry/geometries/polygon.hpp"
#include "boost/geometry/geometries/segment.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/bc_clustering.hpp"
#include "boost/graph/connected_components.hpp"
#include "boost/polygon/voronoi.hpp"
#include "boost/polygon/polygon.hpp"
#include <boost/property_map/property_map.hpp>
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

bool cell_in_admin(vector<PixelData>& pixel_data, size_t idx) {
    if (idx >= pixel_data.size()) {
        return false;  // Some cells aren't associated with input points.
    }
    auto relation{pixel_data.at(idx).overlap};
    return relation == Overlap::in || relation == Overlap::on;
}

using Graph=boost::adjacency_list<boost::listS, boost::vecS>;
using GraphEdge=boost::graph_traits<Graph>::edge_descriptor;
using GraphVertex=boost::graph_traits<Graph>::vertex_descriptor;

struct centrality_done {
    int _partition_cnt;
    int _partition_idx;
    explicit centrality_done(int partition_cnt) : _partition_cnt{partition_cnt}, _partition_idx{0} {}
    centrality_done(const centrality_done&) = default;

    bool operator()(double maximum_centrality, GraphEdge edge_to_remove, const Graph& graph) {
        bool done = (this->_partition_idx >= this->_partition_cnt);
        this->_partition_idx++;
        cout << "centrality " << maximum_centrality << " "
            << edge_to_remove << " " << _partition_idx
            << " partition_cnt " << this->_partition_cnt << " " << done << endl;
        return done;
    }
};


void create_neighbor_graph(vector<PixelData>& settlement_pfpr) {
    array<double, 2> bmin{numeric_limits<double>::max(), numeric_limits<double>::max()};
    array<double, 2> bmax{numeric_limits<double>::min(), numeric_limits<double>::min()};
    for (const auto& data_p: settlement_pfpr) {
        if (data_p.overlap == Overlap::in || data_p.overlap == Overlap::on) {
            auto cpx = data_p.centroid_in.get<0>();
            bmin[0] = min(cpx, bmin[0]);
            bmax[0] = min(cpx, bmax[0]);
            auto cpy = data_p.centroid_in.get<1>();
            bmin[1] = min(cpy, bmin[1]);
            bmax[1] = min(cpy, bmax[1]);
        } else if (data_p.overlap == Overlap::out) {
            auto cpx = data_p.centroid_out.get<0>();
            bmin[0] = min(cpx, bmin[0]);
            bmax[0] = min(cpx, bmax[0]);
            auto cpy = data_p.centroid_out.get<1>();
            bmin[1] = min(cpy, bmin[1]);
            bmax[1] = min(cpy, bmax[1]);
        } else {
            throw runtime_error("data overlap unknown");
        }
    }
    const int max_val{std::numeric_limits<int>::max() - 1};
    array<double, 2> scale = { max_val / (bmax[0] - bmin[0]), max_val / (bmax[1] - bmin[1])};

    std::vector<ipoint> points;
    int settlement_cnt{0};
    for (const auto& sd: settlement_pfpr) {
        if (sd.overlap == Overlap::in || sd.overlap == Overlap::on) {
            auto cpx = sd.centroid_in.get<0>();
            auto cpy = sd.centroid_in.get<1>();
            points.emplace_back(ipoint{
                    static_cast<int>(lround((cpx - bmin[0]) * scale[0])),
                    static_cast<int>(lround((cpy - bmin[1]) * scale[1]))
            });
            ++settlement_cnt;
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
    std::vector<isegment> segments;  // There are no input segments.
    voronoi_diagram<double> vd;
    construct_voronoi(points.begin(), points.end(), segments.begin(), segments.end(), &vd);

    // Make a graph
    Graph connection{settlement_pfpr.size()};

    for (auto& set_edge_unseen: vd.edges()) {
        if (set_edge_unseen.is_primary()) {
            set_edge_unseen.color(0);
        }
    }
    // Walk edges to turn them into graph links.
    int edge_cnt{0};
    for (auto& edge_iter: vd.edges()) {
        if (edge_iter.is_primary() && edge_iter.color() == 0) {
            edge_iter.color(1);
            edge_iter.twin()->color(1);
            auto cell_a = edge_iter.cell()->source_index();
            auto cell_b = edge_iter.twin()->cell()->source_index();
            if (cell_in_admin(settlement_pfpr, cell_a) && cell_in_admin(settlement_pfpr, cell_b)) {
                auto [added_edge, added] = add_edge(cell_a, cell_b, connection);
                if (!added) {
                    cout << "not added " << cell_a << " " << cell_b << endl;
                }
                ++edge_cnt;
            }
        }
    }
    if (edge_cnt == 0 || settlement_cnt == 0) {
        return;
    }
    const int settle_size{500};
    if (settlement_cnt > settle_size) {
        using VertexIndex=map<GraphVertex, int>;
        VertexIndex vertex_index;
        auto[vert_iter, vert_end] = vertices(connection);
        for (int vert_idx = 0; vert_iter != vert_end; ++vert_iter, ++vert_idx) {
            vertex_index[*vert_iter] = vert_idx;
        }
        int split_cnt = settlement_cnt / settle_size;
        auto doneness = centrality_done{split_cnt};
        using CentralityMap=map<GraphEdge, double>;
        CentralityMap centrality;
        boost::associative_property_map<CentralityMap> pcentrality{centrality};
        boost::associative_property_map<VertexIndex> pvertex_index{vertex_index};
        cout << "centrality start " << edge_cnt << " edges "
            << split_cnt << " verts " << settlement_cnt << endl;
        boost::betweenness_centrality_clustering(connection, doneness, pcentrality, pvertex_index);
        cout << "splits " << doneness._partition_idx << " splits" << endl;
    }
}


struct MPDelete {
    void operator()(OGRMultiPolygon* p) { OGRGeometryFactory::destroyGeometry(p); }
};


void CreatePatches(
        OGRMultiPolygon* admin, vector<PixelData>& settlement_pfpr,
        const std::vector<double>& settlement_geo_transform
        )
{
    // Work in projection where units are meters.
    OGREnvelope polygon_bounding_box;
    admin->getEnvelope(&polygon_bounding_box);
    auto [project, unproject] = projection_for_lat_long(polygon_bounding_box.MinY, polygon_bounding_box.MinX);

    // Transform the incoming multipolygon in place, not a copy.
    auto local_admin = unique_ptr<OGRMultiPolygon, MPDelete>(admin->clone()->toMultiPolygon(), MPDelete());
    local_admin->transform(project.get());
    // Use Boost Polygon because it lets us create things and intersect and delaunay them.
    dmpolygon admin_bg = convert_gdal_to_boost(local_admin.get());

    split_patches_retaining_pfpr(settlement_pfpr, settlement_geo_transform, admin_bg, project);

    // Should return a graph with an index into settlement_pfpr.
    create_neighbor_graph(settlement_pfpr);

    // Cluster on the graph, excluding nodes that are outside the polygon.
}

}
