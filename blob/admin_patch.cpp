#include <cmath>
#include <deque>
#include <fstream>
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
#include "component_data.h"
#include "metis_io.h"
#include "on_demand_raster.h"
#include "patch_graph.h"
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


struct centrality_done {
    int _partition_cnt;
    int _partition_idx;
    double _maximum_centrality;
    explicit centrality_done(int partition_cnt, double centrality)
    : _partition_cnt{partition_cnt}, _partition_idx{0}, _maximum_centrality{centrality} {}
    centrality_done(const centrality_done&) = default;

    bool operator()(double maximum_centrality, PatchGraphEdge edge_to_remove, const PatchGraph& graph) {
        bool steps_high = (this->_partition_idx >= this->_partition_cnt);
        bool centrality_enough = (maximum_centrality < this->_maximum_centrality);
        this->_partition_idx++;
        cout << "centrality " << maximum_centrality << " "
            << edge_to_remove << " " << _partition_idx
            << " partition_cnt " << this->_partition_cnt << " " << centrality_enough << endl;
        return steps_high || centrality_enough;
    }
};


template<typename GRAPH>
size_t component_count(const GRAPH& connection) {
    map<PatchGraphVertex, size_t> component_map;
    boost::associative_property_map<map<PatchGraphVertex, size_t>> pcomponent_map{component_map};
    boost::connected_components(connection, pcomponent_map);

    std::unordered_set<size_t> component_cnt;
    for (const auto& how_many: component_map) {
        component_cnt.insert(how_many.second);
    }
    return component_cnt.size();
}


void write_components(const vector<ComponentData>& component_data) {
    fstream csv{"components.csv", fstream::out | fstream::app};
    int component_idx{0};
    for (const auto& cd: component_data) {
        auto sep{", "};
        csv << component_idx << sep << cd.population << sep << cd.settlements << sep
            << cd.pfpr << sep << cd.centroid_lat_long[0] << sep << cd.centroid_lat_long[1] << endl;
        ++component_idx;
    }
}


PatchGraph create_neighbor_graph(vector<PixelData>& settlement_pfpr) {
    array<double, 2> bmin{numeric_limits<double>::max(), numeric_limits<double>::max()};
    array<double, 2> bmax{numeric_limits<double>::min(), numeric_limits<double>::min()};
    for (const auto &data_p: settlement_pfpr) {
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
    const int max_val{std::numeric_limits<int>::max() / 2};
    double scale = max_val / max(bmax[0] - bmin[0], bmax[1] - bmin[1]);

    std::vector<ipoint> points;
    int settlement_cnt{0};
    for (const auto &sd: settlement_pfpr) {
        if (sd.overlap == Overlap::in || sd.overlap == Overlap::on) {
            auto cpx = sd.centroid_in.get<0>();
            auto cpy = sd.centroid_in.get<1>();
            points.emplace_back(ipoint{
                    static_cast<int>(lround((cpx - bmin[0]) * scale)),
                    static_cast<int>(lround((cpy - bmin[1]) * scale))
            });
            ++settlement_cnt;
        }
        // Removing the partial points that are outside the bounds.
        if (sd.overlap == Overlap::out) {
            auto cpx = sd.centroid_out.get<0>();
            auto cpy = sd.centroid_out.get<1>();
            points.emplace_back(ipoint{
                    static_cast<int>(lround((cpx - bmin[0]) * scale)),
                    static_cast<int>(lround((cpy - bmin[1]) * scale))
            });
        }
    }
    std::vector<isegment> segments;  // There are no input segments.
    voronoi_diagram<double> vd;
    construct_voronoi(points.begin(), points.end(), segments.begin(), segments.end(), &vd);

    // Make a graph that knows which settlements are in it.
    map<size_t, size_t> settlement_idx_to_graph_idx;
    size_t remain_idx{0};
    size_t vector_idx{0};
    for (const auto &pixel_read: settlement_pfpr) {
        auto overlap = pixel_read.overlap;
        if (overlap == Overlap::in || overlap == Overlap::on) {
            settlement_idx_to_graph_idx[vector_idx] = remain_idx;
            ++remain_idx;
        }
        ++vector_idx;
    }

    // Walk edges to turn them into graph links.
    for (auto &set_edge_unseen: vd.edges()) {
        if (set_edge_unseen.is_primary()) {
            set_edge_unseen.color(0);
        }
    }
    vector<pair<int, int>> edges;
    int edge_cnt{0};
    for (auto &edge_iter: vd.edges()) {
        if (edge_iter.is_primary() && edge_iter.color() == 0) {
            edge_iter.color(1);
            edge_iter.twin()->color(1);
            auto cell_a = edge_iter.cell()->source_index();
            auto cell_b = edge_iter.twin()->cell()->source_index();
            if (cell_in_admin(settlement_pfpr, cell_a) && cell_in_admin(settlement_pfpr, cell_b)) {
                edges.emplace_back(make_pair(settlement_idx_to_graph_idx[cell_a], settlement_idx_to_graph_idx[cell_b]));
                ++edge_cnt;
            }
        }
    }
    const size_t settlement_in_admin_cnt{settlement_idx_to_graph_idx.size()};
    PatchGraph connection{
            edges.begin(),
            edges.end(),
            settlement_in_admin_cnt,
            edges.size()
    };
    for (const auto &conn_iter: settlement_idx_to_graph_idx) {
        connection[conn_iter.second].index = conn_iter.first;
    }
    return connection;
}

//void blah() {
//
//    const int settle_size{200};
//    if (settlement_cnt > settle_size) {
//        size_t component_cnt = component_count(connection);
//        cout << "graph starts with " << component_cnt << " components" << endl;
//        using VertexIndex=map<PatchGraphVertex, int>;
//        VertexIndex vertex_index;
//        auto[vert_iter, vert_end] = vertices(connection);
//        for (int vert_idx = 0; vert_iter != vert_end; ++vert_iter, ++vert_idx) {
//            vertex_index[*vert_iter] = vert_idx;
//        }
//        int split_cnt = 50;
//        auto doneness = centrality_done{split_cnt, std::pow(settle_size, 2)};
//        using CentralityMap=map<PatchGraphEdge, double>;
//        CentralityMap centrality;
//        boost::associative_property_map<CentralityMap> pcentrality{centrality};
//        boost::associative_property_map<VertexIndex> pvertex_index{vertex_index};
//        cout << "centrality start " << edge_cnt << " edges "
//             << split_cnt << " verts " << settlement_cnt << endl;
//        boost::betweenness_centrality_clustering(connection, doneness, pcentrality, pvertex_index);
//        cout << "splits " << doneness._partition_idx << " splits" << endl;
//    }
//}


vector<vector<size_t>>
connected_connections(PatchGraph& connection) {
    map<PatchGraphVertex, size_t> component_map;
    boost::associative_property_map<map<PatchGraphVertex, size_t>> pcomponent_map{component_map};
    boost::connected_components(connection, pcomponent_map);

    size_t component_cnt{0};
    for (const auto &how_many: component_map) {
        component_cnt = std::max(component_cnt, how_many.second + 1);
    }

    vector<vector<size_t>> components{component_cnt};
    for (const auto &settle: component_map) {
        components.at(settle.second).push_back(settle.first);
    }
    return components;
}


vector<ComponentData>
properties_of_components(
        const vector<vector<size_t>>& components,
        vector<PixelData>& settlement_pfpr
        )
{
    int component_idx{0};
    vector<ComponentData> component_data{components.size()};
    for (const auto& settlements: components) {
        double pfpr{0};
        double pop{0};
        array<double, 2> xy{0, 0};
        for (auto settle_idx: settlements) {
            const PixelData& pd{settlement_pfpr.at(settle_idx)};
            // Because pop is for the whole settlement, including the part outside.
            double add_pop = pd.area_in * pd.pop / (pd.area_out + pd.area_in);
            pop += add_pop;
            pfpr += pd.pfpr * add_pop;
            xy[0] += bg::get<0>(pd.centroid_in) * add_pop;
            xy[1] += bg::get<1>(pd.centroid_in) * add_pop;
        }
        xy[0] /= pop;
        xy[1] /= pop;
        component_data.at(component_idx).population = pop;
        component_data.at(component_idx).pfpr = pfpr / pop;
        component_data.at(component_idx).settlements = settlements.size();
        component_data.at(component_idx).centroid_projected = xy;
        component_data.at(component_idx).centroid_lat_long = xy;  // This will be transformed.

        ++component_idx;
    }
    return component_data;
}


struct MPDelete {
    void operator()(OGRMultiPolygon* p) { OGRGeometryFactory::destroyGeometry(p); }
};


vector<ComponentData>
CreatePatches(
        OGRMultiPolygon* admin, vector<PixelData>& settlement_pfpr,
        const std::vector<double>& settlement_geo_transform, double population_per_patch
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
    auto graph = create_neighbor_graph(settlement_pfpr);
    // Cluster on the graph, excluding nodes that are outside the polygon.
    auto grouped = split_with_metis(graph, settlement_pfpr, population_per_patch);

    auto component_data = properties_of_components(grouped, settlement_pfpr);
    for (auto& pixel_data_transform: component_data) {
        unproject->Transform(1, &pixel_data_transform.centroid_lat_long[0], &pixel_data_transform.centroid_lat_long[1]);
    }
    write_components(component_data);
    return component_data;
}

}
