#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <boost/geometry.hpp>

#include "metis_io.h"

using namespace std;


namespace dd_harp
{
    vector<vector<size_t>>
    split_with_metis(const PatchGraph& graph, std::vector<PixelData>& settlement_pfpr, double population_per_patch)
    {
        const auto graph_name{"graph.metis"s};
        fstream graph_file{graph_name, fstream::out};
        fstream settle_file{"admin_metis.txt", fstream::out};
        const auto sep{" "s};
        const auto vertex_and_edge_weights{"011"s};
        const int one_based{1};
        auto vert_cnt = num_vertices(graph);
        auto edge_cnt = num_edges(graph);
        graph_file << vert_cnt << sep << edge_cnt << sep << vertex_and_edge_weights << endl;

        int settle_idx{1};
        double total_population{0};
        auto [vert, vert_end] = vertices(graph);
        for (; vert != vert_end; ++vert) {
            const PixelData& pd{settlement_pfpr[graph[*vert].index]};
            settle_file << settle_idx << sep << pd.x[0] << sep << pd.x[1] << sep
                << pd.pop << sep << pd.pfpr << sep
                << pd.centroid_in.get<0>() << sep << pd.centroid_in.get<1>() << endl;
            ++settle_idx;
            total_population += pd.pop;

            graph_file << max(1l, lround(pd.pop)) << sep;
            auto [out, out_end] = out_edges(*vert, graph);
            for (; out != out_end; ++out) {
                auto target_vertex = target(*out, graph);
                graph_file << (target_vertex + one_based) << sep;
                double dx = boost::geometry::distance(pd.centroid_in, settlement_pfpr[graph[target_vertex].index].centroid_in);
                const double minimum_distance{30};
                const int maximum_weight{10};
                graph_file << std::max(1l, lround(pow(dx / minimum_distance, -2.0) * maximum_weight)) << sep;
            }
            graph_file << endl;
        }

        size_t partition_cnt = max(1l, lround(total_population / population_per_patch));

        stringstream cmd;
        cmd << "gpmetis " << graph_name << " " << partition_cnt;
        cout << "running " << cmd.str() << endl;
        system(cmd.str().c_str());

        vector<vector<size_t>> groups{partition_cnt};

        stringstream metis_out;
        metis_out << "graph.metis.part." << partition_cnt;
        fstream in_file(metis_out.str(), fstream::in);
        auto [in_vert, in_vert_end] = vertices(graph);
        int settle_input_idx{1};
        for (; in_vert != in_vert_end; ++in_vert) {
            size_t part_idx{0};
            in_file >> part_idx;
            groups.at(part_idx).push_back(graph[*in_vert].index);
            ++settle_input_idx;
        }
        return groups;
    }
}
