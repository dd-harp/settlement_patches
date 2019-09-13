#ifndef BLOB_PATCH_GRAPH_H
#define BLOB_PATCH_GRAPH_H

#include "boost/graph/adjacency_list.hpp"


namespace dd_harp {
    struct PatchDescription {
        size_t index;  // Tells us which settlement this is.
    };

    using PatchGraph=boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, PatchDescription>;
    using PatchGraphEdge=boost::graph_traits<PatchGraph>::edge_descriptor;
    using PatchGraphVertex=boost::graph_traits<PatchGraph>::vertex_descriptor;

}


#endif //BLOB_PATCH_GRAPH_H
