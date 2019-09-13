#ifndef BLOB_METIS_IO_H
#define BLOB_METIS_IO_H

#include <vector>

#include "patch_graph.h"
#include "sparse_settlements.h"


namespace dd_harp
{
    std::vector<std::vector<int>> write_for_metis(const PatchGraph& graph, std::vector<PixelData>& settlement_pfpr);
}

#endif //BLOB_METIS_IO_H
