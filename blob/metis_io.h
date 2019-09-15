#ifndef BLOB_METIS_IO_H
#define BLOB_METIS_IO_H

#include <vector>

#include "patch_graph.h"
#include "sparse_settlements.h"


namespace dd_harp
{
    /*! Uses METIS program to split settlements into separate patches.
     *  This balances the number of people in each group and minimizes
     *  connection among the groups.
     *
     * @param graph Embedded planar graph of settlements next to other settlements.
     *              This ensures we don't connect people across other admin units.
     * @param settlement_pfpr Information on each patch.
     * @param population_per_patch Desired number of people per patch.
     * @return The outer vector is over patches. Each inner vector lists which
     *         settlments are in that patch.
     */
    std::vector<std::vector<size_t>> split_with_metis(
            const PatchGraph& graph, std::vector<PixelData>& settlement_pfpr, double population_per_patch
            );
}

#endif //BLOB_METIS_IO_H
