#ifndef BLOB_ADMIN_PATCH_H
#define BLOB_ADMIN_PATCH_H

#include <array>
#include <map>
#include <vector>

#include "gdal/gdal_priv.h"

class OGRMultiPolygon;

namespace dd_harp
{
    class PixelData;
    /*! Given an admin unit, create patches from the settlement layer.
     *
     * @param admin The admin unit, as a multi polygon, even if it's just one polygon.
     *              This means it can have more than one exterior ring.
     *              This should be in WGS84 geographic units.
     * @param settlement A settlement band has a pixel for each settlement and a value
     *                   for estimated population in that settlement.
     *                   This should be pixels, but in Geographic units.
     * @param PfPR Prevalence of Malaria within that population, to assign to the patch.
     *             This should be pixels, but in Geographic units.
     * @param settlement_geo_transform Raster space (Pixel, Line) to projection coordinate (Xp, Yp)
     *                                 Xp = padfTransform[0] + P*padfTransform[1] + L*padfTransform[2];
     *                                 Yp = padfTransform[3] + P*padfTransform[4] + L*padfTransform[5];
     * @param pfpr_geo_transform Also a transformation from (P,L) to (Xp, Yp).
     */
    void CreatePatches(
            OGRMultiPolygon* admin, std::map<std::array<int, 2>,PixelData>& settlement_pfpr,
            const std::vector<double>& settlement_geo_transform
            );
}

#endif //BLOB_ADMIN_PATCH_H
