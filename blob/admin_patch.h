#ifndef BLOB_ADMIN_PATCH_H
#define BLOB_ADMIN_PATCH_H

#include "gdal/gdal_priv.h"


namespace spacepop
{
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
     */
    void CreatePatches(OGRMultiPolygon* admin, GDALRasterBand* settlement, GDALRasterBand* PfPR);
}

#endif //BLOB_ADMIN_PATCH_H
