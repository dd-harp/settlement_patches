#ifndef BLOB_PROJECTION_H
#define BLOB_PROJECTION_H

#include <memory>
#include <tuple>

#include "gdal/ogr_spatialref.h"

namespace dd_harp {
    std::shared_ptr<OGRCoordinateTransformation>
            reproject(OGRSpatialReference* purely_lat_long_srs);

    /*! Creates projections into and out of Universal transverse mercator (UTM).
     *
     *  The central meridian of each zone is 500,000 meters.
     *
     * @param latitude
     * @param longitude
     * @return
     */
    std::tuple<std::shared_ptr<OGRCoordinateTransformation>,std::shared_ptr<OGRCoordinateTransformation>>
    projection_for_lat_long(double latitude, double longitude);
}

#endif //BLOB_PROJECTION_H
