#ifndef BLOB_PROJECTION_H
#define BLOB_PROJECTION_H

#include <memory>
#include <tuple>

#include "gdal/ogr_spatialref.h"

namespace dd_harp {
    std::shared_ptr<OGRCoordinateTransformation>
            reproject(OGRSpatialReference* purely_lat_long_srs);

    std::tuple<std::shared_ptr<OGRCoordinateTransformation>,std::shared_ptr<OGRCoordinateTransformation>>
    projection_for_lat_long(double latitude, double longitude);
}

#endif //BLOB_PROJECTION_H
