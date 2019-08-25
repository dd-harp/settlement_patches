#ifndef BLOB_PROJECTION_H
#define BLOB_PROJECTION_H

#include "gdal/ogr_spatialref.h"

namespace spacepop {
    void reproject(OGRSpatialReference* purely_lat_long_srs);
}

#endif //BLOB_PROJECTION_H
