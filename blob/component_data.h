#ifndef BLOB_COMPONENT_DATA_H
#define BLOB_COMPONENT_DATA_H

#include <array>

namespace dd_harp {
    struct ComponentData {
        double population;
        double settlements;
        double pfpr;
        std::array<double, 2> centroid_projected;
        std::array<double, 2> centroid_lat_long;
    };
}
#endif //BLOB_COMPONENT_DATA_H
