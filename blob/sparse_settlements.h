#ifndef BLOB_SPARSE_SETTLEMENTS_H
#define BLOB_SPARSE_SETTLEMENTS_H

#include <array>
#include <vector>
#include "boost/geometry/geometries/point_xy.hpp"

class OGRMultiPolygon;

namespace dd_harp {

    class OnDemandRaster;

/*! Finds the minimum and maximum of a polygon within a particular projection transform.
 *
 * @param admin A OGRGeometry.
 * @param settlement_geo_transform Describes how coordinates relate to pixels for GDAL projections.
 * @return
 */
    std::pair <std::array<int, 2>, std::array<int, 2>>
    gdal_geometry_min_max(const OGRMultiPolygon *admin, const std::vector<double> &settlement_geo_transform);

    enum class Overlap {
        unknown, in, out, on
    };

    struct PixelData {
        std::array<int, 2> x;
        double pfpr;
        double pop;
        Overlap overlap;
        boost::geometry::model::d2::point_xy<double> centroid_in;
        boost::geometry::model::d2::point_xy<double> centroid_out;
        double area_in;
        double area_out;
        PixelData(std::array<int, 2> x, double pfpr, double pop);
    };


    std::vector<PixelData>
    sparse_settlements(
            OnDemandRaster &settlement_arr,
            OnDemandRaster &pfpr_arr,
            const OGRMultiPolygon *admin,
            const std::vector<double> &settlement_geo_transform,
            double cutoff
    );
}

#endif //BLOB_SPARSE_SETTLEMENTS_H
