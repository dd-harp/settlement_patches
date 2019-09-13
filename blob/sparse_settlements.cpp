#include <iostream>

#include "gdal/gdal_priv.h"
#include "gdal/ogr_geometry.h"
#include "gdal/ogr_spatialref.h"

#include "gdal_raster.h"
#include "on_demand_raster.h"
#include "sparse_settlements.h"

using namespace std;

namespace dd_harp {

    const int X{0}, Y{1};

/*! Finds the minimum and maximum of a polygon within a particular projection transform.
 *
 * @param admin A OGRGeometry.
 * @param settlement_geo_transform Describes how coordinates relate to pixels for GDAL projections.
 * @return
 */
    pair<array<int, 2>, array<int, 2>>
    gdal_geometry_min_max(const OGRMultiPolygon *admin, const std::vector<double> &settlement_geo_transform) {

        OGREnvelope polygon_bounding_box;
        admin->getEnvelope(&polygon_bounding_box);
        std::cout << "polygon bbox ((" << polygon_bounding_box.MinX << ", " << polygon_bounding_box.MaxX << "), (("
                  << polygon_bounding_box.MinY << ", " << polygon_bounding_box.MaxY << "))" << std::endl;

        std::vector<std::array<int, 2>> bounding_pixels;
        for (auto ex: {polygon_bounding_box.MinX, polygon_bounding_box.MaxX}) {
            for (auto ey: {polygon_bounding_box.MinY, polygon_bounding_box.MaxY}) {
                bounding_pixels.push_back(pixel_containing({ex, ey}, settlement_geo_transform));
            }
        }
        std::array<int, 2> settlement_min{std::numeric_limits<int>::max(), std::numeric_limits<int>::max()};
        std::array<int, 2> settlement_max{0, 0};
        for (auto bpix: bounding_pixels) {
            for (auto mcoord: {X, Y}) {
                settlement_min[mcoord] = std::min(settlement_min[mcoord], bpix[mcoord]);
                settlement_max[mcoord] = std::max(settlement_max[mcoord], bpix[mcoord]);
            }
        }
        return {settlement_min, settlement_max};
    }

    PixelData::PixelData(std::array<int, 2> x_, double pfpr_, double pop_)
    : x{x_}, pfpr{pfpr_}, pop{pop_}, overlap{Overlap::unknown}, centroid_in{}, centroid_out{},
    area_in{0}, area_out{0}
    { }


    vector<PixelData>
    sparse_settlements(
            OnDemandRaster &settlement_arr,
            OnDemandRaster &pfpr_arr,
            const OGRMultiPolygon *admin,
            const vector<double> &settlement_geo_transform,
            const double cutoff
    ) {
        auto settlement_min_max = gdal_geometry_min_max(admin, settlement_geo_transform);
        // Clip bounds because the admin units can be outside the given settlements.
        // The person supplying data will see the bounds of that data.
        const auto&[settlement_min, settlement_max] = settlement_arr.clip_to_bounds(settlement_min_max);
        vector<PixelData> settlement_pfpr;
        for (int pixel_y = settlement_min[Y]; pixel_y < settlement_max[Y]; ++pixel_y) {
            for (int pixel_x = settlement_min[X]; pixel_x < settlement_max[X]; ++pixel_x) {
                double pixel_pop = settlement_arr.at({pixel_x, pixel_y});
                if (pixel_pop > cutoff) {
                    auto settlement_coord = pixel_coord<array<double, 2>>({pixel_x, pixel_y}, settlement_geo_transform);
                    PixelData pd{{pixel_x, pixel_y}, pfpr_arr.at_coord(settlement_coord[0], settlement_coord[1]), pixel_pop};
                    settlement_pfpr.push_back(pd);
                }
            }
        }
        return settlement_pfpr;
    }

}
