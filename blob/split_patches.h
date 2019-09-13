#ifndef BLOB_SPLIT_PATCHES_H
#define BLOB_SPLIT_PATCHES_H

#include <array>
#include <map>
#include <memory>
#include <vector>

using dmpolygon = boost::geometry::model::multi_polygon<boost::geometry::model::polygon<boost::geometry::model::d2::point_xy<double>>>;
class OGRMultiPolygon;
class OGRCoordinateTransformation;


namespace dd_harp {

    class PixelData;
/*! Convert a GDAL OGRMultiPolygon into a Boost multi_polygon
 *
 * @param gdal_poly Reads from the polygon, possibly transforming it in
 *                  the process, in case it has a different winding order.
 * @return Boost multi_polygon
 */
    dmpolygon convert_gdal_to_boost(OGRMultiPolygon* gdal_poly);

    void split_patches_retaining_pfpr(
            std::vector<PixelData> &settlement_pfpr,
            const std::vector<double> &settlement_geo_transform,
            const dmpolygon &admin_bg,
            std::shared_ptr <OGRCoordinateTransformation> &project
    );
}
#endif //BLOB_SPLIT_PATCHES_H
