#include <iterator>

#include "boost/geometry.hpp"
#include "boost/geometry/geometries/point_xy.hpp"
#include "boost/geometry/geometries/polygon.hpp"
#include "boost/geometry/geometries/segment.hpp"
#include "boost/polygon/voronoi.hpp"
#include "gdal/gdal_priv.h"
#include "gdal/ogr_geometry.h"
#include "gdal/ogr_spatialref.h"

#include "gdal_raster.h"
#include "sparse_settlements.h"
#include "split_patches.h"


using namespace boost::geometry;
using namespace std;

namespace bg = boost::geometry;
using dpoint = model::d2::point_xy<double>;
using dpolygon = model::polygon<dpoint>;

namespace dd_harp {

//! Copies a OGRLineString to a Boost linestring, reversing points if necessary.
template<typename RING>
void copy_linear_ring(RING& ring, OGRLinearRing* linear_ring) {
    for (const auto& gd_pt: linear_ring) {
        append(ring, dpoint{gd_pt.getX(), gd_pt.getY()});
    }
    if (linear_ring->isClockwise() == (point_order<RING>::value == order_selector::clockwise)) {
        boost::geometry::reverse(ring);
    }
}


/*! Convert a GDAL OGRMultiPolygon into a Boost multi_polygon
 *
 * @param gdal_poly Reads but does not write this. Uses counter-clockwise winding
 *                  of points on outer rings.
 * @return Boost multi-polygon. Uses clockwise winding of points on outer rings.
 */
dmpolygon convert_gdal_to_boost(OGRMultiPolygon* gdal_poly) {
    dmpolygon dpoly;
    dpoly.resize(gdal_poly->getNumGeometries());
    int poly_idx{0};
    for (const auto& outer_polygon: gdal_poly) {
        dpolygon& single{dpoly[poly_idx]};
        OGRLinearRing* exterior = outer_polygon->getExteriorRing();
        copy_linear_ring(single.outer(), exterior);

        single.inners().resize(outer_polygon->getNumInteriorRings());
        for (int inner_idx=0; inner_idx < outer_polygon->getNumInteriorRings(); ++inner_idx) {
            OGRLinearRing* interior = outer_polygon->getInteriorRing(inner_idx);
            copy_linear_ring(single.inners()[inner_idx], interior);
        }
        ++poly_idx;
    }
    return dpoly;
}


void split_patches_retaining_pfpr(
        map<array<int, 2>,PixelData>& settlement_pfpr,
        const vector<double>& settlement_geo_transform,
        const dmpolygon& admin_bg,
        shared_ptr<OGRCoordinateTransformation>& project
) {
    for (auto& pixel_iter: settlement_pfpr) {
        dpolygon pixel_poly;
        // Apply the projection before making the new polygon.
        auto bounds = pixel_bounds<dpoint>(pixel_iter.first, settlement_geo_transform);
        for (auto& bound: bounds) {
            double bx{get<0>(bound)};
            double by{get<1>(bound)};
            project->Transform(1, &bx, &by);
            boost::geometry::set<0>(bound, bx);
            boost::geometry::set<1>(bound, by);
        }
        append(pixel_poly, bounds);

        // If the polygon is split by the side, get the centroid of the inner pieces
        // and the centroid of the outer pieces. That's enough for making patches.
        PixelData& pd{pixel_iter.second};

        double intersect_area{0};
        double total_area{area(pixel_poly)};
        dpoint pix_centroid{0, 0};
        centroid(pixel_poly, pix_centroid);
        dpoint outside_centroid{pix_centroid};
        dpoint total_centroid(pix_centroid);
        bool does_overlap = overlaps(pixel_poly, admin_bg);
        if (does_overlap) {
            bool is_within = within(pixel_poly, admin_bg);
            if (is_within) {
                pd.overlap = Overlap::in;
                pd.area_in = total_area;
                pd.centroid_in = pix_centroid;
            } else {
                deque<dpolygon> output;
                intersection(pixel_poly, admin_bg, output);
                double cx{0};
                double cy{0};
                for (const auto& geom: output) {
                    double part_area = area(geom);
                    intersect_area += part_area;
                    dpoint part_centroid;
                    centroid(geom, part_centroid);
                    cx += bg::get<0>(part_centroid) * part_area;
                    cy += bg::get<1>(part_centroid) * part_area;
                }
                bg::set<0>(pix_centroid, cx / intersect_area);
                bg::set<1>(pix_centroid, cy / intersect_area);

                if (total_area - intersect_area > 0.01 * total_area) {
                    // The centroid of the outside can be computed from the total centroid and inside centroid.
                    bg::set<0>(
                            outside_centroid,
                            (bg::get<0>(total_centroid) * total_area
                             - bg::get<0>(pix_centroid) * intersect_area) / (total_area - intersect_area)
                    );
                    bg::set<1>(
                            outside_centroid,
                            (bg::get<1>(total_centroid) * total_area
                             - bg::get<1>(pix_centroid) * intersect_area) / (total_area - intersect_area)
                    );
                    pd.overlap = Overlap::on;
                    pd.area_in = intersect_area;
                    pd.centroid_in = pix_centroid;
                    pd.area_out = total_area - intersect_area;
                    pd.centroid_out = outside_centroid;
                } else {
                    pd.overlap = Overlap::in;
                    pd.area_in = total_area;
                    pd.centroid_in = pix_centroid;
                }
            }
        } else {
            intersect_area = 0;
            pd.overlap = Overlap::out;
            pd.area_out = total_area;
            pd.centroid_out = pix_centroid;
        }
    }
}
}
