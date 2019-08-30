#include <cmath>
#include <deque>
#include <iostream>

#include "boost/geometry.hpp"
#include "boost/geometry/geometries/point_xy.hpp"
#include "boost/geometry/geometries/polygon.hpp"
#include "boost/geometry/geometries/segment.hpp"
#include "boost/polygon/voronoi.hpp"
#include "boost/polygon/polygon.hpp"

#include "admin_patch.h"
#include "projection.h"

using namespace boost::geometry;
using namespace std;

namespace bg = boost::geometry;
using dpoint = model::d2::point_xy<double>;
using dpolygon = model::polygon<dpoint>;
using dmpolygon = model::multi_polygon<dpolygon>;

using ipoint = model::d2::point_xy<int>;
using isegment = model::segment<ipoint>;

using boost::polygon::voronoi_builder;
using boost::polygon::voronoi_diagram;

namespace boost::polygon {

    template <>
    struct geometry_concept<ipoint> { typedef point_concept type; };

    template <>
    struct point_traits<ipoint> {
        typedef int coordinate_type;

        static inline coordinate_type get(const ipoint& point, orientation_2d orient) {
            return (orient == HORIZONTAL) ? bg::get<0>(point) : bg::get<1>(point);
        }
    };

    template <>
    struct geometry_concept<isegment> { typedef segment_concept type; };

    template <>
    struct segment_traits<isegment> {
        typedef int coordinate_type;
        typedef ipoint point_type;

        static inline point_type get(const isegment& segment, direction_1d dir) {
            if (dir.to_int()) {
                return {bg::get<1, 0>(segment), bg::get<1, 1>(segment)};
            } else {
                return {bg::get<0,0>(segment), bg::get<0,1>(segment)};
            }
        }
    };
}

namespace spacepop {
// These templates convince boost::geometry to treat
// an array of two doubles as a Point type.
//namespace boost::geometry::traits {
//    template<> struct tag<std::array<double, 2>>
//    { typedef point_tag type; };
//    template<> struct coordinate_type<std::array<double, 2>>
//    { typedef double type; };
//    template<> struct coordinate_system<std::array<double, 2>>
//    { typedef cs::cartesian type; };
//    template<> struct dimension<std::array<double, 2>> : boost::mpl::int_<2> {};
//    template<>
//    struct access<std::array<double, 2>, 0> {
//        static double get(std::array<double, 2> const &p) {
//            return p[0];
//        }
//
//        static void set(std::array<double, 2> &p, double value) {
//            p[0] = value;
//        }
//    };
//
//    template<>
//    struct access<std::array<double, 2>, 1> {
//        static double get(std::array<double, 2> const &p) {
//            return p[1];
//        }
//
//        static void set(std::array<double, 2> &p, double value) {
//            p[1] = value;
//        }
//    };
//}

const int X{0}, Y{1};

//! Convert a lat-long into a specific pixel.
std::array<int, 2> pixel_containing(std::array<double, 2> coord, const std::vector<double>& transform) {
    return {
            static_cast<int>(std::lround(std::floor((1 / transform[1]) * (coord[0] - transform[0])))),
            static_cast<int>(std::lround(std::floor((1 / transform[5]) * (coord[1] - transform[3]))))
    };
}

//! Convert a pixel corner into a lat-long
template<typename POINTISH>
POINTISH pixel_coord(std::array<int, 2> pixel, const std::vector<double>& transform) {
    return {
        transform[0] + pixel[0] * transform[1] + pixel[1] * transform[2],
        transform[3] + pixel[1] * transform[4] + pixel[5]
    };
}


//! Convert a pixel into its four corners as lat-long.
//  Boost::polygon likes to be clockwise, so these are clockwise.
template<typename POINTISH>
std::vector<POINTISH> pixel_bounds(std::array<int, 2> pixel, const std::vector<double>& transform) {
    return {
        pixel_coord<POINTISH>(pixel, transform),
        pixel_coord<POINTISH>({pixel[0], pixel[1] + 1}, transform),
        pixel_coord<POINTISH>({pixel[0] + 1, pixel[1] + 1}, transform),
        pixel_coord<POINTISH>({pixel[0] + 1, pixel[1]}, transform)
    };
}


/*! Convert a GDAL OGRMultiPolygon into a Boost multi_polygon
 *
 * @param gdal_poly Reads but does not write this.
 * @return Boost multi_polygon
 */
dmpolygon convert(OGRMultiPolygon const* gdal_poly) {
    dmpolygon dpoly;
    dpoly.resize(gdal_poly->getNumGeometries());
    int poly_idx{0};
    for (const auto& child_poly: gdal_poly) {
        dpolygon& single{dpoly[poly_idx]};
        OGRLinearRing const* exterior = child_poly->getExteriorRing();
        bool outer_clock = exterior->isClockwise();  // XXX track clockwise or CCW
        for (const auto& gd_pt: exterior) {
            append(single.outer(), dpoint{gd_pt.getX(), gd_pt.getY()});
        }
        single.inners().resize(child_poly->getNumInteriorRings());
        for (int inner=0; inner < child_poly->getNumInteriorRings(); ++inner) {
            OGRLinearRing const* interior = child_poly->getInteriorRing(inner);
            for (const auto& in_pt: interior) {
                append(single.inners()[inner], dpoint{in_pt.getX(), in_pt.getY()});
            }
        }
        ++poly_idx;
    }
    return dpoly;
}


/*! Read GDAL Raster data set block-by-block, as pixels are requested.
 *  GDAL has a least-recently-used cache. This asks for the blocks it
 *  needs when it needs them. It doesn't throw out any blocks because there
 *  is LRU underneath.
 */
class OnDemandRaster
{
    const int x{0}, y{1};  // To make notation clearer.

    GDALRasterBand* _band;
    const std::vector<double>& _transform;
    std::map<std::array<int,2>,std::vector<double>> _buffer;
public:
    OnDemandRaster(GDALRasterBand* band, const std::vector<double>& geo_transform)
    : _band(band), _transform(geo_transform) {
        this->_band->GetBlockSize(&_block_size[x], &_block_size[y]);
        this->_size[x] = this->_band->GetXSize();  // Track long, lat vs lat, long.
        this->_size[y] = this->_band->GetYSize();
        for (auto coord: {x, y}) {
            this->_block_cnt[coord] = (this->_size[coord] + _block_size[coord] - 1) / _block_size[coord];
        }
    }
    double at_coord(double lat_coord, double long_coord) {
        auto ix = pixel_containing({lat_coord, long_coord}, this->_transform);
        return this->at(ix);
    }

    double at(std::array<int, 2> ix) {
        for (auto check_bounds: {x, y}) {
            if (ix[check_bounds] < 0 || ix[check_bounds] >= this->_size[check_bounds]) {
                std::cout << "Out of bounds for raster (" << ix[x] << ", " << ix[y] << ")" << std::endl;
                return 0.0;
            }
        }
        std::array<int, 2> block{0, 0};
        for (auto bidx: {x, y}) {
            block[bidx] = ix[bidx] / this->_block_size[bidx];
        }
        auto buffer = this->_buffer.find(block);
        if (buffer == this->_buffer.end()) {
            auto buffer_cnt = this->_block_size[x] * this->_block_size[y];
            auto insert = this->_buffer.emplace(block, std::vector<double>(buffer_cnt));
            buffer = insert.first;
            auto read_succeed = this->_band->ReadBlock(block[x], block[y], &buffer->second);
            if (read_succeed != CE_None) {
                std::cout << "Could not read block (" << block[x] << ", " << block[y] << ")" << std::endl;
            }
        }
        return (buffer->second[
                    ix[x] - _block_size[x] * block[x] +
                    _block_size[x] * (ix[y] - _block_size[y] * block[y])
                    ]);
    }
private:
    std::array<int,2> _size;  // Total pixels in raster as x, y
    std::array<int,2> _block_size;  // Size of each raster block.
    std::array<int,2> _block_cnt;  // Number of raster blocks in x and y.
};


/*! Finds the minimum and maximum of a polygon within a particular projection transform.
 *
 * @param admin A OGRGeometry.
 * @param settlement_geo_transform Describes how coordinates relate to pixels for GDAL projections.
 * @return
 */
pair<array<int, 2>, array<int, 2>>
gdal_geometry_min_max(OGRMultiPolygon* admin, const std::vector<double>& settlement_geo_transform) {

    OGREnvelope polygon_bounding_box;
    admin->getEnvelope(&polygon_bounding_box);
    std::cout << "polygon bbox ((" << polygon_bounding_box.MinX << ", " << polygon_bounding_box.MaxX << "), (("
              << polygon_bounding_box.MinY << ", " << polygon_bounding_box.MaxY << "))" << std::endl;

    std::vector<std::array<int,2>> bounding_pixels;
    for (auto ex: {polygon_bounding_box.MinX, polygon_bounding_box.MaxX}) {
        for (auto ey: {polygon_bounding_box.MinY, polygon_bounding_box.MaxY}) {
            bounding_pixels.push_back(pixel_containing({ex, ey}, settlement_geo_transform));
        }
    }
    std::array<int,2> settlement_min{std::numeric_limits<int>::max(), std::numeric_limits<int>::max()};
    std::array<int,2> settlement_max{0, 0};
    for (auto bpix: bounding_pixels) {
        for (auto mcoord: {X, Y}) {
            settlement_min[mcoord] = std::min(settlement_min[mcoord], bpix[mcoord]);
            settlement_max[mcoord] = std::max(settlement_max[mcoord], bpix[mcoord]);
        }
    }
    return {settlement_min, settlement_max};
}


enum class Overlap { in, out, on };
struct PixelData {
    double pfpr;
    double pop;
    Overlap overlap;
    dpoint centroid_in;
    dpoint centroid_out;
    double area_in;
    double area_out;
};


map<array<int, 2>,PixelData>
sparse_settlements(
        OnDemandRaster& settlement_arr,
        OnDemandRaster& pfpr_arr,
        const pair<array<int, 2>, array<int, 2>>& settlement_min_max,
        const vector<double>& settlement_geo_transform,
        const double cutoff
        ) {
    const auto& [settlement_min, settlement_max] = settlement_min_max;
    map<array<int, 2>,PixelData> settlement_pfpr;
    for (int pixel_y=settlement_min[Y]; pixel_y < settlement_max[Y]; ++pixel_y) {
        for (int pixel_x = settlement_min[X]; pixel_x < settlement_max[X]; ++pixel_x) {
            double pixel_pop = settlement_arr.at({pixel_x, pixel_y});
            if (pixel_pop > cutoff) {
                auto settlement_coord = pixel_coord<array<double, 2>>({pixel_x, pixel_y}, settlement_geo_transform);
                PixelData pd{pfpr_arr.at_coord(pixel_x, pixel_y), pixel_pop};
                array<int, 2> create_pixel{pixel_x, pixel_y};
                settlement_pfpr.insert(make_pair(create_pixel, pd));
            }
        }
    }
    return settlement_pfpr;
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


void create_neighbor_graph(map<array<int, 2>,PixelData>& settlement_pfpr) {
    array<double, 2> bmin{numeric_limits<double>::max(), numeric_limits<double>::max()};
    array<double, 2> bmax{numeric_limits<double>::min(), numeric_limits<double>::min()};
    for (const auto& [bounds_p, data_p]: settlement_pfpr) {
        if (data_p.overlap == Overlap::in || data_p.overlap == Overlap::on) {
            auto cpx = data_p.centroid_in.get<0>();
            bmin[0] = min(cpx, bmin[0]);
            bmax[0] = min(cpx, bmax[0]);
            auto cpy = data_p.centroid_in.get<1>();
            bmin[1] = min(cpx, bmin[1]);
            bmax[1] = min(cpx, bmax[1]);
        }
        if (data_p.overlap == Overlap::out || data_p.overlap == Overlap::on) {
            auto cpx = data_p.centroid_out.get<0>();
            bmin[0] = min(cpx, bmin[0]);
            bmax[0] = min(cpx, bmax[0]);
            auto cpy = data_p.centroid_out.get<1>();
            bmin[1] = min(cpx, bmin[1]);
            bmax[1] = min(cpx, bmax[1]);
        }
    }
    array<double, 2> scale = { (1<<30) / (bmax[0] - bmin[0]), (1<<30) / (bmax[1] - bmin[1])};

    std::vector<ipoint> points;
    for (const auto& [pix_key, sd]: settlement_pfpr) {
        if (sd.overlap == Overlap::in || sd.overlap == Overlap::on) {
            auto cpx = sd.centroid_in.get<0>();
            auto cpy = sd.centroid_in.get<1>();
            points.emplace_back(ipoint{
                    static_cast<int>(lround((cpx - bmin[0]) * scale[0])),
                    static_cast<int>(lround((cpy - bmin[1]) * scale[1]))
            });
        }
        if (sd.overlap == Overlap::out || sd.overlap == Overlap::on) {
            auto cpx = sd.centroid_out.get<0>();
            auto cpy = sd.centroid_out.get<1>();
            points.emplace_back(ipoint{
                    static_cast<int>(lround((cpx - bmin[0]) * scale[0])),
                    static_cast<int>(lround((cpy - bmin[1]) * scale[1]))
            });
        }
    }
    std::vector<isegment> segments;
    voronoi_diagram<double> vd;
    construct_voronoi(points.begin(), points.end(), segments.begin(), segments.end(), &vd);

    // Make a graph

    // Walk edges to turn them into graph links.
    for (auto edge_iter = vd.edges().begin(); edge_iter != vd.edges().end(); ++edge_iter) {
        if (edge_iter->is_primary()) {
            ;
        }
    }
}


void CreatePatches(
        OGRMultiPolygon* admin, GDALRasterBand* settlement, GDALRasterBand* pfpr,
        const std::vector<double>& settlement_geo_transform,
        const std::vector<double>& pfpr_geo_transform
        )
{
    const int x{0}, y{1};
    auto settlement_min_max = gdal_geometry_min_max(admin, settlement_geo_transform);
    auto settlement_arr = OnDemandRaster(settlement, settlement_geo_transform);
    auto pfpr_arr = OnDemandRaster(pfpr, pfpr_geo_transform);


    const double cutoff = 0.1;
    map<array<int, 2>,PixelData> settlement_pfpr = sparse_settlements(
            settlement_arr, pfpr_arr, settlement_min_max, settlement_geo_transform, cutoff
            );

    // Work in projection where units are meters.
    OGREnvelope polygon_bounding_box;
    admin->getEnvelope(&polygon_bounding_box);
    auto [project, unproject] = projection_for_lat_long(polygon_bounding_box.MinY, polygon_bounding_box.MinX);

    // Transform the incoming multipolygon in place, not a copy.
    admin->transform(project.get());
    // Use Boost Polygon because it lets us create things and intersect and delaunay them.
    auto admin_bg = convert(admin);

    split_patches_retaining_pfpr(settlement_pfpr, settlement_geo_transform, admin_bg, project);

    // Should return a graph with an index into settlement_pfpr.
    create_neighbor_graph(settlement_pfpr);

    // Cluster on the graph, excluding nodes that are outside the polygon.
}

}
