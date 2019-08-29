#include <cmath>
#include <iostream>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

#include "admin_patch.h"
#include "projection.h"

using namespace boost::geometry;

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
        this->_size[x] = this->_band->GetXSize();
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


void CreatePatches(
        OGRMultiPolygon* admin, GDALRasterBand* settlement, GDALRasterBand* pfpr,
        const std::vector<double>& settlement_geo_transform,
        const std::vector<double>& pfpr_geo_transform
        )
{
    const int x{0}, y{1};
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
        for (auto mcoord: {x, y}) {
            settlement_min[mcoord] = std::min(settlement_min[mcoord], bpix[mcoord]);
            settlement_max[mcoord] = std::max(settlement_max[mcoord], bpix[mcoord]);
        }
    }

    auto [project, unproject] = projection_for_lat_long(polygon_bounding_box.MinY, polygon_bounding_box.MinX);

    auto settlement_arr = OnDemandRaster(settlement, settlement_geo_transform);
    auto pfpr_arr = OnDemandRaster(pfpr, pfpr_geo_transform);
    const double cutoff = 0.1;
    using point = model::d2::point_xy<double>;
    point p = {2.0, 3.7};
    using polygon = model::polygon<point>;
    for (int pixel_y=settlement_min[y]; pixel_y < settlement_max[y]; ++pixel_y) {
        for (int pixel_x=settlement_min[x]; pixel_x < settlement_max[x]; ++pixel_x) {
            if (settlement_arr.at({pixel_x, pixel_y}) > cutoff) {
                // if the settlement pop > cutoff.
                polygon pixel_poly;
                append(pixel_poly, pixel_bounds<point>({pixel_x, pixel_y}, settlement_geo_transform));
            }
        }
    }
}

}
