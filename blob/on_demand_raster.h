#ifndef BLOB_ON_DEMAND_RASTER_H
#define BLOB_ON_DEMAND_RASTER_H

#include <array>
#include <map>
#include <vector>

class GDALRasterBand;

namespace dd_harp {

/*! Read GDAL Raster data set block-by-block, as pixels are requested.
 *  GDAL has a least-recently-used cache. This asks for the blocks it
 *  needs when it needs them. It doesn't throw out any blocks because there
 *  is LRU underneath.
 */
    class OnDemandRaster {
        const int _X{0}, _Y{1};  // To make notation clearer.
        GDALRasterBand *_band;
        const std::vector<double> &_transform;
        std::map <std::array<int, 2>, std::vector<double>> _buffer;

    public:
        OnDemandRaster(GDALRasterBand *band, const std::vector<double> &geo_transform);
        double at_coord(double lat_coord, double long_coord);
        double at(std::array<int, 2> ix);

    private:
        std::array<int, 2> _size;  // Total pixels in raster as x, y
        std::array<int, 2> _block_size;  // Size of each raster block.
        std::array<int, 2> _block_cnt;  // Number of raster blocks in x and y.
    };
}
#endif //BLOB_ON_DEMAND_RASTER_H
