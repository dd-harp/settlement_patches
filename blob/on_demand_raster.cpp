#include <iostream>

#include "gdal/gdal_priv.h"

#include "gdal_raster.h"
#include "on_demand_raster.h"


namespace spacepop {

OnDemandRaster::OnDemandRaster(GDALRasterBand* band, const std::vector<double>& geo_transform)
        : _band(band), _transform(geo_transform) {
    this->_band->GetBlockSize(&_block_size[x], &_block_size[y]);
    this->_size[x] = this->_band->GetXSize();  // Track long, lat vs lat, long.
    this->_size[y] = this->_band->GetYSize();
    for (auto coord: {x, y}) {
        this->_block_cnt[coord] = (this->_size[coord] + _block_size[coord] - 1) / _block_size[coord];
    }
}


double OnDemandRaster::at_coord(double lat_coord, double long_coord) {
    auto ix = pixel_containing({lat_coord, long_coord}, this->_transform);
    return this->at(ix);
}


double OnDemandRaster::at(std::array<int, 2> ix) {
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
}
