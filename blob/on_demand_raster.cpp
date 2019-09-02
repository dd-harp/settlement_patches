#include <cassert>
#include <sstream>

#include "gdal/gdal_priv.h"

#include "gdal_raster.h"
#include "on_demand_raster.h"


namespace dd_harp {

OnDemandRaster::OnDemandRaster(GDALRasterBand* band, const std::vector<double>& geo_transform)
        : _band(band), _transform(geo_transform) {
    this->_band->GetBlockSize(&_block_size[_X], &_block_size[_Y]);
    this->_size[_X] = this->_band->GetXSize();  // Track long, lat vs lat, long.
    this->_size[_X] = this->_band->GetYSize();
    for (auto coord: {_X, _Y}) {
        this->_block_cnt[coord] = (this->_size[coord] + _block_size[coord] - 1) / _block_size[coord];
    }
}


double OnDemandRaster::at_coord(double lat_coord, double long_coord) {
    auto ix = pixel_containing({lat_coord, long_coord}, this->_transform);
    return this->at(ix);
}


double OnDemandRaster::at(std::array<int, 2> ix) {
    for (auto check_bounds: {_X, _Y}) {
        if (ix[check_bounds] < 0 || ix[check_bounds] >= this->_size[check_bounds]) {
            std::stringstream berr;
            berr << "Out of bounds for raster (" << ix[_X] << ", " << ix[_Y] << ")" << std::endl;
            throw std::runtime_error(berr.str());
        }
    }
    std::array<int, 2> block{0, 0};
    for (auto bidx: {_X, _Y}) {
        block[bidx] = ix[bidx] / this->_block_size[bidx];
    }
    auto buffer = this->_buffer.find(block);
    if (buffer == this->_buffer.end()) {
        auto buffer_cnt = this->_block_size[_X] * this->_block_size[_Y];
        auto insert = this->_buffer.emplace(block, std::vector<double>(buffer_cnt));
        buffer = insert.first;
        auto read_succeed = this->_band->ReadBlock(block[_X], block[_Y], &buffer->second[0]);
        if (read_succeed != CPLErr::CE_None) {
            std::stringstream err{"Could not read block ("};
            err << block[_X] << ", " << block[_Y] << ")";
            throw std::runtime_error(err.str());
        }
    }
    auto index_in_block = ix[_X] - _block_size[_X] * block[_X] +
                          _block_size[_X] * (ix[_Y] - _block_size[_Y] * block[_Y]);
    assert(index_in_block >= 0);
    assert(index_in_block < _block_size[_X] * _block_size[_Y]);
    return buffer->second.at(index_in_block);
}
}
