#ifndef BLOB_PIXEL_H
#define BLOB_PIXEL_H

#include <tuple>

//! Use a single integer to identify an (x, y) pixel.
typedef size_t PixelKey;


//! Convert an (x, y) into a single integer.
template<typename Point>
PixelKey PointToPixel(const Point& p, uint32 scan_length) {
    return std::lround(p.y()) * scan_length + std::lround(p.x());
}


template<typename PIXELKEY>
std::array<uint32,2> PixelToPoint(PIXELKEY pixel, uint32 scan_length) {
    uint32 x = pixel % scan_length;
    uint32 y = pixel / scan_length;
    return {x, y};
}

#endif //BLOB_PIXEL_H
