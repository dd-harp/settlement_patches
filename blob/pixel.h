//
// Created by adolgert on 8/7/19.
//

#ifndef BLOB_PIXEL_H
#define BLOB_PIXEL_H

//! Use a single integer to identify an (x, y) pixel.
typedef size_t PixelKey;


//! Convert an (x, y) into a single integer.
template<typename Point>
PixelKey PointToPixel(const Point& p, uint32 scan_length) {
    auto x = p.x();
    auto y = p.y();
    return std::lround(y) * scan_length + std::lround(x);
}

#endif //BLOB_PIXEL_H
