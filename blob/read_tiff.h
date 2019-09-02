#ifndef BLOB_READ_TIFF_H
#define BLOB_READ_TIFF_H

#include <CGAL/Cartesian/Cartesian_base.h>

namespace dd_harp {

template<typename TIFFTYPE>
std::vector<uint32> ImageSizes(TIFFTYPE* tif) {
    uint32 imageWidth=0;
    uint32 imageLength=0;
    uint32 tileWidth=0;
    uint32 tileLength=0;

    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &imageWidth);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imageLength);
    TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileWidth);
    TIFFGetField(tif, TIFFTAG_TILELENGTH, &tileLength);
    return {imageWidth, imageLength, tileWidth, tileLength};
}


/*! Read populations from input.
 *
 * @tparam PointOutputIterator An output iterator over Point objects.
 * @param out The output iterator.
 * @param tif A pointer to an open TIFF file.
 * @param cutoff If the estimated number of people is less than this, ignore the pixel.
 * @param corner Default 0. For debugging, determines how many tiles wide and tall to use.
 */
template<typename PointOutputIterator>
void tiff_input(PointOutputIterator out, TIFF* tif, double cutoff, int corner=0) {
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

    auto sizes = ImageSizes(tif);
    auto xlim = sizes[0];
    auto ylim = sizes[1];
    if (corner > 0) {
        xlim = std::min(sizes[2] * corner, xlim);
        ylim = std::min(sizes[3] * corner, ylim);
    }
    auto buffer = static_cast<double*>(_TIFFmalloc(TIFFTileSize(tif)));
    for (uint32 y=0; y<ylim; y+= sizes[3]) {
        for (uint32 x=0; x<xlim; x+= sizes[2]) {
            // Pull one tile at a time.
            uint32 slice_index=0;  // Only have one slice.
            tsample_t sample_plane=0;  // Data aren't in separate planes.
            TIFFReadTile(tif, buffer, x, y, slice_index, sample_plane);

            // Insert points from that tile into our point set.
            for (uint32 j=0; j<sizes[3]; j++) {
                for (uint32 i=0; i<sizes[2]; i++) {
                    if (buffer[j*sizes[2] + i] > cutoff) {
                        *out++ = Kernel::Point_2{i + x, j + y};
                    }
                }
            }
        }
    }
    _TIFFfree(buffer);
}
}
#endif //BLOB_READ_TIFF_H
