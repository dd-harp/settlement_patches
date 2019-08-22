//
// Created by adolgert on 8/16/19.
//

#ifndef BLOB_GDAL_VECTOR_H
#define BLOB_GDAL_VECTOR_H

#include <exception>
#include <filesystem>

#include "gdal/gdal_priv.h"

namespace spacepop {
    void WriteVector(const std::filesystem::path &shape_filename) {
        // https://gdal.org/tutorials/vector_api_tut.html#writing-to-ogr
        GDALDriver *shapefile_driver = GetGDALDriverManager()->GetDriverByName("GTiff");
        if (shapefile_driver == nullptr) {
            throw std::runtime_error("Could not load GeoTIFF driver in GDAL.");
        }
        GDALDataset* ds = shapefile_driver->Create(
                shape_filename.c_str(), 0, 0, 0, GDT_Unknown, nullptr
                );
        if (ds == nullptr) {
            throw std::runtime_error("Could not create dataset");
        }
        OGRLayer* layer = ds->CreateLayer("patch", nullptr, wkbPoint, nullptr);
        if (layer == nullptr) {
            throw std::runtime_error("Could not create layer.");
        }

    }
}

#endif //BLOB_GDAL_VECTOR_H
