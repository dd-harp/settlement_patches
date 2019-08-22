#ifndef BLOB_GDAL_VECTOR_H
#define BLOB_GDAL_VECTOR_H

#include <exception>
#include <filesystem>

#include "gdal/gdal_priv.h"

namespace spacepop {
    void OpenShapefile(const std::filesystem::path &shapefile_path);
    void WriteVector(const std::filesystem::path &shape_filename);
}

#endif //BLOB_GDAL_VECTOR_H
