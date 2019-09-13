#ifndef BLOB_GDAL_VECTOR_H
#define BLOB_GDAL_VECTOR_H

#include <exception>

#include "boost/filesystem.hpp"

#include "gdal/gdal_priv.h"

namespace dd_harp {
    void OpenShapefile(const boost::filesystem::path &shapefile_path);

    class ComponentData;
    void WriteVector(const boost::filesystem::path &shape_filename, const std::vector<ComponentData>& patch);
}

#endif //BLOB_GDAL_VECTOR_H
