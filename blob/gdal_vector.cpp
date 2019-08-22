#include <iomanip>
#include <iostream>

#include "gdal/ogrsf_frmts.h"

#include "gdal_vector.h"


using namespace std;
namespace fs = std::filesystem;

namespace spacepop {
void OpenShapefile(const std::filesystem::path &shapefile_path)
{
    auto dataset = static_cast<GDALDataset *>(GDALOpenEx(
            shapefile_path.c_str(),
            GDAL_OF_VECTOR | GDAL_OF_READONLY | GDAL_OF_VERBOSE_ERROR,
            nullptr,
            nullptr,
            nullptr
    ));
    if (dataset == nullptr) {
        cout << "Could not load dataset from " << shapefile_path << endl;
        return;
    }
    for (auto layer: dataset->GetLayers()) {
        cout << "layer " << layer->GetName() << endl;
        for (auto& feature: layer) {
            bool comma{false};
            for (auto&& oField: *feature) {
                if (comma) {
                    cout << ",";
                } else {
                    comma = true;
                }
                cout << oField.GetName() << "=";
                switch (oField.GetType()) {
                    case OFTInteger:
                        cout << oField.GetInteger();
                        break;
                    case OFTInteger64:
                        cout << oField.GetInteger64();
                        break;
                    case OFTReal:
                        cout << setprecision(3) << oField.GetDouble();
                        break;
                    case OFTString:
                        cout << oField.GetString();
                        break;
                    default:
                        cout << oField.GetString();
                        break;
                }
            }
            cout << endl;
        }

        for (auto& geom_feature: layer) {
            auto geometry = geom_feature->GetGeometryRef();
            if (geometry != nullptr) {
                cout << "geometry " << geometry->getGeometryName() << " " << geometry->getGeometryType() << endl;
                if (geometry->getGeometryType() == wkbLineString) {
                    auto line_string = dynamic_cast<OGRLineString*>(geometry);
                    std::cout << "is closed: " << line_string->get_IsClosed() << endl;
                    for (auto& point: line_string) {
                        std::cout << "(" << point.getX() << "," << point.getY() << ") ";
                    }
                } else {
                    cout << "unknown geometry type" << endl;
                }
            } else {
                cout << "No geometry for feature" << endl;
            }
        }
    }
}

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