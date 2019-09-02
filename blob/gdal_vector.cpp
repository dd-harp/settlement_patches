#include <iomanip>
#include <iostream>

#include "gdal/ogrsf_frmts.h"

#include "gdal_vector.h"


using namespace std;
namespace fs = boost::filesystem;

namespace dd_harp {
    void print_layer_features(OGRLayer *layer) {
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
    }

    void print_layer_geometry(OGRLayer *layer) {
        for (auto& geom_feature: layer) {
            auto geometry = geom_feature->GetGeometryRef();
            if (geometry != nullptr) {
                cout << "geometry " << geometry->getGeometryName() << " " << geometry->getGeometryType() << endl;
                if (geometry->getGeometryType() == wkbLineString) {
                    auto line_string = dynamic_cast<OGRLineString *>(geometry);
                    cout << "is closed: " << line_string->get_IsClosed() << endl;
                    for (auto &point: line_string) {
                        cout << "(" << point.getX() << "," << point.getY() << ") ";
                    }
                    cout << endl;
                } else if (geometry->getGeometryType() == wkbPolygon) {
                    auto polygon = dynamic_cast<OGRPolygon *>(geometry);
                    cout << "has curve " << polygon->hasCurveGeometry() << " "
                         << "exterior points " << polygon->getExteriorRing()->getNumPoints() << endl;
                } else if (geometry->getGeometryType() == wkbMultiPolygon){
                    auto multi_polygon = dynamic_cast<OGRMultiPolygon*>(geometry);
                    int total_points = 0;
                    for (auto& sub_poly: multi_polygon) {
                        total_points += sub_poly->getExteriorRing()->getNumPoints();
                    }
                    cout << "multipolygon has points " << total_points << endl;
                } else {
                    cout << "unknown geometry type" << endl;
                }
            } else {
                cout << "No geometry for feature" << endl;
            }
        }
    }

    void OpenShapefile(const boost::filesystem::path &shapefile_path)
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
            auto feature_definition = layer->GetLayerDefn();
            for (int field_idx=0; field_idx < feature_definition->GetFieldCount(); ++field_idx) {
                OGRFieldDefn* field_definition = feature_definition->GetFieldDefn(field_idx);
                cout << "\tfeature " << field_definition->GetNameRef() << " "
                    << OGRFieldDefn::GetFieldTypeName(field_definition->GetType()) << endl;
            }

            print_layer_features(layer);
    //        print_layer_geometry(layer);
        }
    }


    void WriteVector(const boost::filesystem::path &shape_filename) {
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