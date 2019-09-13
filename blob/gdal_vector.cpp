#include <iomanip>
#include <iostream>

#include "gdal/ogrsf_frmts.h"

#include "component_data.h"
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


    struct FeatureDeleter
    {
        void operator()(OGRFeature* f) { OGRFeature::DestroyFeature(f); }
    };


    void WriteVector(const boost::filesystem::path &shape_filename, const std::vector<ComponentData>& patch) {
        // https://gdal.org/tutorials/vector_api_tut.html#writing-to-ogr
        const char *pszDriverName = "ESRI Shapefile";
        GDALDriver *shapefile_driver = GetGDALDriverManager()->GetDriverByName(pszDriverName);
        if (shapefile_driver == nullptr) {
            throw std::runtime_error("Could not load GeoTIFF driver in GDAL.");
        }
        if (shape_filename.extension() != ".shp") {
            throw std::runtime_error("Shape filename does not end in shp");
        }
        auto ds = GDALDatasetUniquePtr(shapefile_driver->Create(
                shape_filename.c_str(), 0, 0, 0, GDT_Unknown, nullptr
        ), GDALDatasetUniquePtrDeleter());
        if (ds == nullptr) {
            throw std::runtime_error("Could not create dataset");
        }
        OGRLayer* layer = ds->CreateLayer("patch", nullptr, wkbPoint, nullptr);
        if (layer == nullptr) {
            throw std::runtime_error("Could not create layer.");
        }

        OGRFieldDefn population_field("Population", OFTReal);
        auto pop_create = layer->CreateField(&population_field);
        if (pop_create != OGRERR_NONE) {
            throw std::runtime_error("Could not create population field");
        }
        OGRFieldDefn pfpr_field("PfPR", OFTReal);
        auto pfpr_create = layer->CreateField(&pfpr_field);
        if (pfpr_create != OGRERR_NONE) {
            throw std::runtime_error("Could not create PfPR field");
        }

        for (const auto& patch_description: patch) {
            OGRPoint pt;
            pt.setX(patch_description.centroid_lat_long[0]);
            pt.setY(patch_description.centroid_lat_long[1]);

            auto feature = std::unique_ptr<OGRFeature, FeatureDeleter>(
                    static_cast<OGRFeature*>(OGRFeature::CreateFeature(layer->GetLayerDefn()))
            );
            feature->SetGeometry(&pt);
            if (layer->CreateFeature(feature.get()) != OGRERR_NONE) {
                throw runtime_error("Could not create feature in layer.");
            }
        }
    }
}
