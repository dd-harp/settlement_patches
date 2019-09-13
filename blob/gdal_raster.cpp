#include <iostream>
#include <memory>

#include "gdal_raster.h"

using namespace std;
namespace fs = boost::filesystem;

namespace dd_harp {

    struct DatasetClose {
        void operator()(GDALDataset* ds) const {
            GDALClose(ds);
        }
    };

    shared_ptr<GDALDataset> OpenGeoTiff(const boost::filesystem::path &geotiff_filename) {

        auto parent_directory = geotiff_filename.parent_path();
        auto file_stem = geotiff_filename.stem().string();
        auto without_extension = parent_directory / file_stem;
        cout << "Opening " << geotiff_filename << endl;
        const char *allowed_drivers[3];
        allowed_drivers[0] = "GTiff";
        allowed_drivers[1] = "GeoTIFF";
        allowed_drivers[2] = nullptr;

        auto parent_iter = fs::directory_iterator(parent_directory);
        file_stem.append(".");
        cout << "stem is " << file_stem << endl;
        auto compare_len = file_stem.size();
        const auto sibling_cnt = count_if(begin(parent_iter), end(parent_iter), [&](auto dirent) {
            return dirent.path().filename().string().compare(0, compare_len, file_stem) == 0;
        }) - 1;  // for the file itself.
        cout << "sibling count " << sibling_cnt << endl;
        char *siblings[sibling_cnt + 1];
        size_t sibling_idx = 0;
        for (auto &sibling : fs::directory_iterator(parent_directory)) {
            bool matches = sibling.path().filename().string().compare(0, compare_len, file_stem) == 0;
            bool same_file = sibling.path() == geotiff_filename;
            if (matches && !same_file) {
                auto name = sibling.path().string();
                siblings[sibling_idx] = new char[name.size() + 1];
                strcpy(siblings[sibling_idx++], &name[0]);
                cout << "sibling " << name << endl;
            }
        }
        assert(sibling_idx == sibling_cnt);
        siblings[sibling_cnt] = nullptr;

        auto dataset = static_cast<GDALDataset *>(GDALOpenEx(
                geotiff_filename.c_str(),
                GDAL_OF_RASTER | GDAL_OF_READONLY | GDAL_OF_VERBOSE_ERROR,
                allowed_drivers,
                nullptr,
                siblings
        ));
        if (dataset == nullptr) {
            if (!fs::exists(geotiff_filename)) {
                std::stringstream msg{"File path not found "};
                msg << geotiff_filename;
                throw std::runtime_error(msg.str());
            }
            std::stringstream msg{"Could not load the file "};
            msg << geotiff_filename << " Did the code call GDALAllRegister?";
            throw std::runtime_error(msg.str());
        }
        return shared_ptr<GDALDataset>(dataset, DatasetClose());
    }

    std::array<int, 2> pixel_containing(std::array<double, 2> coord, const std::vector<double>& transform) {
        return {
                static_cast<int>(std::lround(std::floor((1 / transform[1]) * (coord[0] - transform[0])))),
                static_cast<int>(std::lround(std::floor((1 / transform[5]) * (coord[1] - transform[3]))))
        };
    }

}
