#include <iostream>

#include "gdal/gdal_priv.h"

#include "projection.h"

using namespace std;

namespace spacepop {

    struct CoordinateTransformClose {
        void operator()(OGRCoordinateTransformation* ct) const {
            OGRCoordinateTransformation::DestroyCT(ct);
        }
    };

    shared_ptr<OGRCoordinateTransformation> reproject(OGRSpatialReference* purely_lat_long_srs) {
        // srs = Spatial Reference System
        OGRSpatialReference uganda_projected_srs;

        uganda_projected_srs.SetProjCS("UTM 36N (WGS84) for Uganda, Kenya");
        uganda_projected_srs.SetWellKnownGeogCS("WGS84");
        int South{0};
        uganda_projected_srs.SetUTM(17, South);
        cout << "linear units for projection " <<
            uganda_projected_srs.GetLinearUnits(nullptr) << " [m]" << endl;

        auto transform_in = OGRCreateCoordinateTransformation(purely_lat_long_srs, &uganda_projected_srs);
        return shared_ptr<OGRCoordinateTransformation>(transform_in, CoordinateTransformClose());
    }
}
