#include <iostream>
#include <sstream>

#include "gdal/gdal_priv.h"

#include "projection.h"

using namespace std;

namespace dd_harp {

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
        // You can see the zones on https://upload.wikimedia.org/wikipedia/commons/e/ed/Utm-zones.jpg
        // but go to georepository.com to get the exact zone, here 36N.
        uganda_projected_srs.SetUTM(36, South);
        cout << "linear units for projection " <<
            uganda_projected_srs.GetLinearUnits(nullptr) << " [m]" << endl;

        auto transform_in = OGRCreateCoordinateTransformation(purely_lat_long_srs, &uganda_projected_srs);
        return shared_ptr<OGRCoordinateTransformation>(transform_in, CoordinateTransformClose());
    }

    tuple<shared_ptr<OGRCoordinateTransformation>,shared_ptr<OGRCoordinateTransformation>>
    projection_for_lat_long(double latitude, double longitude) {
        OGRSpatialReference pure_lat_long;
        // GDAL 3.0 has an axis mapping strategy to allow long-lat.
        pure_lat_long.SetWellKnownGeogCS("WGS84");
        cout << "treats lat as long " << pure_lat_long.EPSGTreatsAsLatLong() << endl;
        int target_utm = static_cast<int>(lround(ceil((longitude + 180) / 6)));
        int North = (latitude > 0) ? 1 : 0;
        OGRSpatialReference target_srs;
        stringstream utm_name;
        utm_name << "UTM" << target_utm << ((North == 1) ? "North" : "South");
        cout << "projection " << utm_name.str() << endl;
        target_srs.SetProjCS(utm_name.str().c_str());
        target_srs.SetUTM(target_utm, North);
        auto transform_in = shared_ptr<OGRCoordinateTransformation>(
                OGRCreateCoordinateTransformation(&pure_lat_long, &target_srs));
        auto transform_out = shared_ptr<OGRCoordinateTransformation>(
                OGRCreateCoordinateTransformation(&target_srs, &pure_lat_long));
        return {transform_in, transform_out};
    }
}
