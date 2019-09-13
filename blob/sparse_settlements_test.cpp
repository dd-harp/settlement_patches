#include <sstream>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include "gdal/ogr_geometry.h"
#include "gtest/gtest.h"

#include "sparse_settlements.h"

using namespace std;
using namespace dd_harp;

TEST(GdalGeometryMinMax,Baseline) {
    namespace geom = boost::geometry;
    typedef geom::model::d2::point_xy<double> point_type;
    geom::model::polygon<point_type> polygon;
    geom::append(geom::exterior_ring(polygon), geom::make<point_type>(32, 1));
    geom::append(geom::exterior_ring(polygon), geom::make<point_type>(32, 2));
    geom::append(geom::exterior_ring(polygon), geom::make<point_type>(33, 2));
    geom::append(geom::exterior_ring(polygon), geom::make<point_type>(32, 2));

    OGRSpatialReference lat_long_srs;
    lat_long_srs.SetWellKnownGeogCS("WGS84");
    std::stringstream poly_text;
    poly_text << boost::geometry::wkt(polygon);
    OGRGeometry* pg;
    OGRGeometryFactory::createFromWkt(poly_text.str().c_str(), &lat_long_srs, &pg);
    ASSERT_EQ(pg->getGeometryType(), wkbPolygon);

    cout << "Created geometry " <<
        pg->getGeometryName() << " and type " << pg->getGeometryType() <<
        " from " << poly_text.str() <<endl;

    OGREnvelope bounds;
    pg->getEnvelope(&bounds);
    EXPECT_GT(bounds.MinX, 10.0);
    EXPECT_LT(bounds.MinY, 10.0);

    OGRMultiPolygon mpg;
    mpg.addGeometry(pg);
    ASSERT_EQ(mpg.getGeometryType(), wkbMultiPolygon);
    vector<double> geo_transform = {29.5735, 0.000277778, 0, 4.22813, 0, -0.000277778};
    auto min_max = gdal_geometry_min_max(&mpg, geo_transform);
    cout << "min " << min_max.first[0] << " " << min_max.first[1] << endl;
    cout << "max " << min_max.second[0] << " " << min_max.second[1] << endl;
    OGRGeometryFactory::destroyGeometry(pg);
}
