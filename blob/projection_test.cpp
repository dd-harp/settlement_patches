#include <vector>

#include "gtest/gtest.h"

#include "projection.h"

using namespace std;
using namespace dd_harp;


/*! Using a Wikipedia example as a test.
 *  https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system#Locating_a_position_using_UTM_coordinates
 *  "The CN Tower is at 43°38′33.24″N 79°23′13.7″W, which is in UTM zone 17,
 *  and the grid position is 630084 m east, 4833438 m north."
 */
TEST(ProjectionForLatLong,Example) {
    double latitude{43 + 38 / 60.0 + 33.24 / 3600.0};  // North, so positive.
    double longitude{-(79 + 23 / 60.0  + 13.7 / 3600.0)};  // West, so minus.
    auto [transform_in, transform_out] = projection_for_lat_long(latitude, longitude);
    double x{longitude}, y{latitude};
    transform_in->Transform(1, &x, &y);
    cout << "takes (" << latitude << ", " << longitude << ") to (" << x << ", " << y << ")" << endl;
    EXPECT_NEAR(x, 630084, 10);  // East
    EXPECT_NEAR(y, 4833438, 10); // North
}


TEST(ProjectionForLatLong,SameOrdering) {
    double latitude{37};  // North, so positive.
    double longitude{1};  // West, so minus.
    auto [transform_in, transform_out] = projection_for_lat_long(latitude, longitude);

    // Make a clockwise square.
    vector<double> longitudes = {37, 37, 37.1, 37.1};
    vector<double> latitudes = {1, 1.1, 1.1, 1};

    transform_in->Transform(4, &longitudes[0], &latitudes[0]);
    for (int i=0; i<4; ++i) {
        cout << longitudes[i] << " " << latitudes[i] << endl;
    }

    // Check that it didn't flip orientation.
    EXPECT_GT(longitudes[3], longitudes[0]);
    EXPECT_GT(latitudes[1], latitudes[0]);
}
