#include <vector>

#include "gtest/gtest.h"

#include "gdal_raster.h"

using namespace std;
using namespace dd_harp;


TEST(PixelContaining,Baseline) {
    vector<double> geo_transform = {29.5735, 0.000277778, 0, 4.22813, 0, -0.000277778};
    auto p = pixel_containing(
            {29.5735, 4.22813}, geo_transform);
    EXPECT_EQ(p[0], 0);
    EXPECT_EQ(p[1], 0);

    auto p1 = pixel_containing(
            {29.5735+0.00027, 4.22813-0.00027}, geo_transform);
    EXPECT_EQ(p1[0], 0);
    EXPECT_EQ(p1[1], 0);

    auto p4 = pixel_containing(
            {29.5735, 4.22813 - 0.00028}, geo_transform);
    EXPECT_EQ(p4[0], 0);
    EXPECT_EQ(p4[1], 1);
}


TEST(PixelContaining,Negative)
{
    vector<double> geo_transform = {29.5735, 0.000277778, 0, 4.22813, 0,
                                    -0.000277778};
    auto p = pixel_containing(
            {29.5735-0.00027, 4.22813}, geo_transform);
    EXPECT_EQ(p[0], -1);
    EXPECT_EQ(p[1], 0);

    auto p1 = pixel_containing(
            {29.5735, 4.22813+0.00027}, geo_transform);
    EXPECT_EQ(p1[0], 0);
    EXPECT_EQ(p1[1], -1);
}


TEST(PixelCoord,Translate)
{
    vector<double> geo_transform = {5, 1, 0, 9, 0, 1};
    auto p0 = pixel_coord<array<double, 2>>({3, 7}, geo_transform);
    EXPECT_FLOAT_EQ(p0[0], 8);
    EXPECT_FLOAT_EQ(p0[1], 16);
}


TEST(PixelCoord,Scale)
{
    vector<double> geo_transform = {0, 3, 0, 0, 0, 2};
    auto p0 = pixel_coord<array<double, 2>>({3, 7}, geo_transform);
    EXPECT_FLOAT_EQ(p0[0], 9);
    EXPECT_FLOAT_EQ(p0[1], 14);
}


TEST(PixelCoord,Shear)
{
    vector<double> geo_transform = {0, 1, 4, 0, 3, 1};
    auto p0 = pixel_coord<array<double, 2>>({3, 7}, geo_transform);
    EXPECT_FLOAT_EQ(p0[0], 31);
    EXPECT_FLOAT_EQ(p0[1], 16);
}
