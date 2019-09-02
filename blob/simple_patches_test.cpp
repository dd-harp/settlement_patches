#include <array>
#include <vector>

#include "gtest/gtest.h"
#include "simple_patches.h"

using namespace std;
using namespace dd_harp;


TEST(PixelsToPolygons, HandlesNullInput) {
    vector<size_t> points;
    auto poly = pixels_to_polygon(points, 100, 4);
    ASSERT_EQ(poly.size(), 0);
}

//
//TEST(PixelsToPolygons, SinglePixel) {
//    vector<size_t> points{0};
//    auto poly = pixels_to_polygon(points, 100, 4);
//    ASSERT_EQ(poly.size(), 4);
//}
//
//
//TEST(PixelsToPolygons, NeighboringLeftRightPixels) {
//    vector<size_t> points{0, 1};
//    auto poly = pixels_to_polygon(points, 100, 4);
//    ASSERT_EQ(poly.size(), 6);
//}
//
//
//TEST(PixelsToPolygons, NeighboringUpDownPixels) {
//    vector<size_t> points{0, 100};
//    auto poly = pixels_to_polygon(points, 100, 4);
//    ASSERT_EQ(poly.size(), 6);
//}
//
//
//TEST(PixelsToPolygons, HasAHole) {
//    vector<size_t> points;
//    int n{20};
//    for (int horizontal=0; horizontal < n; horizontal++) {
//        points.push_back(horizontal);
//        points.push_back((n - 1) * 100 + horizontal);
//    }
//
//    for (int vertical=1; vertical < n - 1; vertical++) {
//        points.push_back(vertical * 100);
//        points.push_back(vertical * 100 + n - 1);
//    }
//    auto poly = pixels_to_polygon(points, 100, 8);
//}
