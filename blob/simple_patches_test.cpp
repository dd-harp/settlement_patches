#include <array>
#include <vector>

#include "gtest/gtest.h"
#include "simple_patches.h"

using namespace std;
using namespace spacepop;


TEST(PixelsToPolygons, HandlesNullInput) {
    vector<size_t> points;
    auto poly = pixels_to_polygon(points, 100, 4);
    ASSERT_EQ(poly.size(), 0);
}


TEST(PixelsToPolygons, SinglePixel) {
    vector<size_t> points{0};
    auto poly = pixels_to_polygon(points, 100, 4);
    ASSERT_EQ(poly.size(), 4);
}


TEST(PixelsToPolygons, NeighboringLeftRightPixels) {
    vector<size_t> points{0, 1};
    auto poly = pixels_to_polygon(points, 100, 4);
    ASSERT_EQ(poly.size(), 6);
}


TEST(PixelsToPolygons, NeighboringUpDownPixels) {
    vector<size_t> points{0, 100};
    auto poly = pixels_to_polygon(points, 100, 4);
    ASSERT_EQ(poly.size(), 6);
}
