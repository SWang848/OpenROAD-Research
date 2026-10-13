#pragma once

#include "HBTTypes.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace hbta {

// 表示一个轴对齐矩形区域，供四叉树划分与范围过滤使用。
struct HBTRegion {
    double x_min = 0.0;
    double y_min = 0.0;
    double x_max = 0.0;
    double y_max = 0.0;

    HBTRegion() = default;

    HBTRegion(double x0, double y0, double x1, double y1)
        : x_min(std::min(x0, x1)),
          y_min(std::min(y0, y1)),
          x_max(std::max(x0, x1)),
          y_max(std::max(y0, y1)) {}

    bool isValid() const {
        return x_min <= x_max && y_min <= y_max;
    }

    double width() const {
        return x_max - x_min;
    }

    double height() const {
        return y_max - y_min;
    }

    double centerX() const {
        return (x_min + x_max) * 0.5;
    }

    double centerY() const {
        return (y_min + y_max) * 0.5;
    }

    bool containsPoint(const Point& point) const {
        return point.x >= x_min
            && point.x <= x_max
            && point.y >= y_min
            && point.y <= y_max;
    }

    bool containsRegion(const HBTRegion& other) const {
        return other.x_min >= x_min
            && other.x_max <= x_max
            && other.y_min >= y_min
            && other.y_max <= y_max;
    }

    bool overlaps(const HBTRegion& other) const {
        return !(other.x_max < x_min
              || other.x_min > x_max
              || other.y_max < y_min
              || other.y_min > y_max);
    }
};

inline std::array<HBTRegion, 4> splitHBTRegion(const HBTRegion& region) {
    const double x_mid = region.centerX();
    const double y_mid = region.centerY();
    const double left_x_max
        = std::nextafter(x_mid, -std::numeric_limits<double>::infinity());
    const double bottom_y_max
        = std::nextafter(y_mid, -std::numeric_limits<double>::infinity());

    return {
        HBTRegion(region.x_min, y_mid, left_x_max, region.y_max),        // top-left
        HBTRegion(x_mid, y_mid, region.x_max, region.y_max),             // top-right
        HBTRegion(region.x_min, region.y_min, left_x_max, bottom_y_max), // bottom-left
        HBTRegion(x_mid, region.y_min, region.x_max, bottom_y_max)       // bottom-right
    };
}

} // namespace hbta
