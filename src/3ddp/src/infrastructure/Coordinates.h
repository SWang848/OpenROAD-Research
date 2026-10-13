// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>

namespace tdp {

// ============================================================================
// TypedCoordinate实现
// ============================================================================

template <typename Tag>
struct TypedCoordinate {
    int64_t v;
    
    TypedCoordinate() : v(0) {}
    TypedCoordinate(int64_t value) : v(value) {}
    
    // Arithmetic operators
    TypedCoordinate operator+(const TypedCoordinate& other) const {
        return TypedCoordinate(v + other.v);
    }
    
    TypedCoordinate operator-(const TypedCoordinate& other) const {
        return TypedCoordinate(v - other.v);
    }
    
    TypedCoordinate operator*(int64_t scalar) const {
        return TypedCoordinate(v * scalar);
    }
    
    TypedCoordinate operator/(int64_t scalar) const {
        return TypedCoordinate(v / scalar);
    }
    
    // Comparison operators
    bool operator==(const TypedCoordinate& other) const { return v == other.v; }
    bool operator!=(const TypedCoordinate& other) const { return v != other.v; }
    bool operator<(const TypedCoordinate& other) const { return v < other.v; }
    bool operator<=(const TypedCoordinate& other) const { return v <= other.v; }
    bool operator>(const TypedCoordinate& other) const { return v > other.v; }
    bool operator>=(const TypedCoordinate& other) const { return v >= other.v; }
    
    // Compound assignment operators
    TypedCoordinate& operator+=(const TypedCoordinate& other) {
        v += other.v;
        return *this;
    }
    
    TypedCoordinate& operator-=(const TypedCoordinate& other) {
        v -= other.v;
        return *this;
    }
    
    TypedCoordinate& operator*=(int64_t scalar) {
        v *= scalar;
        return *this;
    }
    
    TypedCoordinate& operator/=(int64_t scalar) {
        v /= scalar;
        return *this;
    }
    
    // Unary operators
    TypedCoordinate operator-() const {
        return TypedCoordinate(-v);
    }
    
    // Conversion to int64_t (explicit)
    explicit operator int64_t() const { return v; }
    
    // Min and max values
    static TypedCoordinate min() {
        return TypedCoordinate(std::numeric_limits<int64_t>::min());
    }
    
    static TypedCoordinate max() {
        return TypedCoordinate(std::numeric_limits<int64_t>::max());
    }
    
    static TypedCoordinate lowest() {
        return TypedCoordinate(std::numeric_limits<int64_t>::lowest());
    }
};

// Type tags for different coordinate systems (仿照DPL)
struct GridXType;
struct GridYType;
struct DbuXType;
struct DbuYType;
struct DbuZType;  // TDP扩展：3D Z坐标

// Typed coordinate definitions (仿照DPL)
using GridX = TypedCoordinate<GridXType>;
using GridY = TypedCoordinate<GridYType>;
using DbuX = TypedCoordinate<DbuXType>;
using DbuY = TypedCoordinate<DbuYType>;
using DbuZ = TypedCoordinate<DbuZType>;  // TDP扩展：3D Z坐标

// ============================================================================
// 2D坐标结构（与DPL兼容）
// ============================================================================

// 2D点（仿照DPL的DbuPt）
struct Point {
    DbuX x;
    DbuY y;
    
    Point() : x(0), y(0) {}
    Point(DbuX x_, DbuY y_) : x(x_), y(y_) {}
    
    bool operator==(const Point& other) const {
        return x.v == other.x.v && y.v == other.y.v;
    }
    
    bool operator!=(const Point& other) const {
        return !(*this == other);
    }
};

// // 2D矩形（仿照DPL的DbuRect）
// struct Rect {
//     DbuX xMin;
//     DbuY yMin;
//     DbuX xMax;
//     DbuY yMax;
    
//     Rect() : xMin(0), yMin(0), xMax(0), yMax(0) {}
//     Rect(DbuX x1, DbuY y1, DbuX x2, DbuY y2)
//         : xMin(std::min(x1, x2)), yMin(std::min(y1, y2)),
//           xMax(std::max(x1, x2)), yMax(std::max(y1, y2)) {}
    
//     DbuX width() const { return xMax - xMin; }
//     DbuY height() const { return yMax - yMin; }
    
//     Point center() const {
//         return Point((xMin + xMax) / 2, (yMin + yMax) / 2);
//     }
    
//     bool contains(const Point& p) const {
//         return p.x >= xMin && p.x <= xMax && p.y >= yMin && p.y <= yMax;
//     }
    
//     bool overlaps(const Rect& other) const {
//         return !(xMax <= other.xMin || xMin >= other.xMax ||
//                  yMax <= other.yMin || yMin >= other.yMax);
//     }
// };

// Grid point（仿照DPL的GridPt）
struct GridPt {
    GridX x;
    GridY y;
    
    GridPt() : x(0), y(0) {}
    GridPt(GridX x_, GridY y_) : x(x_), y(y_) {}
    
    bool operator==(const GridPt& other) const {
        return x.v == other.x.v && y.v == other.y.v;
    }
};

// ============================================================================
// 3D坐标结构（TDP扩展）
// ============================================================================

// 3D点（TDP扩展：支持bottom层=0, upper层=1）
struct Pointdp3D {
    DbuX x;
    DbuY y;
    DbuZ z;  // 0 = bottom层, 1 = upper层
    
    Pointdp3D() : x(0), y(0), z(0) {}
    Pointdp3D(DbuX x_, DbuY y_, DbuZ z_) : x(x_), y(y_), z(z_) {}
    
    bool operator==(const Pointdp3D& other) const {
        return x.v == other.x.v && y.v == other.y.v && z.v == other.z.v;
    }
    
    bool operator!=(const Pointdp3D& other) const {
        return !(*this == other);
    }
};

// 3D矩形（TDP扩展）
struct Rect3D {
    DbuX xMin;
    DbuY yMin;
    DbuZ zMin;  // 层信息
    DbuX xMax;
    DbuY yMax;
    DbuZ zMax;  // 层信息
    
    Rect3D() : xMin(0), yMin(0), zMin(0), xMax(0), yMax(0), zMax(0) {}
    Rect3D(DbuX x1, DbuY y1, DbuZ z1, DbuX x2, DbuY y2, DbuZ z2)
        : xMin(std::min(x1, x2)), yMin(std::min(y1, y2)), zMin(std::min(z1, z2)),
          xMax(std::max(x1, x2)), yMax(std::max(y1, y2)), zMax(std::max(z1, z2)) {}
    
    DbuX width() const { return xMax - xMin; }
    DbuY height() const { return yMax - yMin; }
    DbuZ depth() const { return zMax - zMin; }
    
    Pointdp3D center() const {
        return Pointdp3D((xMin + xMax) / 2, (yMin + yMax) / 2, (zMin + zMax) / 2);
    }
    
    bool contains(const Pointdp3D& p) const {
        return p.x >= xMin && p.x <= xMax &&
               p.y >= yMin && p.y <= yMax &&
               p.z >= zMin && p.z <= zMax;
    }
    
    bool overlaps(const Rect3D& other) const {
        return !(xMax <= other.xMin || xMin >= other.xMax ||
                 yMax <= other.yMin || yMin >= other.yMax ||
                 zMax <= other.zMin || zMin >= other.zMax);
    }
};

// ============================================================================
// 工具函数（仿照DPL）
// ============================================================================

// 除法函数（仿照DPL的divRound, divCeil, divFloor）
template <typename T>
inline T divRound(T dividend, T divisor) {
    return T((dividend.v + divisor.v / 2) / divisor.v);
}

template <typename T>
inline T divCeil(T dividend, T divisor) {
    return T((dividend.v + divisor.v - 1) / divisor.v);
}

template <typename T>
inline T divFloor(T dividend, T divisor) {
    return T(dividend.v / divisor.v);
}

} // namespace tdp
