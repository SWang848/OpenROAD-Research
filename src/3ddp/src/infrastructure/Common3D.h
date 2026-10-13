// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#pragma once

#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>

#include "Coordinates.h"  // 使用类型安全的坐标系统

namespace tdp {

// 3D常用类型定义（与Coordinates.h保持一致）
using DBU = int64_t;  // 与TypedCoordinate::v类型一致
using FloatDBU = double;

// ============================================================================
// 3D工具函数（使用类型安全坐标）
// ============================================================================

// 3D曼哈顿距离
inline DBU manhattanDistance3D(const Pointdp3D& p1, const Pointdp3D& p2) {
    return std::abs(p1.x.v - p2.x.v) + std::abs(p1.y.v - p2.y.v) + std::abs(p1.z.v - p2.z.v);
}

// 3D HPWL计算
inline DBU calculateHPWL3D(const std::vector<Pointdp3D>& points) {
    if (points.empty()) return 0;
    
    DBU xMin = points[0].x.v;
    DBU xMax = points[0].x.v;
    DBU yMin = points[0].y.v;
    DBU yMax = points[0].y.v;
    DBU zMin = points[0].z.v;
    DBU zMax = points[0].z.v;
    
    for (const auto& p : points) {
        xMin = std::min(xMin, p.x.v);
        xMax = std::max(xMax, p.x.v);
        yMin = std::min(yMin, p.y.v);
        yMax = std::max(yMax, p.y.v);
        zMin = std::min(zMin, p.z.v);
        zMax = std::max(zMax, p.z.v);
    }
    
    return (xMax - xMin) + (yMax - yMin) + (zMax - zMin);
}

// 2D工具函数（向后兼容）
inline DBU manhattanDistance(const Point& p1, const Point& p2) {
    return std::abs(p1.x.v - p2.x.v) + std::abs(p1.y.v - p2.y.v);
}

inline DBU calculateHPWL(const std::vector<Point>& points) {
    if (points.empty()) return 0;
    
    DBU xMin = points[0].x.v;
    DBU xMax = points[0].x.v;
    DBU yMin = points[0].y.v;
    DBU yMax = points[0].y.v;
    
    for (const auto& p : points) {
        xMin = std::min(xMin, p.x.v);
        xMax = std::max(xMax, p.x.v);
        yMin = std::min(yMin, p.y.v);
        yMax = std::max(yMax, p.y.v);
    }
    
    return (xMax - xMin) + (yMax - yMin);
}

// DBU转微米
inline double dbuToMicron(DBU dbu, int dbu_per_micron) {
    return static_cast<double>(dbu) / dbu_per_micron;
}

// 微米转DBU
inline DBU micronToDBU(double micron, int dbu_per_micron) {
    return static_cast<DBU>(micron * dbu_per_micron + 0.5);
}

// ============================================================================
// 日志工具 - 使用OpenROAD的Logger
// ============================================================================

// // 前向声明OpenROAD的Logger
// namespace utl {
//     class Logger;
// }

// // 日志辅助宏（仿照DPL风格）
// // 在实际使用中，应该通过构造函数传入Logger指针
// #define TDP_LOG_INFO(logger, msg) \
//     if (logger) logger->info(msg)
// #define TDP_LOG_WARN(logger, msg) \
//     if (logger) logger->warn(msg)
// #define TDP_LOG_ERROR(logger, msg) \
//     if (logger) logger->error(msg)
// #define TDP_LOG_DEBUG(logger, msg) \
//     if (logger) logger->debug(msg)

// ============================================================================
// 3D配置结构
// ============================================================================

struct Config3D {
    // 3D特定参数
    DBU inter_layer_spacing = 0;      // 层间间距（从LEF文件提取）
    DBU hbt_pitch = 0;                // 混合键合终端间距（从LEF文件提取，不再硬编码）
    DBU alignment_tolerance = 10;     // 对齐容差
    
    // 层配置
    struct LayerConfig {
        std::string name;
        std::string technology;
        DBU z_offset;                 // Z轴偏移（0=bottom层, 1=upper层）
        DBU row_height;               // 行高度
        DBU site_width;               // site宽度
        bool is_bottom_layer;         // 是否为底层
    };
    
    std::vector<LayerConfig> layers;
    
    // 算法参数
    int maxIterations = 100;
    double temperature = 1.0;
    double coolingRate = 0.95;
    
    // 文件路径
    std::string bottom_lef_file;      // 下层工艺LEF文件（如ASAP7）
    std::string upper_lef_file;       // 上层工艺LEF文件（如NanGate45）
    std::string hbt_lef_file;         // HBT（混合键合）LEF文件（包含hb_layer定义）
    std::string def_file;             // 输入DEF文件
    std::string output_def_file;      // 输出DEF文件
    
    // 清空配置
    void clear() {
        inter_layer_spacing = 0;
        hbt_pitch = 0;
        alignment_tolerance = 10;
        layers.clear();
        maxIterations = 100;
        temperature = 1.0;
        coolingRate = 0.95;
        bottom_lef_file.clear();
        upper_lef_file.clear();
        hbt_lef_file.clear();
        def_file.clear();
        output_def_file.clear();
    }
    
    // 加载配置（TODO: 实现）
    bool loadFromFile(const std::string& filename) {
        // TODO: 从文件加载配置
        return false;
    }
    
    // 打印配置
    // void print() const {
    //     std::cout << "=== 3D Configuration ===" << std::endl;
    //     std::cout << "Inter-layer spacing: " << inter_layer_spacing << " DBU" << std::endl;
    //     std::cout << "HBT pitch: " << hbt_pitch << " DBU" << std::endl;
    //     std::cout << "Alignment tolerance: " << alignment_tolerance << " DBU" << std::endl;
    //     std::cout << "Layers: " << layers.size() << std::endl;
    //     for (const auto& layer : layers) {
    //         std::cout << "  - " << layer.name << " (" << layer.technology 
    //                   << "), z=" << layer.z_offset << ", height=" << layer.row_height 
    //                   << ", width=" << layer.site_width 
    //                   << ", bottom=" << (layer.is_bottom_layer ? "yes" : "no") << std::endl;
    //     }
    //     std::cout << "Files:" << std::endl;
    //     std::cout << "  - Bottom LEF: " << (bottom_lef_file.empty() ? "(not set)" : bottom_lef_file) << std::endl;
    //     std::cout << "  - Upper LEF: " << (upper_lef_file.empty() ? "(not set)" : upper_lef_file) << std::endl;
    //     std::cout << "  - HBT LEF: " << (hbt_lef_file.empty() ? "(not set)" : hbt_lef_file) << std::endl;
    //     std::cout << "  - Input DEF: " << (def_file.empty() ? "(not set)" : def_file) << std::endl;
    //     std::cout << "  - Output DEF: " << (output_def_file.empty() ? "(not set)" : output_def_file) << std::endl;
    // }
};

} // namespace tdp