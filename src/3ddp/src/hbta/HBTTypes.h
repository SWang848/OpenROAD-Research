#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <cmath>
#include <cstdlib>

namespace hbta {

// 基本点结构
struct Point {
    double x;
    double y;
    
    Point() : x(0.0), y(0.0) {}
    Point(double x_val, double y_val) : x(x_val), y(y_val) {}
    
    double distanceTo(const Point& other) const {
        double dx = x - other.x;
        double dy = y - other.y;
        return std::sqrt(dx * dx + dy * dy);
    }
    
    double manhattanDistanceTo(const Point& other) const {
        return std::abs(x - other.x) + std::abs(y - other.y);
    }
};

// HBT位置信息
struct HBTInfo {
    std::string id;           // 如 "hbt_0_0"
    int row;                  // 行索引
    int col;                  // 列索引
    Point location;           // 中心位置（DBU）
    bool used;                // 是否已被使用
    std::string connected_net; // 连接的网络名称
    
    HBTInfo() : row(-1), col(-1), used(false) {}
    HBTInfo(const std::string& id_val, int r, int c, const Point& loc)
        : id(id_val), row(r), col(c), location(loc), used(false) {}
};

// 跨层网络连接
struct InterLayerConnection {
    // 网络信息
    std::string net_name;
    bool is_constant_net = false;
    
    // 驱动单元信息（如果有）
    std::string driver_cell_name;
    Point driver_location;
    std::string driver_layer;  // "bottom" or "upper"
    
    // 负载单元信息
    std::vector<std::string> load_cell_names;
    std::vector<Point> load_locations;
    std::vector<std::string> load_layers;
    
    // HBT信息
    Point suggested_hbt_location;  // 建议的HBT位置
    bool hbt_allocated = false;
    std::string allocated_hbt_id;  // 分配的HBT ID
    
    // 历史 HBT 位置标记及坐标，用于计算迁移代价。
    bool has_previous_hbt_location = false;
    Point previous_hbt_location;

    // 计算方法
    double calculateTotalWirelength() const {
        double length = 0.0;
        
        if (!driver_cell_name.empty()) {
            length += driver_location.manhattanDistanceTo(suggested_hbt_location);
        }
        
        for (size_t i = 0; i < load_cell_names.size(); ++i) {
            length += suggested_hbt_location.manhattanDistanceTo(load_locations[i]);
        }
        
        return length;
    }
    
    // 检查是否有跨层连接
    bool hasInterLayerConnection() const {
        // 检查驱动和负载是否在不同层
        if (!driver_cell_name.empty()) {
            for (const auto& load_layer : load_layers) {
                if (driver_layer != load_layer) {
                    return true;
                }
            }
        }
        
        // 检查负载之间是否有跨层连接
        if (load_layers.size() > 1) {
            std::string first_layer = load_layers[0];
            for (size_t i = 1; i < load_layers.size(); ++i) {
                if (first_layer != load_layers[i]) {
                    return true;
                }
            }
        }
        
        return false;
    }
};

} // namespace hbta
