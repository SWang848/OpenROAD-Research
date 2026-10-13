#pragma once

#include "HBTTypes.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <cmath>
#include <algorithm>
#include <functional>

// 为 std::pair<int, int> 提供哈希函数
namespace std {
    template<>
    struct hash<std::pair<int, int>> {
        size_t operator()(const std::pair<int, int>& p) const {
            // 使用简单的哈希组合
            return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
        }
    };
}

namespace hbta {

class HBTArray {
public:
    HBTArray();
    ~HBTArray() = default;
    
    // 从JSON文件加载HBT阵列
    bool loadFromJSON(const std::string& json_file, int dbu_per_micron = 1000);
    
    // 获取HBT信息
    const HBTInfo* getHBT(const std::string& id) const;
    const HBTInfo* getHBT(int row, int col) const;
    
    // 获取所有HBT
    const std::vector<HBTInfo>& getAllHBTs() const { return hbts_; }
    
    // 获取可用HBT数量
    int getAvailableCount() const;
    
    // 获取总HBT数量
    int getTotalCount() const { return hbts_.size(); }
    
    // 查找最近的可用HBT
    const HBTInfo* findNearestAvailableHBT(const Point& location) const;
    
    // 查找指定范围内的可用HBT
    std::vector<const HBTInfo*> findAvailableHBTsInRange(
        const Point& center, double max_distance) const;
    
    // 标记HBT为已使用
    bool markHBTUsed(const std::string& hbt_id, 
                     const std::string& connected_net = "");
    
    // 标记HBT为可用
    bool markHBTAvailable(const std::string& hbt_id);
    
    // 检查HBT是否可用
    bool isHBTAvailable(const std::string& hbt_id) const;
    
    // 获取阵列信息
    int getRows() const { return rows_; }
    int getColumns() const { return columns_; }
    double getPitchDBU() const { return pitch_dbu_; }
    
    // 重置所有HBT状态
    void resetAllHBTs();
    
    // 导出当前状态到JSON
    std::string exportToJSON() const;
    
private:
    std::vector<HBTInfo> hbts_;
    std::unordered_map<std::string, size_t> id_to_index_;
    std::unordered_map<std::pair<int, int>, size_t> coord_to_index_;
    
    int rows_ = 0;
    int columns_ = 0;
    double pitch_dbu_ = 0.0;
    
    // 辅助方法
    void buildIndexMaps();
    double calculateDistance(const Point& p1, const Point& p2) const {
        return p1.distanceTo(p2);
    }
    
    // JSON解析辅助方法
    bool parseJSONFile(const std::string& json_file, int dbu_per_micron);
    bool parseHBTArrayFromJSON(const std::string& json_content, int dbu_per_micron);
};

} // namespace hbta