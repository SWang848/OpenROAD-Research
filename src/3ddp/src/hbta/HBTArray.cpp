#include "HBTArray.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <limits>
#include <cstdlib>
#include <cctype>

// 简单的JSON解析（为了简化，我们使用简单解析）
// 在实际项目中应该使用nlohmann/json或其他JSON库

namespace hbta {

HBTArray::HBTArray() : rows_(0), columns_(0), pitch_dbu_(0.0) {}

bool HBTArray::loadFromJSON(const std::string& json_file, int dbu_per_micron) {
    // 清空现有数据
    hbts_.clear();
    id_to_index_.clear();
    coord_to_index_.clear();
    rows_ = 0;
    columns_ = 0;
    pitch_dbu_ = 0.0;
    
    // 读取JSON文件（简化解析，未使用第三方 JSON 库）
    std::ifstream file(json_file);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open HBT array file: " << json_file << std::endl;
        return false;
    }
    
    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string json_content = buffer.str();
    file.close();
    
    // 解析JSON内容 - 期望新格式包含: "hbt_array" 和 "hbt_positions" 字段。
    // 注：这里是轻量字符串解析，对字段顺序和格式较敏感。
    
    // 查找 hbt_array 字段（包含 rows/columns/pitch_um）
    size_t hbt_array_pos = json_content.find("\"hbt_array\"");
    if (hbt_array_pos == std::string::npos) {
        std::cerr << "Error: Cannot find 'hbt_array' field in JSON" << std::endl;
        return false;
    }
    
    // 查找 hbt_positions 字段（包含每个 HBT 的坐标与索引）
    size_t hbt_positions_pos = json_content.find("\"hbt_positions\"");
    if (hbt_positions_pos == std::string::npos) {
        std::cerr << "Error: Cannot find 'hbt_positions' field in JSON" << std::endl;
        return false;
    }
    
    // 从 hbt_array 中提取行数、列数和 pitch（微米）
    size_t rows_pos = json_content.find("\"rows\"", hbt_array_pos);
    size_t cols_pos = json_content.find("\"columns\"", hbt_array_pos);
    size_t pitch_pos = json_content.find("\"pitch_um\"", hbt_array_pos);
    
    if (rows_pos == std::string::npos || cols_pos == std::string::npos || pitch_pos == std::string::npos) {
        std::cerr << "Error: Missing required fields in hbt_array" << std::endl;
        return false;
    }
    
    // 提取行数（rows）
    rows_pos = json_content.find(":", rows_pos) + 1;
    rows_ = std::stoi(json_content.substr(rows_pos, json_content.find(",", rows_pos) - rows_pos));
    
    // 提取列数（columns）
    cols_pos = json_content.find(":", cols_pos) + 1;
    columns_ = std::stoi(json_content.substr(cols_pos, json_content.find(",", cols_pos) - cols_pos));
    
    // 提取 pitch（微米）并转换为 DBU
    pitch_pos = json_content.find(":", pitch_pos) + 1;
    double pitch_um = std::stod(json_content.substr(pitch_pos, json_content.find(",", pitch_pos) - pitch_pos));
    pitch_dbu_ = pitch_um * dbu_per_micron;
    
    // 提取 HBT 位置数组（hbt_positions）。
    hbt_positions_pos = json_content.find("[", hbt_positions_pos);
    size_t end_pos = json_content.find("]", hbt_positions_pos);
    
    if (hbt_positions_pos == std::string::npos || end_pos == std::string::npos) {
        std::cerr << "Error: Cannot find HBT positions array in JSON" << std::endl;
        return false;
    }
    
    std::string hbts_str = json_content.substr(hbt_positions_pos + 1, end_pos - hbt_positions_pos - 1);
    
    // 解析每个 HBT 条目：id/row/col/x_um/y_um。
    size_t start = 0;
    int hbt_count = 0;
    
    while (start < hbts_str.length()) {
        size_t obj_start = hbts_str.find("{", start);
        if (obj_start == std::string::npos) break;
        
        size_t obj_end = hbts_str.find("}", obj_start);
        if (obj_end == std::string::npos) break;
        
        std::string hbt_obj = hbts_str.substr(obj_start, obj_end - obj_start + 1);
        
        // 提取 HBT 信息 - 新格式使用 x_um 和 y_um
        size_t id_pos = hbt_obj.find("\"id\"");
        size_t row_pos = hbt_obj.find("\"row\"");
        size_t col_pos = hbt_obj.find("\"col\"");
        size_t x_pos = hbt_obj.find("\"x_um\"");
        size_t y_pos = hbt_obj.find("\"y_um\"");
        
        if (id_pos != std::string::npos && row_pos != std::string::npos && 
            col_pos != std::string::npos && x_pos != std::string::npos && 
            y_pos != std::string::npos) {
            
            // 提取 ID（稳健处理空格和引号，避免得到 `"hbt_0_0` 这种前导引号）
            id_pos = hbt_obj.find(":", id_pos) + 1;
            size_t id_end = hbt_obj.find(",", id_pos);
            if (id_end == std::string::npos) id_end = hbt_obj.find("}", id_pos);
            std::string id = hbt_obj.substr(id_pos, id_end - id_pos);

            auto trim_inplace = [](std::string& s) {
                while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) {
                    s.erase(s.begin());
                }
                while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) {
                    s.pop_back();
                }
            };
            trim_inplace(id);
            if (!id.empty() && id.front() == '"' && id.back() == '"' && id.size() >= 2) {
                id = id.substr(1, id.size() - 2);
            }
            
            // 提取 row
            row_pos = hbt_obj.find(":", row_pos) + 1;
            size_t row_end = hbt_obj.find(",", row_pos);
            if (row_end == std::string::npos) row_end = hbt_obj.find("}", row_pos);
            int row = std::stoi(hbt_obj.substr(row_pos, row_end - row_pos));
            
            // 提取 col
            col_pos = hbt_obj.find(":", col_pos) + 1;
            size_t col_end = hbt_obj.find(",", col_pos);
            if (col_end == std::string::npos) col_end = hbt_obj.find("}", col_pos);
            int col = std::stoi(hbt_obj.substr(col_pos, col_end - col_pos));
            
            // 提取 X 坐标（微米，字段名 x_um）
            x_pos = hbt_obj.find(":", x_pos) + 1;
            size_t x_end = hbt_obj.find(",", x_pos);
            if (x_end == std::string::npos) x_end = hbt_obj.find("}", x_pos);
            double x_um = std::stod(hbt_obj.substr(x_pos, x_end - x_pos));
            
            // 提取 Y 坐标（微米，字段名 y_um）
            y_pos = hbt_obj.find(":", y_pos) + 1;
            size_t y_end = hbt_obj.find(",", y_pos);
            if (y_end == std::string::npos) y_end = hbt_obj.find("}", y_pos);
            double y_um = std::stod(hbt_obj.substr(y_pos, y_end - y_pos));
            
            // 转换为 DBU（HBTArray 内部统一使用 DBU）
            double x_dbu = x_um * dbu_per_micron;
            double y_dbu = y_um * dbu_per_micron;
            
            // 内部统一存 DBU，避免后续和放置数据库单位混用。
            HBTInfo hbt(id, row, col, Point(x_dbu, y_dbu));
            hbts_.push_back(hbt);
            hbt_count++;
        }
        
        start = obj_end + 1;
    }
    
    // 构建索引映射（id->index, (row,col)->index）
    buildIndexMaps();
    
    std::cout << "Loaded HBT array: " << rows_ << "x" << columns_ 
              << " (" << hbt_count << " HBTs), pitch: " << pitch_dbu_ << " DBU" << std::endl;
    
    return true;
}

void HBTArray::buildIndexMaps() {
    id_to_index_.clear();
    coord_to_index_.clear();
    
    // 建双索引：按 id 访问与按(row,col)访问都走 O(1) 哈希。
    for (size_t i = 0; i < hbts_.size(); ++i) {
        const HBTInfo& hbt = hbts_[i];
        id_to_index_[hbt.id] = i;
        coord_to_index_[{hbt.row, hbt.col}] = i;
    }
}

const HBTInfo* HBTArray::getHBT(const std::string& id) const {
    auto it = id_to_index_.find(id);
    if (it != id_to_index_.end() && it->second < hbts_.size()) {
        return &hbts_[it->second];
    }
    return nullptr;
}

const HBTInfo* HBTArray::getHBT(int row, int col) const {
    auto it = coord_to_index_.find({row, col});
    if (it != coord_to_index_.end() && it->second < hbts_.size()) {
        return &hbts_[it->second];
    }
    return nullptr;
}

int HBTArray::getAvailableCount() const {
    int count = 0;
    for (const auto& hbt : hbts_) {
        if (!hbt.used) {
            count++;
        }
    }
    return count;
}

const HBTInfo* HBTArray::findNearestAvailableHBT(const Point& location) const {
    const HBTInfo* nearest = nullptr;
    double min_distance = std::numeric_limits<double>::max();
    
    for (const auto& hbt : hbts_) {
        if (!hbt.used) {
            double distance = calculateDistance(location, hbt.location);
            if (distance < min_distance) {
                min_distance = distance;
                nearest = &hbt;
            }
        }
    }
    
    return nearest;
}

std::vector<const HBTInfo*> HBTArray::findAvailableHBTsInRange(
    const Point& center, double max_distance) const {
    
    std::vector<const HBTInfo*> result;
    
    // 只收集未被占用且在给定半径内的 HBT
    for (const auto& hbt : hbts_) {
        if (!hbt.used) {
            double distance = calculateDistance(center, hbt.location);
            if (distance <= max_distance) {
                result.push_back(&hbt);
            }
        }
    }
    
    // 按与 center 的距离升序排序（最近的排在前面）
    std::sort(result.begin(), result.end(),
              [this, &center](const HBTInfo* a, const HBTInfo* b) {
                  return calculateDistance(center, a->location) < 
                         calculateDistance(center, b->location);
              });
    
    return result;
}

bool HBTArray::markHBTUsed(const std::string& hbt_id, 
                          const std::string& connected_net) {
    auto it = id_to_index_.find(hbt_id);
    if (it == id_to_index_.end() || it->second >= hbts_.size()) {
        return false;
    }
    
    HBTInfo& hbt = hbts_[it->second];
    if (hbt.used) {
        return false;  // HBT已经被使用
    }
    
    // 标记占用并记录连接的网络名称
    hbt.used = true;
    hbt.connected_net = connected_net;
    return true;
}

bool HBTArray::markHBTAvailable(const std::string& hbt_id) {
    auto it = id_to_index_.find(hbt_id);
    if (it == id_to_index_.end() || it->second >= hbts_.size()) {
        return false;
    }
    
    HBTInfo& hbt = hbts_[it->second];
    hbt.used = false;
    hbt.connected_net.clear();
    return true;
}

bool HBTArray::isHBTAvailable(const std::string& hbt_id) const {
    auto it = id_to_index_.find(hbt_id);
    if (it == id_to_index_.end() || it->second >= hbts_.size()) {
        return false;
    }
    
    return !hbts_[it->second].used;
}

void HBTArray::resetAllHBTs() {
    for (auto& hbt : hbts_) {
        hbt.used = false;
        hbt.connected_net.clear();
    }
}

std::string HBTArray::exportToJSON() const {
    std::stringstream ss;
    ss << "{\n";
    ss << "  \"rows\": " << rows_ << ",\n";
    ss << "  \"columns\": " << columns_ << ",\n";
    ss << "  \"pitch\": " << (pitch_dbu_ / 1000.0) << ",\n"; // 假设1000 DBU/微米
    ss << "  \"hbts\": [\n";
    
    for (size_t i = 0; i < hbts_.size(); ++i) {
        const HBTInfo& hbt = hbts_[i];
        ss << "    {\n";
        ss << "      \"id\": \"" << hbt.id << "\",\n";
        ss << "      \"row\": " << hbt.row << ",\n";
        ss << "      \"col\": " << hbt.col << ",\n";
        ss << "      \"x\": " << (hbt.location.x / 1000.0) << ",\n"; // 转换为微米
        ss << "      \"y\": " << (hbt.location.y / 1000.0) << ",\n";
        ss << "      \"used\": " << (hbt.used ? "true" : "false") << ",\n";
        ss << "      \"connected_net\": \"" << hbt.connected_net << "\"\n";
        ss << "    }";
        
        if (i < hbts_.size() - 1) {
            ss << ",";
        }
        ss << "\n";
    }
    
    ss << "  ]\n";
    ss << "}\n";
    
    return ss.str();
}

} // namespace hbta
