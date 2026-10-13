#include "HBTACostModel.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace hbta {

namespace {

constexpr double kInfCost = 1e50;
struct Terminal {
    Point location;
    std::string layer;
};

double manhattanDistance(const Point& lhs, const Point& rhs) {
    return lhs.manhattanDistanceTo(rhs);
}

bool sameRoutingLayer(const std::string& lhs, const std::string& rhs) {
    return !lhs.empty() && !rhs.empty() && lhs == rhs;
}

std::vector<Terminal> collectTerminals(const InterLayerConnection& connection) {
    std::vector<Terminal> terminals;
    terminals.reserve((connection.driver_cell_name.empty() ? 0 : 1)
                      + connection.load_locations.size());

    if (!connection.driver_cell_name.empty()) {
        terminals.push_back({connection.driver_location, connection.driver_layer});
    }

    for (size_t i = 0; i < connection.load_locations.size(); ++i) {
        std::string layer;
        if (i < connection.load_layers.size()) {
            layer = connection.load_layers[i];
        }
        terminals.push_back({connection.load_locations[i], layer});
    }

    return terminals;
}

} // namespace

double HBTACostModel::computeRoutingAwareMSTCost(const InterLayerConnection& connection,
                                                 const Point& hbt_location) {
    const auto terminals = collectTerminals(connection);
    if (terminals.empty()) {
        return 0.0;
    }

    // 构图规则：
    // 1) 同层 pin 之间使用 Manhattan 距离作为边权；
    // 2) 跨层 pin 之间不建立直接边；
    // 3) HBT 伪节点与每个 pin 建边，边权为 Manhattan 距离。
    // 在该图上执行 Prim 算法，得到路由感知线长估计。
    const size_t hbt_index = terminals.size();
    const size_t node_count = terminals.size() + 1;
    std::vector<bool> selected(node_count, false);
    std::vector<double> best_edge(node_count, kInfCost);

    best_edge[0] = 0.0;
    double total_cost = 0.0;

    for (size_t iter = 0; iter < node_count; ++iter) {
        size_t best_node = node_count;
        for (size_t node = 0; node < node_count; ++node) {
            if (!selected[node]
                && (best_node == node_count || best_edge[node] < best_edge[best_node])) {
                best_node = node;
            }
        }

        if (best_node == node_count || best_edge[best_node] >= kInfCost) {
            // 若图不可连通，返回极大代价以表示该候选不可取。
            return kInfCost;
        }

        selected[best_node] = true;
        total_cost += best_edge[best_node];

        for (size_t to = 0; to < node_count; ++to) {
            if (selected[to]) {
                continue;
            }

            double edge_cost = kInfCost;
            if (best_node == hbt_index) {
                edge_cost = manhattanDistance(hbt_location, terminals[to].location);
            } else if (to == hbt_index) {
                edge_cost = manhattanDistance(terminals[best_node].location, hbt_location);
            } else if (sameRoutingLayer(terminals[best_node].layer, terminals[to].layer)) {
                edge_cost = manhattanDistance(terminals[best_node].location, terminals[to].location);
            }

            if (edge_cost < best_edge[to]) {
                best_edge[to] = edge_cost;
            }
        }
    }

    return total_cost;
}

double HBTACostModel::computeMatchingCost(const InterLayerConnection& connection,
                                          const Point& hbt_location) {
    return computeRoutingAwareMSTCost(connection, hbt_location);
}

double HBTACostModel::criticalityFromSlack(double slack,
                                           double worst_slack,
                                           double critical_threshold) {
    // criticality 越大，网络在重分配候选集合中的优先级越高。
    // 这里的 critical_threshold 直接作为“critical 裕量上界”使用。
    // 例如 threshold=10ps 时：
    //   slack >= 10ps -> criticality = 0
    //   slack = worst_slack -> criticality = 1
    //   中间线性插值。
    if (!std::isfinite(slack) || !std::isfinite(worst_slack)) {
        return 0.0;
    }

    if (slack >= critical_threshold) {
        return 0.0;
    }

    const double denominator = critical_threshold - worst_slack;
    if (denominator <= 0.0) {
        return 0.0;
    }

    const double criticality = (critical_threshold - slack) / denominator;
    return std::clamp(criticality, 0.0, 1.0);
}

} // namespace hbta
