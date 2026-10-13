#pragma once

#include "HBTTypes.h"

namespace hbta {

class HBTACostModel {
public:
    // 路由感知 Prim/MST 线长估计：
    static double computeRoutingAwareMSTCost(const InterLayerConnection& connection,
                                             const Point& hbt_location);

    // 基础匹配代价等于 routing-aware physical cost。
    static double computeMatchingCost(const InterLayerConnection& connection,
                                      const Point& hbt_location);

    // 根据 slack 与最差 slack 的相对比例计算时序关键度。
    static double criticalityFromSlack(double slack,
                                       double worst_slack,
                                       double critical_threshold);
};

} // namespace hbta
