#include "HBTQuadTree.h"

#include <algorithm>

namespace hbta {

bool compareHBTCandidatePointers(const HBTInfo* lhs, const HBTInfo* rhs) {
    if (lhs == rhs) {
        return false;
    }
    if (lhs == nullptr) {
        return false;
    }
    if (rhs == nullptr) {
        return true;
    }
    if (lhs->row != rhs->row) {
        return lhs->row < rhs->row;
    }
    if (lhs->col != rhs->col) {
        return lhs->col < rhs->col;
    }
    return lhs->id < rhs->id;
}

namespace {

std::vector<const HBTInfo*> normalizeHBTCandidateOrdering(
    const std::vector<const HBTInfo*>& hbts) {
    // 过滤空指针并按 row/col/id 排序
    std::vector<const HBTInfo*> ordered_hbts;
    ordered_hbts.reserve(hbts.size());
    for (const HBTInfo* hbt : hbts) {
        if (hbt != nullptr) {
            ordered_hbts.push_back(hbt);
        }
    }
    std::sort(ordered_hbts.begin(), ordered_hbts.end(), compareHBTCandidatePointers);
    return ordered_hbts;
}

HBTRegion computeHBTCandidateBounds(const std::vector<const HBTInfo*>& hbts) {
    if (hbts.empty()) {
        return HBTRegion();
    }

    // 根节点区域只包住候选 HBT 的中心点范围，后续子区域划分都基于这个 bbox。
    double min_x = hbts.front()->location.x;
    double max_x = hbts.front()->location.x;
    double min_y = hbts.front()->location.y;
    double max_y = hbts.front()->location.y;

    for (const HBTInfo* hbt : hbts) {
        min_x = std::min(min_x, hbt->location.x);
        max_x = std::max(max_x, hbt->location.x);
        min_y = std::min(min_y, hbt->location.y);
        max_y = std::max(max_y, hbt->location.y);
    }

    if (min_x == max_x) {
        min_x -= 0.5;
        max_x += 0.5;
    }
    if (min_y == max_y) {
        min_y -= 0.5;
        max_y += 0.5;
    }

    return HBTRegion(min_x, min_y, max_x, max_y);
}

bool shouldSplitNode(const HBTRegion& region, const HBTQuadTreeConfig& config) {
    return region.width() >= config.min_region_wh_dbu
        && region.height() >= config.min_region_wh_dbu;
}

int findPointChildIndex(const Point& location,
                        const std::array<HBTRegion, 4>& child_regions) {
    // splitHBTRegion 使用 nextafter 避免中线处被多个子区域同时包含；
    // 因此每个 HBT 中心最多只会落入一个 child。
    for (int i = 0; i < static_cast<int>(child_regions.size()); ++i) {
        if (child_regions[i].containsPoint(location)) {
            return i;
        }
    }
    return -1;
}

} // namespace

HBTQuadTree::HBTQuadTree(const std::vector<const HBTInfo*>& hbts,
                         const HBTQuadTreeConfig& config) {
    build(hbts, config);
}

bool HBTQuadTree::build(const std::vector<const HBTInfo*>& hbts,
                        const HBTQuadTreeConfig& config) {
    config_ = config;
    stats_ = HBTQuadTreeStats{};
    root_.reset();

    const auto ordered_hbts = normalizeHBTCandidateOrdering(hbts);
    if (ordered_hbts.empty()) {
        bounds_ = HBTRegion();
        return false;
    }

    bounds_ = computeHBTCandidateBounds(ordered_hbts);
    if (!bounds_.isValid()) {
        return false;
    }

    root_ = buildNode(bounds_, ordered_hbts, 0);
    return root_ != nullptr;
}

std::unique_ptr<HBTQuadTree::Node> HBTQuadTree::buildNode(
    const HBTRegion& region,
                                    const std::vector<const HBTInfo*>& hbts,
                                    int depth) {
    auto node = std::make_unique<Node>(region);
    node->candidate_hbts = hbts;

    stats_.nodes++;
    stats_.max_depth = std::max(stats_.max_depth, depth);
    stats_.max_node_candidates = std::max(
        stats_.max_node_candidates,
        static_cast<int>(node->candidate_hbts.size()));
    stats_.total_node_candidates +=
        static_cast<long long>(node->candidate_hbts.size());

    if (node->candidate_hbts.empty() || !shouldSplitNode(region, config_)) {
        stats_.leaf_nodes++;
        return node;
    }

    const auto child_regions = splitHBTRegion(region);
    std::array<std::vector<const HBTInfo*>, 4> child_hbts;
    for (const HBTInfo* hbt : node->candidate_hbts) {
        // 每个候选 HBT 只被下推到包含其中心点的子区域；不复制到多个 child。
        const int child_index = findPointChildIndex(hbt->location, child_regions);
        if (child_index >= 0) {
            child_hbts[child_index].push_back(hbt);
        }
    }

    node->has_children = true;
    for (size_t i = 0; i < child_regions.size(); ++i) {
        // 即使某个 child 没有 HBT，也保留空节点，递归求解时可直接判空/不可行。
        node->children[i] = buildNode(child_regions[i], child_hbts[i], depth + 1);
    }

    return node;
}

} // namespace hbta
