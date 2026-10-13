#pragma once

#include "HBTRegion.h"
#include "HBTTypes.h"

#include <array>
#include <memory>
#include <vector>

namespace hbta {

bool compareHBTCandidatePointers(const HBTInfo* lhs, const HBTInfo* rhs);

struct HBTQuadTreeConfig {
    double min_region_wh_dbu = 50000.0;
};

struct HBTQuadTreeStats {
    int nodes = 0;
    int leaf_nodes = 0;
    int max_depth = 0;
    int max_node_candidates = 0;
    long long total_node_candidates = 0;
};

class HBTQuadTree {
public:
    struct Node {
        explicit Node(const HBTRegion& region_) : region(region_) {}

        HBTRegion region;
        std::vector<const HBTInfo*> candidate_hbts;
        std::array<std::unique_ptr<Node>, 4> children;
        bool has_children = false;
    };

    HBTQuadTree() = default;
    HBTQuadTree(const std::vector<const HBTInfo*>& hbts,
                const HBTQuadTreeConfig& config = HBTQuadTreeConfig{});

    bool build(const std::vector<const HBTInfo*>& hbts,
               const HBTQuadTreeConfig& config = HBTQuadTreeConfig{});

    bool empty() const { return root_ == nullptr; }
    const Node* root() const { return root_.get(); }
    const HBTRegion& bounds() const { return bounds_; }
    const HBTQuadTreeStats& stats() const { return stats_; }

private:
    std::unique_ptr<Node> buildNode(const HBTRegion& region,
                                    const std::vector<const HBTInfo*>& hbts,
                                    int depth);

    HBTQuadTreeConfig config_;
    HBTRegion bounds_;
    HBTQuadTreeStats stats_;
    std::unique_ptr<Node> root_;
};

} // namespace hbta
