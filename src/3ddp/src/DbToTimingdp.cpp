// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#include "tdp/Timingdp.h"
#include "infrastructure/Architecture3D.h"
#include "infrastructure/Network3D.h"
#include "infrastructure/Common3D.h"

#include <algorithm>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "odb/db.h"
#include "odb/dbTypes.h"
#include "odb/geom.h"
#include "utl/Logger.h"

namespace tdp {

using utl::TDP;
// ============================================================================
// 数据库转换实现（仿照DPL的dbToOpendp.cpp）
// ============================================================================

bool Timingdp::designset(){
    
    block_ = db_->getChip()->getBlock();
    if (!db_ || !block_) {
        logger_->error(TDP, 10000, "Database not ready for import");
        return false;
    }

    importClear();
    
    // 2. 获取核心区域
    core_ = block_->getCoreArea();
    if (core_.isInverted()) {
        logger_->error(TDP, 3000, "Core area is inverted");
        return false;
    }
    
    return true;
}

bool Timingdp::importDb()
{
    // block_ = db_->getChip()->getBlock();
    // if (!db_ || !block_) {
    //     logger_->error(TDP, 1, "Database not ready for import");
    //     return false;
    // }
    
    try {
        logger_->info(TDP, 2, "Importing database to 3D representation...");
        
        // // 1. 清除现有数据
        // importClear();
        
        // // 2. 获取核心区域
        // core_ = block_->getCoreArea();
        // if (core_.isInverted()) {
        //     logger_->error(TDP, 3, "Core area is inverted");
        //     return false;
        // }
        
        // // 3. 设置层架构（仿照DPL的createArchitecture）
        // if (!setupLayerArchitecture()) {
        //     //logger_->error("Failed to setup layer architecture");
        //     return false;
        // }
        
        // // 4. 转换行（仿照DPL的dbToOpendp::createArchitecture）
        // if (!convertRows()) {
        //     //logger_->error("Failed to convert rows");
        //     return false;
        // }
        
        // 5. 转换实例（仿照DPL的dbToOpendp::createNetwork）
        if (!convertInstances()) {
            //logger_->error("Failed to convert instances");
            return false;
        }
        
        // 6. 转换网络（仿照DPL的dbToOpendp::createNetwork）
        if (!convertNets()) {
            //logger_->error("Failed to convert nets");
            return false;
        }
        
        // 7. 设置布局组（仿照DPL的setUpPlacementGroups）
        if (!setUpPlacementGroups()) {
            //logger_->warn("Failed to setup placement groups, continuing...");
        }
        
        // 8. 提取PDK参数（仿照DPL从数据库获取）：此阶段仅刷新间距。
        if (!extractPDKParametersFromDatabase("")) {
            //logger_->warn("Failed to extract PDK parameters, using defaults...");
        }
        
        // 9. 缓存设计数据
        cacheDesignData();

        //初始化段管理器
        initializeSegmentManager();
        
        // 10. 标记为已转换
        is_converted_ = true;
        is_loaded_ = true;
        
        logger_->info(TDP, 400, "Successfully imported database: " +
                      std::to_string(all_insts_.size()) + " instances, " +
                      std::to_string(nets_.size()) + " nets");
        return true;
        
    } catch (const std::exception& e) {
        //logger_->error("Failed to import database: " + std::string(e.what()));
        return false;
    }
}

void Timingdp::importClear()
{
    //logger_->info("Clearing 3D conversion data...");
    
    // 清除3D架构数据（类似DPL的importClear）
    if (arch_3d_) {
        arch_3d_->clear();
    }
    
    // 清除3D网络数据
    if (network_3d_) {
        network_3d_->clear();
    }
    
    // 清除缓存数据
    cleanup();
    
    // 清除核心区域
    core_ = odb::Rect();
    
    // 重置转换状态
    is_converted_ = false;
    is_loaded_ = false;
    
    //logger_->info("3D conversion data cleared");
}

// ============================================================================
// 行转换方法（仿照DPL的dbToOpendp::createArchitecture）
// ============================================================================

bool Timingdp::convertRows(const std::string& layer)
{
    if (!block_ || !arch_3d_) {
        logger_->error(TDP, 401, "Block or architecture not available for row conversion");
        return false;
    }
    
    try {
        //logger_->info("Converting rows to 3D architecture...");
        
        // 计算最小行高（仿照DPL处理混合工艺）
        DBU min_row_height = std::numeric_limits<DBU>::max();
        for (odb::dbRow* row : block_->getRows()) {
            if (row->getSite()) {
                DBU row_height = static_cast<DBU>(row->getSite()->getHeight());
                min_row_height = std::min(min_row_height, row_height);
            }
        }
        
        // 跳过列表（处理多高单元）
        std::map<DBU, std::unordered_set<std::string>> skip_list;
        
        int row_count = 0;
        for (odb::dbRow* db_row : block_->getRows()) {
            if (!db_row->getSite()) {
                //logger_->warn("Row has no site, skipping");
                continue;
            }
            
            // 跳过PAD行（仿照DPL）
            if (db_row->getSite()->getClass() == odb::dbSiteClass::PAD) {
                continue;
            }
            
            // 只处理水平行（仿照DPL）
            if (db_row->getDirection() != odb::dbRowDir::HORIZONTAL) {
                continue;
            }
            
            odb::dbSite* site = db_row->getSite();
            
            // 跳过高度大于最小行高的site（多高单元处理，仿照DPL）
            if (site->getHeight() > min_row_height) {
                skip_list[site->getHeight()].insert(site->getName());
                continue;
            }
            
            // 层类型由调用方（newrow 的"建一层 site -> 抓一层"顺序）显式给定，
            // 避免同构（上下层同名 site）下靠 lib 名/高度猜测产生的误判。
            Architecture3D::LayerType layer_type = Architecture3D::stringToLayer(layer);
            
            // 创建3D行（仿照DPL的createRowFromDatabase）
            Architecture3D::Row3D* row_3d = arch_3d_->createAndAddRow(layer_type);
            if (!row_3d) {
                //logger_->warn("Failed to create 3D row, skipping");
                continue;
            }
            
            // 设置行ID
            row_3d->setId(row_count++);
            
            // 设置行参数（仿照DPL的dbToOpendp）
            odb::Point origin = db_row->getOrigin();
            // row_3d->setBottom(DbuY{origin.y() - core_.yMin()});
            // row_3d->setSubRowOrigin(DbuX{origin.x() - core_.xMin()});
            row_3d->setBottom(DbuY{origin.y()});
            row_3d->setSubRowOrigin(DbuX{origin.x()});
            row_3d->setSiteSpacing(DbuX{db_row->getSpacing()});
            row_3d->setNumSites(db_row->getSiteCount());
            row_3d->setSiteWidth(DbuX{site->getWidth()});
            row_3d->setHeight(DbuY{site->getHeight()});
            row_3d->setOrient(db_row->getOrient());
            
            // 设置对称性（从site获取，仿照DPL）
            unsigned symmetry = 0;
            if (site->getSymmetryX()) symmetry |= 0x01;
            if (site->getSymmetryY()) symmetry |= 0x02;
            if (site->getSymmetryR90()) symmetry |= 0x04;
            row_3d->setSymmetry(symmetry);
            
            //logger_->debugCheck("Converted row " + std::to_string(row_count) +
            //              " to layer " + std::to_string(layer_type) +
            //              ", sites=" + std::to_string(db_row->getSiteCount()));
        }
        
        // 设置边界（仿照DPL）
        for (int layer = 0; layer < Architecture3D::LAYER_COUNT; layer++) {
            Architecture3D::LayerType layer_type = static_cast<Architecture3D::LayerType>(layer);
            
            if (arch_3d_->hasLayer(layer_type)) {
                DbuX xmin = DbuX::max();
                DbuX xmax = DbuX::lowest();
                DbuY ymin = DbuY::max();
                DbuY ymax = DbuY::lowest();
                
                for (int r = 0; r < arch_3d_->getNumRows(layer_type); r++) {
                    auto* row = arch_3d_->getRow(layer_type, r);
                    if (!row) continue;
                    
                    xmin = std::min(xmin, row->getLeft());
                    xmax = std::max(xmax, row->getRight());
                    ymin = std::min(ymin, row->getBottom());
                    ymax = std::max(ymax, row->getTop());
                }
                
                arch_3d_->setMinX(layer_type, xmin);
                arch_3d_->setMaxX(layer_type, xmax);
                arch_3d_->setMinY(layer_type, ymin);
                arch_3d_->setMaxY(layer_type, ymax);
            }
        }
        
        //logger_->info("Successfully converted " + std::to_string(row_count) + " rows");
        return true;
        
    } catch (const std::exception& e) {
        //logger_->error("Failed to convert rows: " + std::string(e.what()));
        return false;
    }
}

// ============================================================================
// 实例转换方法（仿照DPL的dbToOpendp::createNetwork）
// ============================================================================

bool Timingdp::convertInstances()
{
    if (!block_ || !network_3d_ || !arch_3d_) {
        //logger_->error("Block, network or architecture not available for instance conversion");
        return false;
    }
    
    try {
        logger_->info(TDP, 10, "Converting instances to 3D network...");
        
        // 按名称排序（仿照DPL）
        auto block_insts = block_->getInsts();
        std::vector<odb::dbInst*> insts(block_insts.begin(), block_insts.end());
        std::sort(insts.begin(), insts.end(), 
                  [](odb::dbInst* a, odb::dbInst* b) {
                      return a->getName() < b->getName();
                  });
        
        // 缓存每个 master 的层类型，避免对每个 instance 重复做字符串转换
        std::unordered_map<odb::dbMaster*, Architecture3D::LayerType> master_layer_cache;

        int instance_count = 0;
        for (odb::dbInst* db_inst : insts) {
            // 跳过不可放置的实例（仿照DPL）
            if (!db_inst->getMaster()->isCoreAutoPlaceable()) {
                continue;
            }
            
            // 确定层类型（通过 master 缓存）：统一依据标准单元名后缀。
            odb::dbMaster* master = db_inst->getMaster();
            auto layer_it = master_layer_cache.find(master);
            if (layer_it == master_layer_cache.end()) {
                Architecture3D::LayerType layer_type = Architecture3D::layerOfMaster(master);
                layer_it = master_layer_cache.emplace(master, layer_type).first;
            }
            Architecture3D::LayerType layer_type = layer_it->second;
            
            // 添加主单元到网络（仿照DPL的addMaster）
            network_3d_->addMaster(master, arch_3d_.get());
            
            // 添加实例到网络（仿照DPL的addNode）
            network_3d_->addNode(db_inst);
            
            // 获取添加的节点并设置层信息
            Node3D* node = network_3d_->getCell(db_inst);
            if (node) {
                node->layer = (layer_type == Architecture3D::LAYER_BOTTOM) ? "bottom" : "upper";
            } else {
                //logger_->warn("Failed to get node after adding: " + std::string(db_inst->getName()));
                continue;
            }
            
            instance_count++;
            
            //logger_->debugCheck("Converted instance " + std::string(db_inst->getName()) +
            //              " to layer " + std::to_string(layer_type));
        }
        
        // ============================================
        // 处理IO Ports (BTerm) - 参考DPL
        // ============================================
        int terminal_count = 0;
        for (odb::dbBTerm* bterm : block_->getBTerms()) {
            // 跳过电源网络
            odb::dbNet* net = bterm->getNet();
            if (!net || net->getSigType().isSupply()) {
                continue;
            }
            // 检查BTerm是否有有效的边界框
            if (bterm->getBBox().isInverted()) {
                //logger_->warn("BTerm " + std::string(bterm->getName()) + " has no shapes, skipping");
                continue;
            }
            // 添加BTerm为Terminal节点（层归属由 pin 几何相对 hb_layer 判定）
            network_3d_->addNode(bterm);
            terminal_count++;
        }
        //logger_->info("Successfully converted " + std::to_string(terminal_count) + " IO terminals");

        // 批量添加完成后一次性重新分类 movable/fixed
        network_3d_->updateCellCategories();

        //logger_->info("Successfully converted " + std::to_string(instance_count) + " instances");
        return true;
        
    } catch (const std::exception& e) {
        //logger_->error("Failed to convert instances: " + std::string(e.what()));
        return false;
    }
}

// ============================================================================
// 网络转换方法（仿照DPL的dbToOpendp::createNetwork）
// ============================================================================

bool Timingdp::convertNets()
{
    if (!block_ || !network_3d_) {
        //logger_->error("Block or network not available for net conversion");
        return false;
    }
    
    try {
        logger_->info(TDP, 11, "Converting nets to 3D network...");
        
        int net_count = 0;
        for (odb::dbNet* db_net : block_->getNets()) {
            // 跳过电源网络（仿照DPL）
            if (db_net->getSigType().isSupply()) {
                continue;
            }
            
            // 添加网络到3D网络（仿照DPL的addEdge）
            try {
                network_3d_->addEdge(db_net);
            } catch (const std::exception& e) {
                //logger_->warn("Failed to add net to network: " + std::string(db_net->getName()) + " - " + e.what());
                continue;
            }
            
            // 获取刚刚创建的Net3D对象
            Net3D* net_3d = network_3d_->getNet(db_net);
            if (!net_3d) {
                //logger_->warn("Failed to get Net3D object for net: " + std::string(db_net->getName()));
                continue;
            }
            
            // 填充connected_cells和pins：遍历网络连接的所有实例引脚
            for (odb::dbITerm* iterm : db_net->getITerms()) {
                odb::dbInst* inst = iterm->getInst();
                if (!inst) continue;
                
                // 获取对应的Node3D对象
                Node3D* node = network_3d_->getCell(inst);
                if (node) {
                    // 添加到网络的connected_cells
                    net_3d->connected_cells.push_back(node);
                    // 反向映射：node → db_net
                    node->addConnectedNet(db_net);
                    
                    // 创建Pin并连接到网络
                    Pin3D* pin = network_3d_->addPin(iterm);
                    if (pin) {
                        // 连接Pin到Node
                        network_3d_->connect(pin, node);
                        // 连接Pin到Net
                        network_3d_->connect(pin, net_3d);
                    }
                }
            }
            
            // 填充connected_cells和pins：遍历网络连接的所有终端
            for (odb::dbBTerm* bterm : db_net->getBTerms()) {
                // 获取对应的Node3D对象
                Node3D* node = network_3d_->getNode(bterm);
                if (node) {
                    // 添加到网络的connected_cells
                    net_3d->connected_cells.push_back(node);
                    // 反向映射：node → db_net
                    node->addConnectedNet(db_net);
                    
                    // 创建Pin并连接到网络
                    Pin3D* pin = network_3d_->addPin(bterm);
                    if (pin) {
                        // 连接Pin到Node
                        network_3d_->connect(pin, node);
                        // 连接Pin到Net
                        network_3d_->connect(pin, net_3d);
                    }
                }
            }
            
            net_count++;
            
            // logger_->debugCheck("Converted net " + std::string(db_net->getName()) + 
            //          " with " + std::to_string(net_3d->connected_cells.size()) + " connected cells, " +
            //          std::to_string(net_3d->getPins().size()) + " pins");
        }
        
        //logger_->info("Successfully converted " + std::to_string(net_count) + " nets");
        return true;
        
    } catch (const std::exception& e) {
        //logger_->error("Failed to convert nets: " + std::string(e.what()));
        return false;
    }
}

// ============================================================================
// 布局组设置方法（仿照DPL的setUpPlacementGroups）
// ============================================================================

bool Timingdp::setUpPlacementGroups()
{
    if (!block_ || !arch_3d_ || !network_3d_) {
        //logger_->error("Block, architecture or network not available for group setup");
        return false;
    }
    
    try {
        logger_->info(TDP, 12, "Setting up placement groups...");
        
        int group_count = 0;
        auto db_groups = block_->getGroups();
        
        for (auto db_group : db_groups) {
            odb::dbRegion* region = db_group->getRegion();
            if (!region) {
                continue;
            }
            
            // 为每个层创建组（3D扩展）
            for (int layer = 0; layer < Architecture3D::LAYER_COUNT; layer++) {
                Architecture3D::LayerType layer_type = static_cast<Architecture3D::LayerType>(layer);
                
            // 创建3D组
            auto* group_3d_ptr = arch_3d_->createAndAddRegion(layer_type);
            if (!group_3d_ptr) {
                continue;
            }
            
            // 转换为Group3D指针
            Group3D* group_3d = static_cast<Group3D*>(group_3d_ptr);
            group_3d->setId(group_count++);
            
            // 设置组边界（仿照DPL）
            odb::Rect bbox;
            bbox.mergeInit();
            for (odb::dbBox* boundary : region->getBoundaries()) {
                odb::Rect box = boundary->getBox();
                box = box.intersect(core_);
                box.moveDelta(-core_.xMin(), -core_.yMin());
                group_3d->addRect(box);
                bbox.merge(box);
            }
            group_3d->setBoundary(bbox);
            
            // 添加实例到组（仿照DPL）
            for (auto db_inst : db_group->getInsts()) {
                auto* node = network_3d_->getCell(db_inst);
                if (node) {
                    // 检查实例是否属于当前层
                    Architecture3D::LayerType inst_layer = network_3d_->getCellLayer(db_inst);
                    if (inst_layer == layer_type) {
                        node->setGroupId(group_3d->getId());
                        node->setGroup(group_3d);
                        group_3d->addCell(node);
                    }
                }
            }
            
            // logger_->debugCheck("Created group " + std::to_string(group_3d->getId()) +
            //               " for layer " + std::to_string(layer_type));
            }
        }
        
        //logger_->info("Successfully setup " + std::to_string(group_count) + " placement groups");
        return true;
        
    } catch (const std::exception& e) {
        //logger_->error("Failed to setup placement groups: " + std::string(e.what()));
        return false;
    }
}

// ============================================================================
// 状态检查方法
// ============================================================================
// void Timingdp::adjustNodesOrient()
// {
//   for (auto& node : network_3d_->getNode()) {
//     if (node->getType() != Node::CELL) {
//       continue;
//     }
//     auto inst = node->getDbInst();
//     node->adjustCurrOrient(inst->getOrient());
//   }
// }


bool Timingdp::isDatabaseReady() const
{
    return db_ != nullptr && block_ != nullptr;
}

bool Timingdp::isConverted() const
{
    return is_converted_;
}

} // namespace tdp