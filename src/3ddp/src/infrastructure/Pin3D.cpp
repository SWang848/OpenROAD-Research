// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#include "Pin3D.h"
#include "Node3D.h"
#include "Net3D.h"
#include "Common3D.h"

#include <algorithm>

namespace tdp {

Pin3D::Pin3D() {
    // 默认构造函数
}

Pin3D::~Pin3D() {
    // 析构函数
}

void Pin3D::setDirection(int dir) {
    dir_ = dir;
}

int Pin3D::getDirection() const {
    return dir_;
}

void Pin3D::setNode(Node3D* node) {
    node_ = node;
}

void Pin3D::setNet(Net3D* net) {
    net_ = net;
}

Node3D* Pin3D::getNode() const {
    return node_;
}

Net3D* Pin3D::getNet() const {
    return net_;
}

void Pin3D::setOffsetX(DbuX offsetX) {
    offsetX_ = offsetX;
}

DbuX Pin3D::getOffsetX() const {
    return offsetX_;
}

void Pin3D::setOffsetY(DbuY offsetY) {
    offsetY_ = offsetY;
}

DbuY Pin3D::getOffsetY() const {
    return offsetY_;
}

void Pin3D::setOffsetZ(DbuZ offsetZ) {
    offsetZ_ = offsetZ;
}

DbuZ Pin3D::getOffsetZ() const {
    return offsetZ_;
}

void Pin3D::setPinLayer(int layer) {
    pinLayer_ = layer;
}

int Pin3D::getPinLayer() const {
    return pinLayer_;
}

void Pin3D::setPinWidth(DbuX width) {
    pinWidth_ = width;
}

DbuX Pin3D::getPinWidth() const {
    return pinWidth_;
}

void Pin3D::setPinHeight(DbuY height) {
    pinHeight_ = height;
}

DbuY Pin3D::getPinHeight() const {
    return pinHeight_;
}

// 权重管理
void Pin3D::setEarlyWeight(double weight) {
    early_weight_ = weight;
}

double Pin3D::getEarlyWeight() const {
    return early_weight_;
}

void Pin3D::setLateWeight(double weight) {
    late_weight_ = weight;
}

double Pin3D::getLateWeight() const {
    return late_weight_;
}

void Pin3D::setEarlyWeightini(double weight) {
    early_weight_ini_ = weight;
}

double Pin3D::getEarlyWeightini() const {
    return early_weight_ini_;
}

void Pin3D::setLateWeightini(double weight) {
    late_weight_ini_ = weight;
}

double Pin3D::getLateWeightini() const {
    return late_weight_ini_;
}

// 时序信息管理（可选）
void Pin3D::setEarlySlack(double slack) {
    early_slack_ = slack;
}

double Pin3D::getEarlySlack() const {
    return early_slack_;
}

void Pin3D::setLateSlack(double slack) {
    late_slack_ = slack;
}

double Pin3D::getLateSlack() const {
    return late_slack_;
}

void Pin3D::setEarlyArrival(double arrival) {
    early_arrival_ = arrival;
}

double Pin3D::getEarlyArrival() const {
    return early_arrival_;
}

void Pin3D::setLateArrival(double arrival) {
    late_arrival_ = arrival;
}

double Pin3D::getLateArrival() const {
    return late_arrival_;
}

void Pin3D::setEarlyRequired(double required) {
    early_required_ = required;
}

double Pin3D::getEarlyRequired() const {
    return early_required_;
}

void Pin3D::setLateRequired(double required) {
    late_required_ = required;
}

double Pin3D::getLateRequired() const {
    return late_required_;
}

Pointdp3D Pin3D::getLocation3D() const {
    if (!node_) {
        return Pointdp3D(DbuX{0}, DbuY{0}, DbuZ{0});
    }
    
    // // 计算引脚相对于单元中心的位置
    // DbuX cell_center_x = node_->location.x + node_->getWidth() / 2;
    // DbuY cell_center_y = node_->location.y + node_->getHeight() / 2;
    // DbuZ cell_z = node_->location.z;
    
    // // 加上偏移量
    // DbuX pin_x = cell_center_x + offsetX_;
    // DbuY pin_y = cell_center_y + offsetY_;
    // DbuZ pin_z = cell_z + offsetZ_;

    DbuX pin_x = node_->location.x + offsetX_;
    DbuY pin_y = node_->location.y + offsetY_;
    DbuZ pin_z = node_->location.z + offsetZ_;
    
    return Pointdp3D(pin_x, pin_y, pin_z);
}

void Pin3D::setDbITerm(odb::dbITerm* iterm) {
    db_iterm_ = iterm;
    
    if (iterm) {
        // 从ITerm获取偏移信息
        odb::dbMTerm* mterm = iterm->getMTerm();
        if (mterm) {
            // 获取引脚几何信息
            // 注意：odb::dbMPin没有getBoxes()方法，使用getGeometries()或直接访问
            for (odb::dbMPin* mpin : mterm->getMPins()) {
                // 遍历几何形状
                odb::dbSet<odb::dbBox> boxes = mpin->getGeometry();
                for (odb::dbBox* box : boxes) {
                    // 计算相对于单元中心的偏移
                    int offset_x = box->xMin() + (box->xMax() - box->xMin()) / 2;
                    int offset_y = box->yMin() + (box->yMax() - box->yMin()) / 2;
                    
                    
                    offsetX_ = DbuX{offset_x};
                    offsetY_ = DbuY{offset_y};
                    break;
                }
                break;
            }
            
            // 设置引脚尺寸
            pinWidth_ = DbuX{0};  // 需要从实际几何获取
            pinHeight_ = DbuY{0}; // 需要从实际几何获取
            
            // 设置方向
            odb::dbIoType io_type = mterm->getIoType();
            switch (io_type.getValue()) {
                case odb::dbIoType::INPUT:
                    dir_ = Dir_IN;
                    break;
                case odb::dbIoType::OUTPUT:
                    dir_ = Dir_OUT;
                    break;
                case odb::dbIoType::INOUT:
                    dir_ = Dir_INOUT;
                    break;
                default:
                    dir_ = Dir_UNKNOWN;
                    break;
            }
        }
    }
}

void Pin3D::setDbBTerm(odb::dbBTerm* bterm) {
    db_bterm_ = bterm;
    
    if (bterm) {
        // 从BTerm获取位置信息
        odb::Rect bbox = bterm->getBBox();
        if (!bbox.isInverted()) {
            // 设置引脚位置（相对于芯片边界）
            offsetX_ = DbuX{bbox.xMin() + (bbox.xMax() - bbox.xMin()) / 2};
            offsetY_ = DbuY{bbox.yMin() + (bbox.yMax() - bbox.yMin()) / 2};
            offsetZ_ = DbuZ{0};  // 终端通常在边界上
            
            // 设置引脚尺寸
            pinWidth_ = DbuX{bbox.xMax() - bbox.xMin()};
            pinHeight_ = DbuY{bbox.yMax() - bbox.yMin()};
        }
        
        // 设置方向
        odb::dbIoType io_type = bterm->getIoType();
        switch (io_type.getValue()) {
            case odb::dbIoType::INPUT:
                dir_ = Dir_IN;
                break;
            case odb::dbIoType::OUTPUT:
                dir_ = Dir_OUT;
                break;
            case odb::dbIoType::INOUT:
                dir_ = Dir_INOUT;
                break;
            default:
                dir_ = Dir_UNKNOWN;
                break;
        }
    }
}

odb::dbITerm* Pin3D::getDbITerm() const {
    return db_iterm_;
}

odb::dbBTerm* Pin3D::getDbBTerm() const {
    return db_bterm_;
}

std::string Pin3D::getName() const {
    return name_;
}

void Pin3D::setName(const std::string& name) {
    name_ = name;
}

} // namespace tdp