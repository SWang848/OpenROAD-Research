// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#pragma once

#include "tcl.h"

namespace tdp {

class Timingdp;

// The TDP engine is owned by this module rather than ord::OpenRoad so the
// shared OpenRoad class and the existing placement engines remain unchanged.
Timingdp* getTimingdp();
void initTimingdp(Tcl_Interp* tcl_interp);

}  // namespace tdp
