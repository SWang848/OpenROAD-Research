// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#pragma once

#include "tcl.h"

namespace dcts {

class TritonCTS;

// The 3DCTS engine is intentionally owned by this module instead of
// ord::OpenRoad so the original CTS object and public OpenRoad state remain
// unchanged.
TritonCTS* getTriton3DCts();
void init3DCTS(Tcl_Interp* tcl_interp);

}  // namespace dcts
