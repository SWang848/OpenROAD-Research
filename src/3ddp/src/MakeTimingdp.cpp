// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#include "tdp/MakeTimingdp.h"

#include "ord/OpenRoad.hh"
#include "tcl.h"
#include "tdp/Timingdp.h"
#include "utl/decode.h"

extern "C" {
extern int Tdp_Init(Tcl_Interp* interp);
}

namespace tdp {

// Tcl files encoded into strings.
extern const char* tdp_tcl_inits[];

namespace {

// OpenROAD tools have process lifetime.  Keeping TDP here avoids adding a
// plugin-specific member and getter to the shared ord::OpenRoad class.
Timingdp* timingdp = nullptr;

}  // namespace

Timingdp* getTimingdp()
{
  return timingdp;
}

void initTimingdp(Tcl_Interp* tcl_interp)
{
  ord::OpenRoad* openroad = ord::OpenRoad::openRoad();
  if (timingdp == nullptr) {
    timingdp = new Timingdp(openroad->getDb(),
                            openroad->getLogger(),
                            openroad->getSta(),
                            openroad->getResizer(),
                            openroad->getGlobalRouter(),
                            openroad->getOpendp());
  }

  // Define swig TCL commands.
  Tdp_Init(tcl_interp);
  // Eval encoded sta TCL sources.
  utl::evalTclInit(tcl_interp, tdp::tdp_tcl_inits);
}

}  // namespace tdp
