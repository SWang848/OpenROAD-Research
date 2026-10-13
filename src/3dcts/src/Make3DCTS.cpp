// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#include "dcts/Make3DCTS.h"

#include "dcts/TritonCTS.h"
#include "ord/OpenRoad.hh"
#include "tcl.h"
#include "utl/decode.h"

extern "C" {
extern int Dcts_Init(Tcl_Interp* interp);
}

namespace dcts {

extern const char* dcts_tcl_inits[];

namespace {

// OpenROAD owns process-wide modules for the lifetime of the executable.  Keep
// the isolated 3DCTS instance at the same lifetime without adding a member to
// ord::OpenRoad or changing its destructor.
TritonCTS* triton3d_cts = nullptr;

}  // namespace

TritonCTS* getTriton3DCts()
{
  return triton3d_cts;
}

void init3DCTS(Tcl_Interp* tcl_interp)
{
  ord::OpenRoad* openroad = ord::OpenRoad::openRoad();
  if (triton3d_cts == nullptr) {
    triton3d_cts = new TritonCTS(openroad->getLogger(),
                                 openroad->getDb(),
                                 openroad->getDbNetwork(),
                                 openroad->getSta(),
                                 openroad->getSteinerTreeBuilder(),
                                 openroad->getResizer(),
                                 openroad->getEstimateParasitics());
  }

  // Define swig TCL commands.
  Dcts_Init(tcl_interp);
  utl::evalTclInit(tcl_interp, dcts::dcts_tcl_inits);
}

}  // namespace dcts
