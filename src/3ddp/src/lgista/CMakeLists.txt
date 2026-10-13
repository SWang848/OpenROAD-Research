# SPDX-License-Identifier: BSD-3-Clause
#
# local graph ISTA (lgista)
#
# A standalone local-graph incremental timing engine extracted from the lrf
# module. Contains ONLY the local timing graph (PtGraph) + local STA engine
# (LocalSta) + delay calc + parasitics + arrival/required search, plus a small
# write API to commit local results back to the global graph. The sizer
# (LrSizer), the Lagrangian-relaxation machinery, AND the scheduling/operator
# framework (TaskArranger / ParallelVisitor / operators) are intentionally
# excluded -- traversal order and optimization belong to whoever builds an
# optimizer on top of this engine.
#
# NOTE: this module uses its own `lgista` C++ namespace (it was extracted from
# the now-removed `lrf` module; some file-level comments still say "lrf").

include("openroad")

add_library(tdp_lgista_lib
    src/PtPiElmore.cc
    src/PtElmoreCeff.cc
    src/PtGraph.cc
    src/LocalReduceParasitic.cc
    src/LocalParasitics.cc
    src/LocalDmpDelayCalc.cc
    src/LocalElmoreCeffDelayCalc.cc
    src/LocalSearch.cc
    src/LocalPath.cc
    src/LocalSta.cc
    src/ViolationCheck.cc
    src/LrfSearchPred.cc
    src/IstaWriteback.cc
    src/LocalGraphIsta.cc
)

target_include_directories(tdp_lgista_lib
  PUBLIC
    include
  PRIVATE
    ${OPENSTA_HOME}
    ${OPENSTA_HOME}/include/sta
    src
)

target_link_libraries(tdp_lgista_lib
  PUBLIC
    odb
    OpenSTA
    dbSta_lib
    est_lib
    utl_lib
)

# ---------------------------------------------------------------------------
# Minimal C++ test driver: loads a design (.odb + .lib) and runs the local
# graph timer over its instances, printing local timing.  Built only when
# TDP_LGISTA_BUILD_TEST is enabled.  The standalone driver is retained from
# the source bundle for traceability, but is not part of the OpenROAD runtime.
# ---------------------------------------------------------------------------
option(TDP_LGISTA_BUILD_TEST "Build the private lgista standalone test program" OFF)
if(TDP_LGISTA_BUILD_TEST)
  add_executable(lgista_test
      src/LgistaMain.cc
      src/LgistaTests.cc
  )
  target_include_directories(lgista_test
    PRIVATE
      include
      ${OPENSTA_HOME}
      ${OPENSTA_HOME}/include/sta
      src
  )
  target_link_libraries(lgista_test
    PRIVATE
      tdp_lgista_lib
      dbSta_lib
      OpenSTA
      est_lib
      odb
      utl_lib
      ${TCL_LIBRARY}
      ${CMAKE_THREAD_LIBS_INIT}
  )
endif()
