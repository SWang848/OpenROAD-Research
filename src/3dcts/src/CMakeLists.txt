# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2019-2025, The OpenROAD Authors

include("openroad")

# https://github.com/The-OpenROAD-Project/OpenROAD/issues/1186
find_package(LEMON NAMES LEMON lemon REQUIRED)

# OR-Tools for LP solver (GLOP)
# Same pattern as gpl/mpl/par modules
find_package(ortools CONFIG REQUIRED)
message(STATUS "3DCTS: Found OR-Tools: ${ortools_DIR} (version: ${ortools_VERSION})")

add_library(dcts_lib
    Clock.cpp
    TreeBuilder.cpp
    HTreeBuilder.cpp
    SinkClustering.cpp
    TechChar.cpp
    TritonCTS.cpp
    Clustering.cpp
    LatencyBalancer.cpp
    CtsOptions.cpp
    VerilogFFExtractor.cpp
    Cts3DDatabase.cpp
    ClockLatencyEstimator.cpp
    CtsSkewLpSolver.cpp  
    BufSizingLpSolver.cpp
    TarpClustering.cpp   
)

swig_lib(NAME      dcts
         NAMESPACE dcts
         I_FILE    TritonCTS.i
         SCRIPTS   TritonCTS.tcl
)

target_sources(dcts
  PRIVATE
    Make3DCTS.cpp
    CtsGraphics.cpp
)

target_include_directories(dcts_lib
  PUBLIC
    ../include
  PRIVATE
    .
  ${LEMON_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
 )

target_include_directories(dcts
  PUBLIC
    ../include
  PRIVATE
    .
  ${Boost_INCLUDE_DIRS}
 )

target_link_libraries(dcts_lib
  PUBLIC
    rsz_lib
  PRIVATE
    odb
    dbSta_lib
    OpenSTA
    stt_lib
    utl_lib
    ortools::ortools
)

target_link_libraries(dcts
  PRIVATE
    odb
    dcts_lib
    dbSta_lib
    rsz_lib
    OpenSTA
    gui
 )

messages(
  TARGET dcts
  OUTPUT_DIR ..
)
