# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025-2025, The OpenROAD Authors

include("openroad")
find_package(TCL)
find_package(Boost CONFIG)
find_package(LEMON NAMES LEMON lemon REQUIRED)
find_package(OpenMP REQUIRED)

# lgista is a private timing-analysis dependency of TDP.  Keep it below the
# 3ddp source tree instead of registering another top-level OpenROAD module.
add_subdirectory(src/lgista)

swig_lib(NAME         tdp
         NAMESPACE    tdp
         I_FILE       src/Timingdp.i
         SCRIPTS      src/Timingdp.tcl
         SWIG_INCLUDES ${ODB_HOME}/src/swig/common
                       ${ODB_HOME}/src/swig/tcl
                       ${ODB_HOME}/include
)

# tdp_lib is intended to be the functional code of this module with
# without GUI dependencies.  C++ unit test will depend on tdp_lib
# to minimize pulling the many transitive dependencies of the GUI.
add_library(tdp_lib
  src/DbToTimingdp.cpp
  src/Timingdp.cpp
  
  src/infrastructure/Architecture3D.cpp
  src/infrastructure/Master3D.cpp
  src/infrastructure/Net3D.cpp
  src/infrastructure/Network3D.cpp
  src/infrastructure/Node3D.cpp
  src/infrastructure/Pin3D.cpp
  src/infrastructure/Segment3D.cpp
  src/infrastructure/SegmentManager3D.cpp

  src/eva/eva.cpp
  src/eva/TimingEngine3D.cpp
  src/placement/placement.cpp
  src/placement/color.cxx
  src/placement/LSP3D.cpp
  src/placement/MIS3D.cpp
  
  src/hbta/HBTACostModel.cpp
  src/hbta/HBTAllocator.cpp
  src/hbta/HBTArray.cpp
  src/hbta/HBTDatabaseInterface.cpp
  src/hbta/HBTQuadTree.cpp
  src/hbta/HBTQuadTreeAllocator.cpp
)

target_link_libraries(tdp_lib
  PUBLIC
    odb
    dbSta
    grt_lib
    rsz_lib
    tdp_lgista_lib
    dpl_lib
    ifp
  PRIVATE
    utl_lib
    OpenSTA
    est_lib
    stt_lib
    OpenMP::OpenMP_CXX
    
)

target_include_directories(tdp_lib
  PUBLIC
    include
    ${LEMON_INCLUDE_DIRS}
  PRIVATE
    src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/lgista/src
    ${OPENSTA_HOME}
    ${OPENSTA_HOME}/include/sta
)

# tdp is the library that includes the GUI dependent code.  OpenROAD
# dependes on tdp.
target_sources(tdp
  PRIVATE
    src/MakeTimingdp.cpp
)

target_include_directories(tdp
  PUBLIC
    include
  PRIVATE
    src
)

target_link_libraries(tdp
  PRIVATE
    ${TCL_LIBRARY}
    tdp_lib
    gui
)

messages(
  TARGET tdp
)

################################################################

if (Python3_FOUND AND BUILD_PYTHON)
  swig_lib(NAME          tdp_py
           NAMESPACE     tdp
           LANGUAGE      python
           I_FILE        src/Timingdp-py.i
           SWIG_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/include/tdp
                         ${ODB_HOME}/src/swig/common
                         ${ODB_HOME}/src/swig/python
           SCRIPTS       ${CMAKE_CURRENT_BINARY_DIR}/tdp_py.py
  )

  target_link_libraries(tdp_py
    PUBLIC
      tdp_lib
      gui
      odb
  )

endif()

if(ENABLE_TESTS)
  # add_subdirectory(test)
endif()
