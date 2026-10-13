include("openroad")

set(TEST_LIBS
  GTest::gtest
  GTest::gtest_main
  GTest::gmock
  odb
  dft_base_scan_cell_lib
  dft_config_lib
  dft_architect_lib
  dft_clock_domain_lib
  dft_utils_scan_pin_lib
)

add_executable(TestScanArchitect TestScanArchitect.cpp)
target_link_libraries(TestScanArchitect ${TEST_LIBS})
# gtest_discover_tests call removed

add_executable(TestScanArchitectHeuristic TestScanArchitectHeuristic.cpp ScanCellMock.cpp)
target_link_libraries(TestScanArchitectHeuristic ${TEST_LIBS})
# gtest_discover_tests call removed


add_dependencies(build_and_test
  TestScanArchitect
  TestScanArchitectHeuristic
)
