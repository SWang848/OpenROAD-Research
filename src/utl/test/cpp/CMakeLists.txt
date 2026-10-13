include("openroad")

set(TEST_LIBS
  GTest::gtest
  GTest::gtest_main
  GTest::gmock
  utl_lib
)

add_executable(TestCFileUtils TestCFileUtils.cpp)

target_link_libraries(TestCFileUtils ${TEST_LIBS})

# gtest_discover_tests call removed

add_dependencies(build_and_test
  TestCFileUtils
)


add_executable(TestMetrics TestMetrics.cpp)

target_link_libraries(TestMetrics ${TEST_LIBS})

# gtest_discover_tests call removed

add_dependencies(build_and_test
  TestMetrics
)

add_executable(TestAlgorithms TestAlgorithms.cpp)

target_link_libraries(TestAlgorithms ${TEST_LIBS})

# gtest_discover_tests call removed

add_dependencies(build_and_test
  TestAlgorithms
)

add_executable(TestSuppressStdout TestSuppressStdout.cpp)

target_link_libraries(TestSuppressStdout ${TEST_LIBS})

# gtest_discover_tests call removed

add_dependencies(build_and_test
  TestSuppressStdout
)
