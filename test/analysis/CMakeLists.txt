add_executable(
        analysis_test
        test_analysis.cpp
)
target_link_libraries(
        analysis_test
        analysis
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(analysis_test)
