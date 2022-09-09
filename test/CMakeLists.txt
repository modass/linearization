add_executable(
        gtest_test
        test_gtest.cpp
)
target_link_libraries(
        gtest_test
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(gtest_test)