add_executable(
        linearization_test
        test_heuristics.cpp
        test_relaxation.cpp
)
target_link_libraries(
        linearization_test
        linearization
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(linearization_test)
