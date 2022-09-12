cmake_minimum_required(VERSION 3.20)
project(linearization)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
option(COVERAGE "Enable to produce coverage reports" OFF)

add_subdirectory(src/linearization)
add_subdirectory(src/tool)
add_subdirectory(dependencies)

enable_testing()
add_subdirectory(test)

if(COVERAGE)
    message(STATUS "Producing test coverage report")
    include(CodeCoverage)
    setup_target_for_coverage_gcovr_xml(NAME ${PROJECT_NAME}_coverage EXECUTABLE gtest_test)
endif()