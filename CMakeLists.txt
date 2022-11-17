cmake_minimum_required(VERSION 3.20)
project(linearization VERSION 0.0.2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
option(COVERAGE "Enable to produce coverage reports" OFF)
option(BUILD_DOC "Build a doxygen documentation" OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (BUILD_DOC)
    add_subdirectory(doc)
endif ()

add_subdirectory(src/analysis)
add_subdirectory(src/linearization)
add_subdirectory(src/utility)
add_subdirectory(src/tool)
# add last, as the dependencies need to be linked to existing targets, which get defined in the previous includes
add_subdirectory(dependencies)

enable_testing()
add_subdirectory(test)

if (COVERAGE)
    message(STATUS "Producing test coverage report")
    include(CodeCoverage)
    append_coverage_compiler_flags()
    setup_target_for_coverage_lcov(NAME ${PROJECT_NAME}_coverage EXECUTABLE ctest)
endif ()
