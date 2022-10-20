include(ExternalProject)
include(FetchContent)

add_custom_target(dependencies)

find_package(hypro CONFIG REQUIRED)
if (${hypro_FOUND})
    message(STATUS "HyPro found")
    target_link_libraries(linearization PUBLIC hypro-shared)
    target_link_libraries(analysis INTERFACE hypro-shared)
    target_include_directories(analysis INTERFACE ${hypro_INCLUDE_DIR})
endif ()

include(mc++.cmake)
include(gtest.cmake)