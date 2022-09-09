include(ExternalProject)
include(FetchContent)

add_custom_target(dependencies)

find_package(hypro CONFIG REQUIRED)
if (${hypro_FOUND})
    message(STATUS "HyPro found")
    target_link_libraries(linearization PUBLIC hypro-shared)
endif ()

include(mc++.cmake)
include(gtest.cmake)