include(ExternalProject)

find_package(hypro CONFIG REQUIRED)
if (${hypro_FOUND})
    message(STATUS "HyPro found")
    target_link_libraries(linearization hypro-shared)
endif ()

include(mc++.cmake)