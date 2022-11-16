ExternalProject_Add(
        mc++_ep
        LOG_CONFIGURE 1
        LOG_DOWNLOAD 1
        LOG_BUILD 1
        GIT_REPOSITORY https://github.com/coin-or/MCpp
        GIT_SHALLOW 1
        UPDATE_COMMAND ""
        PATCH_COMMAND git apply ${CMAKE_SOURCE_DIR}/dependencies/mcpp.patch
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} -C ${CMAKE_BINARY_DIR}/dependencies/mc++_ep-prefix/src/mc++_ep/MCpp/src install
)
ExternalProject_Get_Property(mc++_ep source_dir)
ExternalProject_Get_Property(mc++_ep binary_dir)

#message(STATUS "MC++ source dir: ${source_dir}")
#message(STATUS "MC++ binary dir: ${binary_dir}")

set(mc++_INCLUDE_DIR "${source_dir}")
#message(STATUS "MC++ include dir: ${mc++_INCLUDE_DIR}")

target_include_directories(linearization PUBLIC $<BUILD_INTERFACE:${mc++_INCLUDE_DIR}>)

add_dependencies(linearization mc++_ep)