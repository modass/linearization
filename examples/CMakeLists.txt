file(GLOB_RECURSE example_sources LIST_DIRECTORIES FALSE RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} example_*.cpp)
message(STATUS "Found example files: ${example_sources}")

foreach (file ${example_sources})
    string(REGEX REPLACE "\.cpp" "" target_name "${file}")
    add_executable(${target_name} ${file})
    target_link_libraries(${target_name} linearization spdlog CLI11)
endforeach ()
