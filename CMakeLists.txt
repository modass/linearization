cmake_minimum_required(VERSION 3.20)
project(linearization)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(src/linearization)
add_subdirectory(src/tool)
add_subdirectory(dependencies)

enable_testing()
add_subdirectory(test)