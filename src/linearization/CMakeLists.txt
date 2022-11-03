add_library(linearization SHARED library.cpp)
target_include_directories(linearization PUBLIC include)
target_link_libraries(linearization PUBLIC utility spdlog)
