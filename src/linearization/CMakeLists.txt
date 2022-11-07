add_library(linearization SHARED uniformRandom.cpp)
target_include_directories(linearization PUBLIC include)
target_link_libraries(linearization PUBLIC utility spdlog)
