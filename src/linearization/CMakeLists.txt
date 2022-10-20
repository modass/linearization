add_library(linearization SHARED library.cpp mcCormick.cpp include/linearization/mcCormick.h)
target_include_directories(linearization PUBLIC include)
target_link_libraries(linearization PUBLIC utility)