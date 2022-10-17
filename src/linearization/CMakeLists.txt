add_library(linearization SHARED library.cpp mcpp/mcCormick.cpp include/mcpp/mcCormick.h)
target_include_directories(linearization PUBLIC include)