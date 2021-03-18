set(MS_SOURCES  ms_rle_string.hpp
ms_pointers.hpp
thresholds_ds.hpp)

add_library(ms OBJECT ${MS_SOURCES})
target_link_libraries(ms common sdsl)
target_include_directories(ms PUBLIC "${r-index_SOURCE_DIR}/internal")
