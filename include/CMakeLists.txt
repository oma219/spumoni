set(COMMON_SOURCES common.hpp)

add_library(common_h OBJECT ${COMMON_SOURCES})
target_link_libraries(common_h)

set(MS_SOURCES  ms_rle_string.hpp  thresholds_ds.hpp 
                spumoni_main.hpp compute_ms_pml.hpp
                doc_array.hpp refbuilder.hpp)

add_library(ms OBJECT ${MS_SOURCES})
target_link_libraries(ms common_h sdsl)
target_include_directories(ms PUBLIC "${r-index_SOURCE_DIR}/internal")
