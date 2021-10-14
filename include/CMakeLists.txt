set(COMMON_SOURCES common.hpp)

add_library(common_h OBJECT ${COMMON_SOURCES})
target_link_libraries(common_h)

set(MS_SOURCES  ms_rle_string.hpp ms_pointers.hpp thresholds_ds.hpp spumoni.hpp spumoni_main.hpp
                run_spumoni.hpp matching_statistics.hpp compute_ms_pml.hpp)

add_library(ms OBJECT ${MS_SOURCES})
target_link_libraries(ms common_h sdsl)
target_include_directories(ms PUBLIC "${r-index_SOURCE_DIR}/internal")
