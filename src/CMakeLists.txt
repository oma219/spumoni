#-------------------------------------------------------------------
# Builds the executable for compressing the dictionary
#-------------------------------------------------------------------
add_executable(compress_dictionary compress_dictionary.cpp)
target_link_libraries(compress_dictionary common_h sdsl malloc_count)
target_include_directories(compress_dictionary PUBLIC  "../include")

FetchContent_GetProperties(r-index)
FetchContent_GetProperties(shaped_slp)
FetchContent_GetProperties(klib)
FetchContent_GetProperties(bigbwt)

set(FOLCA_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/folca)
set(SUX_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/external/sux/sux)

#-------------------------------------------------------------------
# Builds the final executable 
#-------------------------------------------------------------------

## Declare some variables
set(bonsai_SOURCE_DIR "${CMAKE_BINARY_DIR}/_deps/bonsai-src")
set(zlib_SOURCE_DIR "${CMAKE_BINARY_DIR}/_deps/zlib-src")
set(klib_SOURCE_DIR "${CMAKE_BINARY_DIR}/_deps/klib-src")

set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(spumoni spumoni.cpp  compute_ms_pml.cpp doc_array.cpp 
                        refbuilder.cpp emp_null_database.cpp 
                        ks_test.cpp batch_loader.cpp)
target_link_libraries(spumoni sdsl common_h divsufsort divsufsort64 ri pthread zlib bonsai "-fopenmp")
target_include_directories(spumoni PUBLIC
                            "../include"
                            "${shaped_slp_SOURCE_DIR}" "${FOLCA_SOURCE_DIR}" "${SUX_SOURCE_DIR}/function" 
                            "${SUX_SOURCE_DIR}/support" 
                            "${bonsai_SOURCE_DIR}/include/bonsai"
                            "${bonsai_SOURCE_DIR}"
                            "${bonsai_SOURCE_DIR}/hll/include/"
                            "${bonsai_SOURCE_DIR}/include")
target_compile_options(spumoni PUBLIC "-std=c++17" "-fopenmp")

