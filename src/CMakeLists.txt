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
# Builds the two executables for index building
#-------------------------------------------------------------------
# Commented out construction of executables that build index
#add_executable(rlebwt_ms_build rlebwt_ms_build.cpp)
#target_link_libraries(rlebwt_ms_build common_h sdsl divsufsort divsufsort64 ri)
#target_include_directories(rlebwt_ms_build PUBLIC    "../include" 
#                                        "${shaped_slp_SOURCE_DIR}" 
#                                        "${FOLCA_SOURCE_DIR}" 
#                                        "${SUX_SOURCE_DIR}/function" 
#                                        "${SUX_SOURCE_DIR}/support"
#                                        )
#target_compile_options(rlebwt_ms_build PUBLIC "-std=c++17")

#add_executable(build_spumoni build_spumoni.cpp)
#target_link_libraries(build_spumoni common_h sdsl divsufsort divsufsort64 malloc_count ri)
#target_include_directories(build_spumoni PUBLIC    "../include")
#target_compile_options(build_spumoni PUBLIC "-std=c++17")

# Commented out construction of stand-alone executables
#add_executable(ms matching_statistics.cpp ${bigbwt_SOURCE_DIR}/xerrors.c)
#target_link_libraries(ms common_h sdsl divsufsort divsufsort64 malloc_count ri pthread)
#target_include_directories(ms PUBLIC    "../include"  
#                                        "${shaped_slp_SOURCE_DIR}" 
#                                        "${FOLCA_SOURCE_DIR}" 
#                                        "${SUX_SOURCE_DIR}/function" 
#                                        "${SUX_SOURCE_DIR}/support"
#                                        "${bigbwt_SOURCE_DIR}"
#                                        )
#target_compile_options(ms PUBLIC "-std=c++17")

#add_executable(run_spumoni run_spumoni.cpp ${bigbwt_SOURCE_DIR}/xerrors.c)
#target_link_libraries(run_spumoni common sdsl divsufsort divsufsort64 malloc_count ri pthread)
#target_include_directories(run_spumoni PUBLIC    "../../include/ms" "../../include/common" "${bigbwt_SOURCE_DIR}" "${klib_SOURCE_DIR}")
#target_compile_options(run_spumoni PUBLIC "-std=c++17")

#-------------------------------------------------------------------
# Builds the final executable 
#-------------------------------------------------------------------

## Declare some variables
set(bonsai_SOURCE_DIR "${CMAKE_BINARY_DIR}/_deps/bonsai-src")
set(zlib_SOURCE_DIR "${CMAKE_BINARY_DIR}/_deps/zlib-src")
set(klib_SOURCE_DIR "${CMAKE_BINARY_DIR}/_deps/klib-src")
 #"${bigbwt_SOURCE_DIR}" "${bigbwt_SOURCE_DIR}"

# "${bigbwt_SOURCE_DIR}"
# ${bigbwt_SOURCE_DIR}/xerrors.c
# "${bigbwt_SOURCE_DIR}"

add_executable(spumoni spumoni.cpp  compute_ms_pml.cpp doc_array.cpp refbuilder.cpp emp_null_database.cpp ks_test.cpp)
target_link_libraries(spumoni sdsl common_h divsufsort divsufsort64 ri pthread zlib bonsai "-fsanitize=address")
target_include_directories(spumoni PUBLIC 
                            "../include" 
                            "${shaped_slp_SOURCE_DIR}" "${FOLCA_SOURCE_DIR}" "${SUX_SOURCE_DIR}/function" 
                            "${SUX_SOURCE_DIR}/support" 
                            "${bonsai_SOURCE_DIR}/include/bonsai"
                            "${bonsai_SOURCE_DIR}"
                            "${bonsai_SOURCE_DIR}/hll/include/"
                            "${bonsai_SOURCE_DIR}/include")
target_compile_options(spumoni PUBLIC "-std=c++17" "-fsanitize=address")

#set(bonsai_SOURCE_DIR "${CMAKE_BINARY_DIR}_deps/bonsai")
#message(STATUS "include_directory_is = ${bonsai_SOURCE_DIR}/include/bonsai")
#message(STATUS "include_directory_is = ${bonsai_BINARY_DIR}/include/bonsai")
#message(STATUS "include_directory_is = ${SUX_SOURCE_DIR}/support")
#message(STATUS "include_directory_is = ${shaped_slp_SOURCE_DIR}")
#message(STATUS "include_directory_is = ${CMAKE_BINARY_DIR}")

# Old binary tree code ...
#add_executable(bst bst.cpp)
#target_include_directories(bst PUBLIC "../include")

