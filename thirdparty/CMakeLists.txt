include(FetchContent)

#----------------------------------------
# Add malloc_count
#----------------------------------------
FetchContent_Declare(
  malloc_count
  GIT_REPOSITORY https://github.com/bingmann/malloc_count
  )
  
FetchContent_GetProperties(malloc_count)
if(NOT malloc_count_POPULATED)
  FetchContent_Populate(malloc_count)

  add_library(malloc_count OBJECT ${malloc_count_SOURCE_DIR}/malloc_count.c ${malloc_count_SOURCE_DIR}/malloc_count.h)
  target_link_libraries(malloc_count dl)
  target_include_directories(malloc_count PUBLIC "${malloc_count_SOURCE_DIR}")

  add_library(memprofile OBJECT ${malloc_count_SOURCE_DIR}/memprofile.h)
  target_include_directories(memprofile PUBLIC "${malloc_count_SOURCE_DIR}")
endif()

#----------------------------------------
# Add klib
#----------------------------------------
FetchContent_Declare(
  klib
  GIT_REPOSITORY https://github.com/attractivechaos/klib
)

FetchContent_GetProperties(klib)
if(NOT klib_POPULATED)
  FetchContent_Populate(klib)
  
  # add_subdirectory(${klib_SOURCE_DIR} ${klib_BINARY_DIR} EXCLUDE_FROM_ALL)
  add_library(klib INTERFACE)

  target_include_directories(klib INTERFACE ${klib_SOURCE_DIR})
endif()

#----------------------------------------
# Add zlib
#----------------------------------------

FetchContent_Declare(
  zlib
  GIT_REPOSITORY https://github.com/madler/zlib.git
)

FetchContent_GetProperties(zlib)
if(NOT zlib_POPULATED)
  message (STATUS "zlib is populated.")
  FetchContent_Populate(zlib)
  add_library(zlib INTERFACE)
  include_directories(${zlib_SOURCE_DIR})
  #target_include_directories(zlib INTERFACE ${zlib_SOURCE_DIR})
endif()

#----------------------------------------
# Add bonsai
#----------------------------------------
FetchContent_Declare(
  bonsai
  GIT_REPOSITORY https://github.com/dnbaker/bonsai.git
  GIT_TAG        5273b81a92
)

FetchContent_GetProperties(bonsai)
if(NOT bonsai_POPULATED)
  FetchContent_Populate(bonsai)  
  add_library(bonsai INTERFACE)
endif()

#----------------------------------------
# Add Big-BWT
#----------------------------------------
FetchContent_Declare(
  bigbwt
  GIT_REPOSITORY https://github.com/oma219/Big-BWT.git
  )
  
FetchContent_GetProperties(bigbwt)
if(NOT bigbwt_POPULATED)
  FetchContent_Populate(bigbwt)
  add_subdirectory(${bigbwt_SOURCE_DIR} ${bigbwt_BINARY_DIR})

  # add_library(gsacak OBJECT ${bigbwt_SOURCE_DIR}/gsa/gsacak.c ${bigbwt_SOURCE_DIR}/gsa/gsacak.h)
  # target_include_directories(gsacak PUBLIC "${bigbwt_SOURCE_DIR}/gsa")
  
  # add_library(malloc_count OBJECT ${bigbwt_SOURCE_DIR}/malloc_count.c ${bigbwt_SOURCE_DIR}/malloc_count.h)
  # target_link_libraries(malloc_count dl)
  # target_include_directories(malloc_count PUBLIC "${bigbwt_SOURCE_DIR}")
  
  endif()

#----------------------------------------
# Add gsacak
#----------------------------------------
FetchContent_Declare(
  gsacak
  GIT_REPOSITORY https://github.com/felipelouza/gsa-is.git
)

FetchContent_GetProperties(gsacak)
if(NOT gsacak_POPULATED)
    FetchContent_Populate(gsacak)
    add_library(gsacak OBJECT ${gsacak_SOURCE_DIR}/gsacak.c ${gsacak_SOURCE_DIR}/gsacak.h)
    target_include_directories(gsacak PUBLIC "${gsacak_SOURCE_DIR}")
    
    add_library(gsacak64 OBJECT ${gsacak_SOURCE_DIR}/gsacak.c ${gsacak_SOURCE_DIR}/gsacak.h)
    target_include_directories(gsacak64 PUBLIC "${gsacak_SOURCE_DIR}")
    target_compile_options(gsacak64 PUBLIC -DM64)
  # add_subdirectory(${gsacak_SOURCE_DIR} ${gsacak_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

#----------------------------------------
# Add sdsl
#----------------------------------------
FetchContent_Declare(
  sdsl
  # GIT_REPOSITORY https://github.com/elarielcl/sdsl-lite.git
  GIT_REPOSITORY https://github.com/simongog/sdsl-lite
)

FetchContent_GetProperties(sdsl)
if(NOT sdsl_POPULATED)
  FetchContent_Populate(sdsl)

  set(GENERATE_DOC OFF CACHE BOOL "Do not generate doxygen for sdsl-lite")
  
  add_subdirectory(${sdsl_SOURCE_DIR} ${sdsl_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

#----------------------------------------
# Add divsuffsort
#----------------------------------------
FetchContent_Declare(
  divsufsort
  GIT_REPOSITORY https://github.com/simongog/libdivsufsort.git
  GIT_TAG        2.0.1
)

FetchContent_GetProperties(divsufsort)
if(NOT divsufsort_POPULATED)
  FetchContent_Populate(divsufsort)

  set(BUILD_SHARED_LIBS OFF CACHE BOOL "Do not build a shared library for libdivsufsort")
  set(BUILD_EXAMPLES OFF CACHE BOOL "Do not build libdivsufsort example")
  set(BUILD_DIVSUFSORT64 ON CACHE BOOL "Build libdivsufsort in 64-bits mode")

  add_subdirectory(${divsufsort_SOURCE_DIR} ${divsufsort_BINARY_DIR} EXCLUDE_FROM_ALL)

  target_include_directories(divsufsort PUBLIC "${divsufsort_BINARY_DIR}/include")
  target_include_directories(divsufsort64 PUBLIC "${divsufsort_BINARY_DIR}/include")
endif()

#----------------------------------------
# Add r-index
#----------------------------------------
FetchContent_Declare(
  r-index
  GIT_REPOSITORY https://github.com/maxrossi91/r-index.git
)

FetchContent_GetProperties(r-index)
if(NOT r-index_POPULATED)
  FetchContent_Populate(r-index)

  add_subdirectory(${r-index_SOURCE_DIR} ${r-index_BINARY_DIR} )#EXCLUDE_FROM_ALL)
  add_library(ri INTERFACE)
  target_link_libraries(ri INTERFACE klib z)
  target_include_directories(ri INTERFACE ${r-index_SOURCE_DIR}/internal)
endif()

#----------------------------------------
# Add pfp-thresholds
#----------------------------------------
FetchContent_Declare(
  pfp_thresholds
  GIT_REPOSITORY https://github.com/oma219/pfp-thresholds.git
  GIT_TAG master
  )
  
FetchContent_GetProperties(pfp_thresholds)
if(NOT pfp_thresholds_POPULATED)
  FetchContent_Populate(pfp_thresholds)
  add_subdirectory(${pfp_thresholds_SOURCE_DIR} ${pfp_thresholds_BINARY_DIR})
  
endif()

message(STATUS "hello world ....")
message(STATUS $ENV{CC})

message(STATUS $CC)

#----------------------------------------
# Add bigrepair
#----------------------------------------
FetchContent_Declare(
  bigrepair
  GIT_REPOSITORY https://gitlab.com/omaryfekry/bigrepair.git
  )
  
FetchContent_GetProperties(bigrepair)
if(NOT bigrepair_POPULATED)
  FetchContent_Populate(bigrepair)
  # add_subdirectory(${bigrepair_SOURCE_DIR} ${bigrepair_BINARY_DIR})

  execute_process(COMMAND make CC=$ENV{CC} CXX=$ENV{CXX}
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${bigrepair_SOURCE_DIR} 
        COMMAND_ECHO STDOUT)
  if(result)
      message(FATAL_ERROR "CMake step for bigrepair failed: ${result}")
  endif()
  
endif()

#----------------------------------------
# Add ShapedSlp
#----------------------------------------
FetchContent_Declare(
  shaped_slp
  GIT_REPOSITORY https://github.com/koeppl/ShapedSlp.git
  GIT_TAG master
  )
  
FetchContent_GetProperties(shaped_slp)
if(NOT shaped_slp_POPULATED)
  FetchContent_Populate(shaped_slp)
  add_subdirectory(${shaped_slp_SOURCE_DIR} ${shaped_slp_BINARY_DIR})
  set(FOLCA_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/folca)
  set(SUX_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/external/sux/sux)
endif()

message(STATUS "finished the thirdparty code repos ...")


# ## Add pfp_ds
# FetchContent_Declare(
#   pfp_ds
#   GIT_REPOSITORY https://github.com/maxrossi91/pfp-data-structures
# )

# FetchContent_GetProperties(pfp_ds)
# if(NOT pfp_ds_POPULATED)
#   FetchContent_Populate(pfp_ds)

#   add_subdirectory(${pfp_ds_SOURCE_DIR} ${pfp_ds_BINARY_DIR} EXCLUDE_FROM_ALL)
# endif()

# ## Add rlbwt2lcp
# FetchContent_Declare(
#   rlbwt2lcp
#   GIT_REPOSITORY https://github.com/nicolaprezza/rlbwt2lcp
# )

# FetchContent_GetProperties(rlbwt2lcp)
# if(NOT rlbwt2lcp_POPULATED)
#   FetchContent_Populate(rlbwt2lcp)

#   add_subdirectory(${rlbwt2lcp_SOURCE_DIR} ${rlbwt2lcp_BINARY_DIR} )#EXCLUDE_FROM_ALL)
#   # add_library(rlbwt2lcp OBJECT  ${rlbwt2lcp_SOURCE_DIR}/internal/include.hpp 
#   #                               ${rlbwt2lcp_SOURCE_DIR}/internal/lcp.hpp 
#   #                               ${rlbwt2lcp_SOURCE_DIR}/internal/dna_string.hpp 
#   #                               ${rlbwt2lcp_SOURCE_DIR}/internal/dna_string_n.hpp
#   #                               ${rlbwt2lcp_SOURCE_DIR}/internal/dna_bwt.hpp
#   #                               ${rlbwt2lcp_SOURCE_DIR}/internal/dna_bwt_n.hpp
#   #                               )
#   # target_include_directories(rlbwt2lcp PUBLIC "${rlbwt2lcp_BINARY_DIR}/internal")
# endif()

# ## Add bwt2lcp
# FetchContent_Declare(
#   bwt2lcp
#   GIT_REPOSITORY https://github.com/nicolaprezza/bwt2lcp
# )

# FetchContent_GetProperties(bwt2lcp)
# if(NOT bwt2lcp_POPULATED)
#   FetchContent_Populate(bwt2lcp)

#   add_subdirectory(${bwt2lcp_SOURCE_DIR} ${bwt2lcp_BINARY_DIR} )#EXCLUDE_FROM_ALL)
#   # add_library(rlbwt2lcp OBJECT  ${rlbwt2lcp_SOURCE_DIR}/internal/include.hpp 
#   #                               ${rlbwt2lcp_SOURCE_DIR}/internal/lcp.hpp 
#   #                               ${rlbwt2lcp_SOURCE_DIR}/internal/dna_string.hpp 
#   #                               ${rlbwt2lcp_SOURCE_DIR}/internal/dna_string_n.hpp
#   #                               ${rlbwt2lcp_SOURCE_DIR}/internal/dna_bwt.hpp
#   #                               ${rlbwt2lcp_SOURCE_DIR}/internal/dna_bwt_n.hpp
#   #                               )
#   # target_include_directories(rlbwt2lcp PUBLIC "${rlbwt2lcp_BINARY_DIR}/internal")
# endif()

# ## Google benchmark
# FetchContent_Declare(
#   benchmark
#   GIT_REPOSITORY https://github.com/google/benchmark.git
#   GIT_TAG        master
#   )
  
# FetchContent_GetProperties(benchmark)
# if(NOT benchmark_POPULATED)
#   FetchContent_Populate(benchmark)
#   set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Enable testing of the benchmark library." FORCE)
#   add_subdirectory(${benchmark_SOURCE_DIR} ${benchmark_BINARY_DIR} EXCLUDE_FROM_ALL)
# endif()


  # ## Add tlx
  # FetchContent_Declare(
  #   tlx
  #   GIT_REPOSITORY https://github.com/tlx/tlx
  #   GIT_TAG master
  #   )
    
  # FetchContent_GetProperties(tlx)
  # if(NOT tlx_POPULATED)
  #   FetchContent_Populate(tlx)
  #   add_subdirectory(${tlx_SOURCE_DIR} ${tlx_BINARY_DIR})
  # endif()
  
  