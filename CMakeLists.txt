CMAKE_MINIMUM_REQUIRED(VERSION 3.15)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")

# Set a default build type if none was specified
# ------------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

# About this project
# ------------------------------------------------------------------------------
project(moni)
SET(VERSION_MAJOR "1")
SET(VERSION_MINOR "0")
SET(VERSION_PATCH "0")
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

# Set environment
# ------------------------------------------------------------------------------
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install)

find_package(Git)
if(GIT_FOUND)
    message("git found: ${GIT_EXECUTABLE}")
else()
    message(WARNING "git not found. Cloning of submodules will not work.")
endif()

# Configure thirdparty
# ------------------------------------------------------------------------------
set(CMAKE_INSTALL_INCLUDEDIR "include") # This is an hack because include(GUIInstallDirs) doesn't work
add_subdirectory(thirdparty)

# Configure the compiler with the appropriate flags
# ------------------------------------------------------------------------------
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  # using Clang
  include(ConfigureCompilerClang)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  include(ConfigureCompilerGcc)
else ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	message(FATAL_ERROR "Only the compiler gcc and clang are supported")
endif()

# Adding different sub-directories to the CMake build
#-------------------------------------------------------------------------------
add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(utils)

#configure_file(${PROJECT_SOURCE_DIR}/pipeline/moni ${PROJECT_BINARY_DIR}/moni)

# Move BigRepair binaries to be used by pipeline
configure_file(${PROJECT_BINARY_DIR}/_deps/bigrepair-src/procdic ${PROJECT_BINARY_DIR}/bin/procdic COPYONLY)
configure_file(${PROJECT_BINARY_DIR}/_deps/bigrepair-src/postproc ${PROJECT_BINARY_DIR}/bin/postproc COPYONLY)
configure_file(${PROJECT_BINARY_DIR}/_deps/bigrepair-src/largeb_repair/irepair ${PROJECT_BINARY_DIR}/bin/irepair COPYONLY)

# Install pipeline executables
install(TARGETS spumoni DESTINATION ${PROJECT_BINARY_DIR})
install(TARGETS newscanNT.x pscan.x newscan.x pfp_thresholds64 pfp_thresholds 
                compress_dictionary SlpEncBuild rlebwt_ms_build build_spumoni 
        DESTINATION ${PROJECT_BINARY_DIR}/bin)