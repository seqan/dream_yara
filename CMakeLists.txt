# ===========================================================================
#          DREAM-Yara - a distributed version of the Yara read mapper
# ===========================================================================

cmake_minimum_required (VERSION 3.8)

# Messages
string (ASCII 27 Esc)
set (FontBold "${Esc}[1m")
set (FontReset "${Esc}[m")
set (ColourRed "${Esc}[31m")

message (STATUS "${FontBold}Compiler Detection...${FontReset}")

# Define the application name and version.
project (DREAM-Yara CXX)

## BUILD

# ----------------------------------------------------------------------------
# Make "Release" the default cmake build type
# ----------------------------------------------------------------------------
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release Coverage RelWithDebInfo MinSizeRel."
         FORCE)
endif ()

# ----------------------------------------------------------------------------
# Warn if cmake build type is not "Release"
# ----------------------------------------------------------------------------
if (NOT CMAKE_BUILD_TYPE STREQUAL Release)
  message (STATUS "${ColourRed}CMAKE_BUILD_TYPE is not \"Release\", your binaries will be slow.${FontReset}")
endif ()

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# ----------------------------------------------------------------------------
# Begin of dependency detection
# ----------------------------------------------------------------------------
message ("\n${FontBold}Dependency detection${FontReset}")

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/lib/seqan/include/seqan/version.h")
    set (CMAKE_INCLUDE_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/lib/seqan/include
        ${CMAKE_INCLUDE_PATH})
    set (CMAKE_PREFIX_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/lib/seqan/util/cmake
        ${CMAKE_PREFIX_PATH})
    set (CMAKE_MODULE_PATH
       ${CMAKE_CURRENT_SOURCE_DIR}/lib/seqan/util/cmake
       ${CMAKE_MODULE_PATH})
    message (STATUS "Found a local SeqAn library provided with the ${CMAKE_PROJECT_NAME} source code.")
    message (    "   This will be preferred over system global headers.")
endif ()

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/lib/sdsl-lite/include/sdsl/bits.hpp")
    set (SDSL_INCLUDE_DIRS
        ${CMAKE_CURRENT_SOURCE_DIR}/lib/sdsl-lite/include)
    message (STATUS "Found a local SDSL library provided with the ${CMAKE_PROJECT_NAME} source code.")
endif ()

# Dependency: SeqAn3.
find_package (SeqAn3 QUIET REQUIRED HINTS lib/seqan3/build_system)

# Use ccache.
include ("${SEQAN3_CLONE_DIR}/test/cmake/seqan3_require_ccache.cmake")
seqan3_require_ccache ()

# ----------------------------------------------------------------------------
# Add the application.
# Add Yara targets
# ----------------------------------------------------------------------------
message (STATUS "${FontBold}Configuring ${CMAKE_PROJECT_NAME}...${FontReset}")
add_subdirectory (src)
message (STATUS "${FontBold}You can run `make` to build the application.${FontReset}")

## DOCUMENTATION

add_subdirectory (doc EXCLUDE_FROM_ALL)

# ----------------------------------------------------------------------------
## TEST
# Add Tests
# ----------------------------------------------------------------------------

message (STATUS "${FontBold}Configuring ${CMAKE_PROJECT_NAME} Tests...${FontReset}")
enable_testing ()
add_subdirectory (test EXCLUDE_FROM_ALL)
