##############################################
# CMake build script for the BLASR library
##############################################

cmake_policy(SET CMP0048 NEW)
project(BLASR VERSION 5.3.0 LANGUAGES CXX C)
cmake_minimum_required(VERSION 3.2)

set(ROOT_PROJECT_NAME ${PROJECT_NAME} CACHE STRING "root project name")

# Build type
IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile RelWithDebInfo ReleaseWithAssert" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Build-time options
option(BLASR_build_bin    "Build binaries." ON)
option(BLASR_build_tests  "Build BLASR's unit tests." ON)
option(BLASR_use_ccache   "Build BLASR using ccache, if available." ON)

# Main project paths
set(BLASR_RootDir       ${BLASR_SOURCE_DIR})

# Project configuration
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake ${CMAKE_MODULE_PATH})

# Fixed order, do not sort or shuffle
include(blasr-ccache)
include(blasr-releasewithassert)
include(blasr-dependencies)
include(blasr-compilerflags)
include(blasr-gitsha1)
include(blasr-config)

file(WRITE libcpp/pbdata/LibBlasrConfig.h "")

include_directories(
    ${BLASR_RootDir}/libcpp/hdf
    ${BLASR_RootDir}/libcpp/alignment
    ${BLASR_RootDir}/libcpp/pbdata
    ${HDF5_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${LIBLZMA_INCLUDE_DIRS}
    ${BZIP2_INCLUDE_DIRS}
    ${CURL_INCLUDE_DIRS}
    ${CRYPTO_INCLUDE_DIRS}
    ${PacBioBAM_INCLUDE_DIRS}
    ${CMAKE_BINARY_DIR}/generated
)

# Dependencies
file(GLOB HDF5_CPP  "${BLASR_RootDir}/libcpp/hdf/*.cpp")
file(GLOB_RECURSE ALIGNMENT_CPP  "${BLASR_RootDir}/libcpp/alignment/*.cpp")
file(GLOB_RECURSE PBDATA_CPP  "${BLASR_RootDir}/libcpp/pbdata/*.cpp")

add_library(deps STATIC
    ${HDF5_CPP}
    ${ALIGNMENT_CPP}
    ${PBDATA_CPP}
)

target_link_libraries(
    deps
    ${PacBioBAM_LIBRARIES}
    ${hts_LIBRARIES}
    ${ZLIB_LDFLAGS}
    ${HDF5_LDFLAGS}
    ${LIBLZMA_LDFLAGS}
    ${BZIP2_LDFLAGS}
    ${CURL_LDFLAGS}
    ${CRYPTO_LDFLAGS}
    )

# Binaries
function(create_exe exe_name exe_cpp)
    add_executable(${exe_name} ${exe_cpp})

    set_target_properties(${exe_name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    )

    target_link_libraries(${exe_name}
        deps
        ${CMAKE_DL_LIBS}
        Threads::Threads
    )
    if (ARGV2)
        target_link_libraries(${exe_name}
            hdf5
            hdf5_cpp
        )
        GET_PROPERTY(HDF5_LINKER_FLAG_LOCAL GLOBAL PROPERTY HDF5_LINKER_FLAG_GLOBAL)
        if (HDF5_LINKER_FLAG_LOCAL)
            set_target_properties(${exe_name} PROPERTIES LINK_FLAGS ${HDF5_LINKER_FLAG_LOCAL})
        endif()
    endif()

    install(TARGETS ${exe_name} RUNTIME DESTINATION bin)
endfunction()

## BLASR
create_exe(blasr ${BLASR_RootDir}/Blasr.cpp T)

## utils
create_exe(sawriter ${BLASR_RootDir}/utils/SAWriter.cpp)
create_exe(toAfg ${BLASR_RootDir}/utils/ToAfg.cpp T)
