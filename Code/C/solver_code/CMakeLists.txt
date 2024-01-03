# CMake configuration for ECOS

cmake_minimum_required(VERSION 3.5)

project(ecos
        VERSION 2.0.8
        LANGUAGES C)

# Options
option(USE_LONG "Whether to use Long or Int for index type." ON)

# Set compiler flags
set(DEBUG_OPTIONS "")
set(RELEASE_OPTIONS "")

if(MSVC)
    list(APPEND DEBUG_OPTIONS "/Od" "/Wall")
    list(APPEND RELEASE_OPTIONS "/O2")
else()
    list(APPEND DEBUG_OPTIONS "-O0" "-Wall" "-Wextra")
    list(APPEND RELEASE_OPTIONS "-O3")
endif()

# Add the ECOS headers
set(ecos_headers
    "${CMAKE_CURRENT_SOURCE_DIR}/external/SuiteSparse_config/SuiteSparse_config.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/cone.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/ctrlc.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/data.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/ecos.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/ecos_bb.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/equil.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/expcone.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/glblopts.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/kkt.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/spla.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/splamm.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/timer.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/wright_omega.h"
)

set(ecos_sources
    # AMD
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_1.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_2.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_aat.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_control.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_defaults.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_dump.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_global.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_info.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_order.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_post_tree.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_postorder.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_preprocess.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_valid.c"

    # LDL
    "${CMAKE_CURRENT_SOURCE_DIR}/external/ldl/src/ldl.c"

    # ECOS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/cone.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/ctrlc.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/ecos.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/equil.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/expcone.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/kkt.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/preproc.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/spla.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/splamm.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/timer.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/wright_omega.c"

    # ECOS BB
    "${CMAKE_CURRENT_SOURCE_DIR}/ecos_bb/ecos_bb.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/ecos_bb/ecos_bb_preproc.c"
)


set(solver_head "${ecos_headers}" PARENT_SCOPE)
set(solver_src "${ecos_sources}" PARENT_SCOPE)