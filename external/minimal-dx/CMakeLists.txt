# TEB version 4.0.0 (https://github.com/teb-model/teb).
# Copyright 2018 D. Meyer. Licensed under CeCILL version 2.1.

include(ExternalProject)

ExternalProject_Add(minimaldx-project
    GIT_REPOSITORY https://github.com/dmey/minimal-dx.git
    GIT_TAG master
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DUSE_REAL8=${USE_REAL8}
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
)

ExternalProject_Get_Property(minimaldx-project install_dir)

add_library(minimaldx INTERFACE)
add_dependencies(minimaldx minimaldx-project)
target_include_directories(minimaldx INTERFACE ${install_dir}/mod)
target_link_libraries(minimaldx INTERFACE ${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}MinimalDX${CMAKE_STATIC_LIBRARY_SUFFIX})
