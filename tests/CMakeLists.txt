# MinimalDX version 0.1.1 (https://www.github.com/dmey/minimal-dx).
# Copyright 2018 D. Meyer and R. Raustad. Licensed under MIT.

include(ExternalProject)

set (EnergyPlus_GIT_REPO "https://github.com/dmey/EnergyPlusFortran.git")

if(NOT OFFLINE)
    SET(OFFLINE FALSE)
endif()
# Default back to pull repo from github unless specified with flag at configure time
# i.e. -DOFFLINE=ON
UNSET(OFFLINE CACHE)

# detect build config (like "Debug") for single- and multi-configuration generators
# (e.g. Makefile vs Visual Studio)
if(CMAKE_CONFIGURATION_TYPES)
    set(configuration $<CONFIGURATION>)
else()
    set(configuration ${CMAKE_BUILD_TYPE})
endif()

ExternalProject_Add(EnergyPlus-Fortran
    GIT_REPOSITORY ${EnergyPlus_GIT_REPO}
    CMAKE_GENERATOR ${CMAKE_GENERATOR}
    CMAKE_ARGS
        "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
        "-DCMAKE_INSTALL_PREFIX=install"
        "-DCMAKE_BUILD_TYPE=${configuration}"
        "-DUSE_SQLIGHT=OFF"
        "-DDISABLE_PSYCHROMETRICS_CACHE=TRUE"
    INSTALL_COMMAND
        "${CMAKE_COMMAND}" --build . --target install --config ${configuration}
    EXCLUDE_FROM_ALL TRUE
    STEP_TARGETS install
    UPDATE_DISCONNECTED ${OFFLINE}
)

# EnergyPlus library
ExternalProject_Get_Property(EnergyPlus-Fortran binary_dir)
set(ENERGYPLUS_FORTRAN_LIBRARY ${binary_dir}/install/lib/${CMAKE_STATIC_LIBRARY_PREFIX}energyplus${CMAKE_STATIC_LIBRARY_SUFFIX})
include_directories(${binary_dir}/install/mod)

add_subdirectory(cooling)
add_subdirectory(debugging)