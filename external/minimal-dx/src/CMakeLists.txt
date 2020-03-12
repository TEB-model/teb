# MinimalDX version 0.1.1 (https://www.github.com/dmey/minimal-dx).
# Copyright 2018 D. Meyer and R. Raustad. Licensed under MIT.

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library(MinimalDX STATIC
    ../external/EnergyPlus/src/EnergyPlus/eplus_psychro.f90
    ../external/psychrolib/src/fortran/psychrolib.f90
    cooling/minimal_dx_cooling.f90
    cooling/minimal_dx_cooling_driver.f90
    shared/minimal_dx_fan.f90
    shared/psychro_wrapper.f90
)

target_include_directories(MinimalDX
    INTERFACE ${CMAKE_Fortran_MODULE_DIRECTORY}
)

install(TARGETS MinimalDX
    ARCHIVE DESTINATION lib
)

# install Fortran mod files, see https://cmake.org/pipermail/cmake/2010-May/037106.html
install(CODE
  "EXECUTE_PROCESS (COMMAND \"${CMAKE_COMMAND}\" -E copy_directory \"${CMAKE_Fortran_MODULE_DIRECTORY}/\${BUILD_TYPE}\" \"${CMAKE_INSTALL_PREFIX}/mod\")"
)