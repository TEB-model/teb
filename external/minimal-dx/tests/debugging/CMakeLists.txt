# MinimalDX version 0.1.1 (https://www.github.com/dmey/minimal-dx).
# Copyright 2018 D. Meyer and R. Raustad. Licensed under MIT.

# Program used for debugging
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/test_cooling_point)
add_executable(test_cooling_point
    test_cooling_point.f90
    ../cooling/eplus_wrapper_cooling.f90
    ../shared/psychro_wrapper_eplus.f90
    ../../src/cooling/minimal_dx_cooling_driver.f90
    ../../src/cooling/minimal_dx_cooling.f90
    ../../src/shared/minimal_dx_fan.f90
    ../../external/EnergyPlus/src/EnergyPlus/eplus_psychro.f90
    ../../external/psychrolib/src/fortran/psychrolib.f90
)
target_link_libraries(test_cooling_point ${ENERGYPLUS_FORTRAN_LIBRARY})
add_dependencies(test_cooling_point EnergyPlus-Fortran-install)
install(TARGETS test_cooling_point
    RUNTIME DESTINATION bin
)