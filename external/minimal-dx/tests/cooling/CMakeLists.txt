# MinimalDX version 0.1.1 (https://www.github.com/dmey/minimal-dx).
# Copyright 2018 D. Meyer and R. Raustad. Licensed under MIT.

# Common files
set(SRC_TEST_COOLING
    ../../src/cooling/minimal_dx_cooling_driver.f90
    ../../src/cooling/minimal_dx_cooling.f90
    ../../src/shared/minimal_dx_fan.f90
    eplus_wrapper_cooling.f90
    test_cooling.f90
)

# Program to test the new implementation using the PsychroLib library
# against the original EnergyPlus implementation
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/test_cooling)
add_executable(test_cooling
    ${SRC_TEST_COOLING}
    ../../src/shared/psychro_wrapper.f90
    ../../external/EnergyPlus/src/EnergyPlus/eplus_psychro.f90
    ../../external/psychrolib/src/fortran/psychrolib.f90
)
target_link_libraries(test_cooling ${ENERGYPLUS_FORTRAN_LIBRARY})
add_dependencies(test_cooling EnergyPlus-Fortran-install)
install(TARGETS test_cooling
    RUNTIME DESTINATION bin
)

# Program to test the new implementation using the EnergyPlus psychrometric functions
# against the original EnergyPlus implementation
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/test_cooling_eplus)
add_executable(test_cooling_eplus
    ${SRC_TEST_COOLING}
    ../shared/psychro_wrapper_eplus.f90
)
target_link_libraries(test_cooling_eplus ${ENERGYPLUS_FORTRAN_LIBRARY})
add_dependencies(test_cooling_eplus EnergyPlus-Fortran-install)
install(TARGETS test_cooling_eplus
    RUNTIME DESTINATION bin
)

# Program to create test data
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/make_test_data_cooling)
add_executable(make_test_data_cooling make_test_data_cooling.f90)
add_dependencies(make_test_data_cooling EnergyPlus-Fortran-install)
target_link_libraries(make_test_data_cooling ${ENERGYPLUS_FORTRAN_LIBRARY})
install(TARGETS make_test_data_cooling
    RUNTIME DESTINATION bin
)