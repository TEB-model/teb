cmake_minimum_required(VERSION 3.1)
project(TEB Fortran)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Debug")
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type (default Debug)" FORCE)
endif()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # -fimplicit-none
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -fdefault-real-8")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -ffpe-trap=invalid,zero,overflow,underflow -fcheck=all")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    # -implicitnone
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -traceback -fpe-all0 -check all")
else()
    message(FATAL_ERROR "Only gfortran and ifort compilers supported")
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)

set(src_driver
    src/driver/abor1_sfx.F90
    src/driver/add_forecast_to_date_surf.F90
    src/driver/close_file.F90
    src/driver/close_file_asc.F90
    src/driver/driver.F90
    src/driver/modd_arch.F90
    src/driver/modd_forc_atm.F90
    src/driver/modd_reprod_oper.F90
    src/driver/modd_surf_conf.F90
    src/driver/mode_char2real.F90
    src/driver/ol_alloc_atm.F90
    src/driver/ol_read_atm.F90
    src/driver/ol_read_atm_ascii.F90
    src/driver/ol_time_interp_atm.F90
    src/driver/open_close_bin_asc_forc.F90
    src/driver/read_surf_atm.F90
)

file(GLOB_RECURSE src_proxi_SVAT "src/proxi_SVAT/*")
file(GLOB_RECURSE src_solar "src/solar/*")
file(GLOB_RECURSE src_struct "src/struct/*")
file(GLOB_RECURSE src_teb "src/teb/*")

add_library(teb STATIC
    ${src_driver}
    ${src_proxi_SVAT}
    ${src_solar}
    ${src_struct}
    ${src_teb}
)

add_executable(driver src/driver/driver.F90)
target_link_libraries(driver teb)

target_include_directories(teb
 INTERFACE ${CMAKE_Fortran_MODULE_DIRECTORY}
)