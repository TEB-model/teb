<div align="center">

<!-- omit in toc -->
# _Minimal_**DX**

[![Build Status](https://travis-ci.com/dmey/minimal-dx.svg?token=3tqUbxqJuLtozjxqDymC&branch=master)](https://travis-ci.com/dmey/minimal-dx)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3562310.svg)](https://doi.org/10.5281/zenodo.3562310)

[Overview](#overview) | [Prerequisites](#prerequisites) | [Build instructions](#build-instructions) | [Documentation](#documentation) | [Testing](#testing) | [How to cite](#how-to-cite) | [Copyright and License](#copyright-and-license)
</div>

## Overview
 _Minimal_**DX** is a standalone library for simulating the performance of single air-cooled cooling direct expansion (DX) coil -- i.e. a model for a very simple air-to-air air conditioner (AC).  _Minimal_**DX** is based on the [EnergyPlus](https://energyplus.net/) [Single-Speed Electric DX Air Cooling Coil](http://bigladdersoftware.com/epx/docs/8-7/engineering-reference/coils.html#single-speed-electric-dx-air-cooling-coil) model. It is meant to be a more stripped-down version of the EnergyPlus model and its main purpose is to be used as a standalone Fortran library in your own program for simulating basic AC equipment. All the functions used in this library were derived and adapted from [EnergyPlus (Fortran)](https://github.com/NREL/EnergyPlusRelease).

## Prerequisites

- [Git](https://git-scm.com/) >= 2
- [CMake](https://cmake.org/) >= 3.1
- A Recent version of GNU/Intel/Cray Fortran compiler
- [Python](https://www.python.org/) >= 3.6 [*Optional for testing*]

## Build instructions

These instructions will work on Windows (Visual Studio with Intel Fortran compilers), Linux (Intel or GNU compilers), and MacOS (Intel or GNU compilers).
To build the library without tests, clone this repository and run the following commands from your command prompt:

```
> mkdir build && cd build
> cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF ..
> cmake --build .
```

If you want to specify the installation directory and install, use the `CMAKE_INSTALL_PREFIX` flag and run `cmake --build . --target install`. The `BUILD_TESTING` flag is used to allow users to carry out regression tests -- see [Testing](#Testing) below. By default, we set the real type to an 8 byte wide. This behaviour is controlled by the optional `USE_REAL8` flag (default ON) -- the use of real 8 is required when building with `BUILD_TESTING`.


<details>
<summary>Note for Windows Users</summary>

Make sure you have installed the Intel® Visual Studio Integration plugins or CMake will not be able to identify your compiler (No CMAKE_Fortran_COMPILER could be found error). Make sure that you use Intel® Command-Line Window when launching CMake - The Intel® compiler provides a command-line window with the appropriate environment variables already set (see: [Using the Intel® Command-Line Window](https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-using-the-command-line-on-windows)). You may also need to specify the generator flag -G in CMake; for example, if you are using Intel® Command-Line Window for Visual Studio 2010, then the CMake command should now be `cmake -G "Visual Studio 10 2010" ..`. For more information on how to specify generators in CMake see [cmake-generators](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html#visual-studio-generators).
</details>

## Documentation

The main subroutine used to calculate the performance of a single coil air-cooled DX cooling equipment is `MinimalDXCoolingDriver`.

- `MinimalDXCoolingDriver` was adapted from EnergyPlus `CalcDoe2DXCoil`. The wrapper to call is in `src/cooling/minimal_dx_cooling_driver.f90`, the main calculations are carried out in `src/cooling/minimal_dx_cooling.f90`.


The API documentation can be found at [https://dmey.github.io/minimal-dx](https://dmey.github.io/minimal-dx).

The original EnergyPlus documentation can be found in the [EnergyPlus Engineering Reference](https://bigladdersoftware.com/epx/docs/8-7/engineering-reference/):
- [Single-Speed Electric DX Air Cooling Coil](http://bigladdersoftware.com/epx/docs/8-7/engineering-reference/coils.html#single-speed-electric-dx-air-cooling-coil)


## Testing

Regression tests are useful to compare the two implementations over a large number of different conditions. In order to run these tests you will need to have Python 3.5 or above installed on your system and had enabled testing with the `BUILD_TESTING` flag.

All files used to generate the tests can be found under the `tests` folder. `make_test_data_cooling.f90` and `test_cooling.f90` are responsible for generating the data and invoking the functions respectively. The program `run_tests.py` is responsible for calling the `make_test_data_cooling` and `test_cooling` programs and for generating plots.

To run the regression tests, first ensure that you have built _Minimal_**DX** with `BUILD_TESTING=ON` and in debug mode (`CMAKE_BUILD_TYPE=Debug`), then go to the `tests` folder and run the following commands:

```
pip3 install -r requirements.txt
python3 run_tests.py
```

All output data and plots will be written to `build/outputs`.

### Offline build

When enabling build tests `BUILD_TESTING=ON`,  _Minimal_**DX** will require a working internet connection to download up-to-date EnergyPlus dependencies. If these dependencies have already been downloaded and you do not have a working internet connection, you can still build _Minimal_**DX** with tests enabled by passing the the offline flag (i.e. `-DOFFLINE=ON`) at configure time.

## How to cite

If you are using _Minimal_**DX**, please make sure to cite the specific version you are using. For a list of avalable versions and Digital Object Identifiers (DOI), please refer to the [_Minimal_**DX** archive on Zenodo](https://doi.org/10.5281/zenodo.3562310).

## Copyright and License

Copyright 2018-2020 D. Meyer and R. Raustad. Licensed under [MIT](LICENSE.txt).