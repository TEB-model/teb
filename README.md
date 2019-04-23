# Town Energy Balance (TEB) model [![Build Status](https://travis-ci.com/dmey/teb-model.svg?token=3tqUbxqJuLtozjxqDymC&branch=master)](https://travis-ci.com/dmey/teb-model)

## Build instructions for Windows, Linux and macOS

To build the Town Energy Balance (TEB) executable and library, run the following commands:

```bash
$ mkdir build && cd build
$ cmake ..
$ cmake --build .
```

### Note for Windows Users

Make sure you have installed the Intel® Visual Studio Integration plugins or CMake will not be able to identify your compiler (`No CMAKE_Fortran_COMPILER could be found` error).
Make sure that you use Intel® Command-Line Window when launching CMake - The Intel® compiler provides a command-line window with the appropriate environment variables already set (see: [Using the Intel® Command-Line Window](https://software.intel.com/en-us/node/522358)).
You may also need to specify the generator flag `-G` in CMake; for example, if you are using Intel® Command-Line Window for Visual Studio 2017, then the CMake command should now be `cmake -G "Visual Studio 15 2017 Win64" ..`. For more information on how to specify generators in CMake see [cmake-generators](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html#visual-studio-generators)

E.g. on Windows using Intel® Command-Line Window for Visual Studio 2017:

```
> mkdir build && cd build
> cmake -G "Visual Studio 15 2017 Win64" ..
> cmake --build .
```

## Running the CAPITOUL test case

You can run a test case using the CAPITOUL data provided in `tests/CAPITOUL`. From your command prompt, execute the `driver` program from the `test/CAPITOUL` folder as you current working directory.

E.g. on Linux and macOS:

``` bash
cd test/CAPITOUL
../../build/driver
```

The output files are written to `tests/CAPITOUL/output`. You can compare them with the reference outputs in `tests/CAPITOUL/ref_output`. You can modify parameter-specific values or activate/deactivate part of TEB from the `input.nml` file.
