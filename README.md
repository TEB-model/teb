# Town Energy Balance (TEB) model [![Build Status](https://travis-ci.com/teb-model/teb.svg?branch=master)](https://travis-ci.com/teb-model/teb)

This simple driver for the TEB (Town Energy Balance; [Masson, 2000](https://dx.doi.org/10.1023/A:1002463829265) and subsequent papers), is intended to be used by scientists whishing to implement TEB in their own software. This simple driver provides all the source code of TEB alone, making it easier to analyse and to integrate in its own atmospheric model for example.

If you wish to use TEB for physical simulations, without intention to extract all the TEB routines into another software environment, you could use the SURFEX platform. It contains TEB and much more, especially the ISBA scheme for vegetation, or several I/O formats, including NetCDF (see http://www.cnrm.meteo.fr/surfex/).


## Documentation

Please refer to the [SURFEX website](http://www.umr-cnrm.fr/surfex/IMG/pdf/surfex_scidoc_v8.1.pdf) for the complete documentation on TEB. The API documentation is available [here](https://teb-model.github.io/teb).

## Installing


### Prerequisites

- [Git](https://git-scm.com/) version 2.13 or above
- [CMake](https://cmake.org/) version 3.1.0 or above
- GNU/Intel/Cray Fortran compiler
- [Python](https://www.python.org/) 3.6+ [*Optional*: for testing and tutorial]


### Build instructions for Windows, Linux and macOS

To build the Town Energy Balance (TEB) executable and library, run the following commands from the command prompt/terminal window:

```bash
git clone --recurse-submodules https://github.com/teb-model/teb.git
cd teb
```

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```


#### Note for Windows Users

Make sure you have installed the Intel® Visual Studio Integration plugins or CMake will not be able to identify your compiler (`No CMAKE_Fortran_COMPILER could be found` error).
Make sure that you use Intel® Command-Line Window when launching CMake - The Intel® compiler provides a command-line window with the appropriate environment variables already set (see: [Using the Intel® Command-Line Window](https://software.intel.com/en-us/node/522358)).
You may also need to specify the generator flag `-G` in CMake; for example, if you are using Intel® Command-Line Window for Visual Studio 2017, then the CMake command should now be `cmake -G "Visual Studio 15 2017 Win64" ..`. For more information on how to specify generators in CMake see [cmake-generators](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html#visual-studio-generators)

E.g. on Windows using Intel® Command-Line Window for Visual Studio 2017:

```powershell
mkdir build && cd build
cmake -G "Visual Studio 15 2017 Win64" ..
cmake --build .
```

By default, we set the real type to an 8 byte wide. This behaviour is controlled by the optional `USE_REAL8` flag (default ON).

## Testing

Tests are run using the CAPITOUL data provided in `examples/CAPITOUL`. Tests are found in `tests/tests.py` and an overview is given in `tests/README.md`[tests/README.md]. All tests are automatically run at every commit using Continous Integration. If you are looking to run your tests locally on Linux or macOS, first make sure you have installed all the [prerequisites](#prerequisites), then from the command prompt:

``` bash
python -m pip install -r requirements.txt
python tests/test.py --build_type=<BUILD_TYPE>
```

where `<BUILD_TYPE>` is either `Debug` or `Release`. The output files are written to `temp`. Plots are written to `plots`.


## How to cite

When using the TEB software, please cite both model, and software (with version) as follow:

| Model        | Software          | Version*                                              |
| ------------ | ----------------- | ---------------------------------------------------- |
| Masson, 2000 | XX et al, XXXX    | [see Zenodo](https://XXXX) |


The corresponding reference list should be as follows:

> Masson, V., 2000: A Physically-Based Scheme For The Urban Energy Budget In Atmospheric Models. Boundary-Layer Meteorology, 94, 357–397, https://doi.org/10.1023/A:1002463829265.

> XX et al, XXX


*please make sure to cite the same version you are using with the correct DOI. For a list of all avalable versions see the list of versions on [Zenodo](XXXX).


## Copyright and license

Copyright name and year is clearly stated at the top of each source file. All software released under [CeCILL-C version 1](Licence_CeCILL-C_V1-en.txt).
