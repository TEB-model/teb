<div align="center">

<!-- omit in toc -->
# The Town Energy Balance (TEB) model

[![Build Status](https://travis-ci.com/teb-model/teb.svg?branch=master)](https://travis-ci.com/teb-model/teb)

[Overview](#overview) | [Prerequisites](#prerequisites) | [Installation](#installation) | [Documentation](#documentation) | [Example application](#example-application) | [Testing](#testing) | [How to cite](#how-to-cite) | [Contributing](#contributing) | [Copyright and license](#copyright-and-license) | [References](#references)
</div>

## Overview

This enhanced software and platform for TEB (Town Energy Balance; [Masson, 2000](https://dx.doi.org/10.1023/A:1002463829265) and subsequent papers), is intended to help scientists and practitioners wishing to use the TEB model in their research as a standalone software application or as a library to calculate the urban surface energy balance at neighborhood scale assuming a simplified canyon geometry.


## Prerequisites

- [Git](https://git-scm.com/) >= 2
- [CMake](https://cmake.org/) >= 3.1
- A Recent version of GNU/Intel/Cray Fortran compiler
- [Python](https://www.python.org/) >= 3.7 [*Optional for testing and tutorial*]


## Installation

To build the Town Energy Balance (TEB) executable and library on Windows, Linux and macOS, clone this repository and run the following commands from your command prompt:

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```

By default, we set the real type to an 8 byte wide. This behaviour is controlled by the optional `USE_REAL8` flag (default ON).

<details>
<summary>Note for Windows Users</summary>

Make sure you have installed the Intel® Visual Studio Integration plugins or CMake will not be able to identify your compiler (`No CMAKE_Fortran_COMPILER could be found` error).
Make sure that you use Intel® Command-Line Window when launching CMake - The Intel® compiler provides a command-line window with the appropriate environment variables already set (see: [Using the Intel® Command-Line Window](https://software.intel.com/en-us/node/522358)).
You may also need to specify the generator flag `-G` in CMake; for example, if you are using Intel® Command-Line Window for Visual Studio 2017, then the CMake command should now be `cmake -G "Visual Studio 15 2017 Win64" ..`. For more information on how to specify generators in CMake see [cmake-generators](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html#visual-studio-generators)

E.g. on Windows using Intel® Command-Line Window for Visual Studio 2017:

```powershell
mkdir build && cd build
cmake -G "Visual Studio 15 2017 Win64" ..
cmake --build .
```
</details>


## Documentation

This section includes links to software and model documentation. If you are new to TEB and simply looking to easily get started, please refer to the [Example application](#example-application) instead.

### Software

See the [software documentation](docs/software-docs.md) for a general overview. For configuration options available at runtime please refer to [namelist options](docs/namelist-options.md).

### Scientific
The complete scientific documentation is included with the [SURFEX scientific documentation](http://www.umr-cnrm.fr/surfex/IMG/pdf/surfex_scidoc_v8.1.pdf).


### Code browser
The API documentation and code browser is available at [https://teb-model.github.io/teb](https://teb-model.github.io/teb).


## Example application

[`tutorial.ipynb`](examples/CAPITOUL/tutorial.ipynb) provides a simple tutorial detailing how TEB can be used to to estimate the buildings' energy demand for cooling using sample data and parameters from the CAPITOUL (Canopy and Aerosol Particles Interactions in TOulouse Urban Layer; Masson et al., 2008) campaign.


## Testing

Tests are run using the CAPITOUL data provided in `examples/CAPITOUL`. Tests are found in `tests/tests.py` and an overview is given in [`tests/README.md`](tests/README.md). All tests are automatically run at every commit using Continous Integration. If you are looking to run your tests locally on Linux or macOS, first make sure you have installed all the [prerequisites](#prerequisites), then from the command prompt:

``` bash
python -m pip install -r requirements.txt
python tests/test.py --build_type=<BUILD_TYPE> --case=<TEST_CASE>
```

where `<BUILD_TYPE>` is either `Debug` or `Release` and `<TEST_CASE>` is the test case currently supported in `tests/test.py` -- see [`tests/README.md`](tests/README.md) for more information. The output files are written to `temp`. Plots are written to `plots`.


E.g.

``` bash
python -m pip install -r requirements.txt
python tests/test.py --build_type=Debug --case=integration
```

## How to cite

When using the TEB software, please cite both model, and software (with version) as follows:

| Model        | Software | Version*                   |
| ------------ | -------- | -------------------------- |
| Masson, 2000 | XXXX     | [see Zenodo](https://XXXX) |


The corresponding reference list should be as follows:

> Masson, V., 2000: A Physically-Based Scheme For The Urban Energy Budget In Atmospheric Models. Boundary-Layer Meteorology, 94, 357–397, https://doi.org/10.1023/A:1002463829265.

> XXXX *In preparation*


*please make sure to cite the same version you are using with the correct DOI. For a list of all available versions see the list of versions on [Zenodo](XXXX).


## Contributing

Please see the [CONTRIBUTING.md](CONTRIBUTING.md) file.


## Copyright and license

General copyright name and year is clearly stated at the top of each source file. All software released under [CeCILL version 2.1](Licence_CeCILL_V2.1-en.txt).


## References

> Masson, V., 2000: A Physically-Based Scheme For The Urban Energy Budget In Atmospheric Models. Boundary-Layer Meteorology, 94, 357–397, https://doi.org/10.1023/A:1002463829265.

> Masson, V., and Coauthors, 2008: The Canopy and Aerosol Particles Interactions in TOulouse Urban Layer (CAPITOUL) experiment. Meteorol Atmos Phys, 102, 135–157, https://doi.org/10.1007/s00703-008-0289-4.
