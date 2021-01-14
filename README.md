<div align="center">

<!-- omit in toc -->
# The Town Energy Balance (TEB) model

[![GitHub release (latest by date)](https://img.shields.io/github/v/release/TEB-model/teb)](https://github.com/TEB-model/teb/releases/latest) [![CI](https://github.com/TEB-model/teb/workflows/CI/badge.svg)](https://github.com/TEB-model/teb/actions) [![DOI](https://joss.theoj.org/papers/10.21105/joss.02008/status.svg)](https://doi.org/10.21105/joss.02008) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3887080.svg)](https://doi.org/10.5281/zenodo.3887080)

[Overview](#overview) | [Prerequisites](#prerequisites) | [Installation](#installation) | [Documentation](#documentation) | [Example application](#example-application) | [Testing](#testing) | [How to cite](#how-to-cite) | [Contributing](#contributing) | [Copyright and license](#copyright-and-license)
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

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```

By default, we set the real type to an 8 byte wide. This behavior is controlled by the optional `USE_REAL8` flag (default ON).

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

This section includes links to software and model documentation. If you are new to TEB and simply looking to easily get started, please refer to the [example application](#example-application) instead.

### Software

See [software documentation](docs/software-docs.md) for a general overview. For configuration options available at runtime please refer to [namelist options](docs/namelist-options.md).

### Scientific
The complete scientific documentation is included in the [SURFEX scientific documentation](http://www.umr-cnrm.fr/surfex/IMG/pdf/surfex_scidoc_v8.1.pdf).


### Code browser
The code browser is available at [https://teb-model.github.io/teb](https://teb-model.github.io/teb).


## Example application

To get started with TEB, see [`examples/CAPITOUL/README.md`](examples/CAPITOUL/README.md). The folder contains a simple tutorial in [Jupyter Notebook](https://jupyter.org/) to estimate the buildings' power demand for cooling using TEB and data and parameters from the CAPITOUL (Canopy and Aerosol Particles Interactions in TOulouse Urban Layer; [Masson et al., 2008](https://doi.org/10.1007/s00703-008-0289-4)) campaign.


## Testing

Tests are run using the CAPITOUL data provided in `examples/CAPITOUL`. Tests are found in `tests/tests.py` and an overview is given in [`tests/README.md`](tests/README.md). All tests are automatically run at every commit using Continuos Integration. If you are looking to run your tests locally on Linux or macOS, first make sure you have installed all the [prerequisites](#prerequisites), then from the command prompt:

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

| Physical model                                          | This software                                             | Version*                                             |
| ------------------------------------------------------- | --------------------------------------------------------- | ---------------------------------------------------- |
| [Masson, 2000](https://doi.org/10.1023/A:1002463829265) | [Meyer et al., 2020](https://doi.org/10.21105/joss.02008) | [see Zenodo](https://doi.org/10.5281/zenodo.3887080) |


The corresponding reference list should be as follows:

> Masson, V., 2000: A Physically-Based Scheme For The Urban Energy Budget In Atmospheric Models. Boundary-Layer Meteorology, 94, 357–397, https://doi.org/10.1023/A:1002463829265.

> Meyer, D., Schoetter, R., Masson, V., Grimmond, S., 2020: Enhanced software and platform for the Town Energy Balance (TEB) model. Journal of Open Source Software, 5(50), 2008. https://doi.org/10.21105/joss.02008.

*please make sure to cite the same version you are using with the correct DOI. For a list of all available versions see the list of versions on [Zenodo](https://doi.org/10.5281/zenodo.3887080).


## Contributing

Please see the [CONTRIBUTING.md](CONTRIBUTING.md) file.


## Copyright and license

Copyright stated at the top of source files. Software released under [CeCILL version 2.1](Licence_CeCILL_V2.1-en.txt).
Third party files released under MIT or CeCILL-C. See [ThirdPartyNotices.txt](ThirdPartyNotices.txt) for more details.
