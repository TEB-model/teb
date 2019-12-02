# Tests

This folder contains regression and integrations tests used to develop TEB. All tests are based on data from the [CAPITOUL experiment](../examples/CAPITOUL/).

An `integration` test is run at every commit through continuous integration (CI) to check that code changes do not impact the results. However, when bug fixes / code improvements are implemented, results will change and equality from the previous version can no longer be expected. In such case, the PR or commit message will report details about these changes. In the even of major changes, a small test should also be implemented in `tests/test.py` so that such differences can be easily reproduced locally.

The self-contained tests that have been implemented to this date are:

## Make v CMake

Compares change in build system from Perl/Make to CMake/Make on Linux and Debug mode only.
The test passes with no discernable differences.

``` bash
# Only avalable in Debug mode and for Linux only.
python tests/test.py --build_type=Debug --case=make_cmake
```

## Psychrometric functions

TODO:

## MinimalDX

Compares changes due to a new DX cooling coil as implemented in the MinimalDX library.
Changes are large as expected given the different typo of coils and improvements in the calculations.

``` bash
python tests/test.py --build_type=Debug --case=integration
python tests/test.py --build_type=Release --case=integration
```
