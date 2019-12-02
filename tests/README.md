# Tests

This folder contains regression and integrations tests used to develop TEB. All tests are based on data from the [CAPITOUL experiment](../examples/CAPITOUL).

The current tests and results are as follows:

- [x] Make v CMake
- [ ] Psychrometric functions
- [ ] MinimalDX

## Make v CMake

``` bash
 python tests/test.py --build_type=Debug --case=make_cmake
```

Compares change in build system from Perl/Make to CMake/Make on Linux and Debug mode only.
The test passes with no discernable differences. Results are from commit ID: 

## Psychrometric functions


TODO:

## MinimalDX

TODO:
