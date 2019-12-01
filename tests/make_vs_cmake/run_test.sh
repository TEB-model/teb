#!/usr/bin/env bash

# TEB version x.x.x (https://github.com/teb-model/teb).
# Copyright 2018 D. Meyer. Licensed under CeCILL-C version 1.

set -xe

if [ ! -f compare_make_cmake.py ]; then
    echo "Please run this script from the make_cmake test folder"
    exit 1
fi

mkdir -p temp && pushd temp
# Remove artefacts from previous tests.
rm -rf *


# Build TEB make variant
wget -O teb-3_sfx8.1.tar.gz https://codeload.github.com/teb-model/teb/tar.gz/3_sfx8.1
tar xvzf teb-3_sfx8.1.tar.gz

pushd teb-3_sfx8.1
./mkmf.pl -t gfortran_args -p driver.exe src_driver src_struct src_proxi_SVAT src_solar src_teb
make
popd

# Build TEB CMake variant
# TODO: git checkout commit
mkdir teb_cmake && pushd teb_cmake
# Can only compare against Debug as no other congigs were used in
# original Makefile.
cmake -DCMAKE_BUILD_TYPE=Debug ../../../.. -LA
VERBOSE=1 make
popd
popd

# CAPITOUL case
cp -r ../CAPITOUL temp/CAPITOUL

python compare_make_cmake.py
