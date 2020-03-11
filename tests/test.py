# TEB version 4.0.0 (https://github.com/teb-model/teb).
# Copyright 2018 D. Meyer. Licensed under CeCILL version 2.1.

import sys
import os
from pathlib import Path
import shutil
import argparse


THIS_DIR = Path(__file__).parent.resolve()
TEMP_DIR = THIS_DIR / 'temp'
PROJ_DIR = THIS_DIR.parent
sys.path.append(str(PROJ_DIR / 'tools'))

from helpers import compare


def main(build_type, case, allow_failure):

    # Make sure there are no previous versions installed
    shutil.rmtree(PROJ_DIR / 'build', ignore_errors=True)
    shutil.rmtree(PROJ_DIR / 'temp', ignore_errors=True)
    shutil.rmtree(PROJ_DIR / 'plots', ignore_errors=True)

    if case == 'integration':
        trial = '__THIS__integration'
        ref = 'master'
        case_name = 'CAPITOUL'
        patch_nml_trial = {'parameters': {'CCOOL_COIL': 'MinimalDX', 'LSOLAR_PANEL': True}}
        patch_nml_ref = {'parameters': {'CCOOL_COIL': 'MinimalDX', 'LSOLAR_PANEL': True}}
        patch_nml = [patch_nml_trial, patch_nml_ref]
        make = [False, False]

    elif case == 'make_cmake':
        # This test only works on Linux as TEB make only supported gcc.
        if (not os.name == 'posix') or (build_type != 'Debug'):
            print('test_make_cmake: This test only works on Linux -- skipping')
            return
        else:
            trial = '2aa67d0e961890704516567571a84e3e74776b59' # Last commit before psychrometric fixes
            ref = 'https://github.com/teb-model/teb/archive/3_sfx8.1.tar.gz'
            case_name = 'CAPITOUL'
            # In the old TEB, XTSTEP_SURF was hardcoded in the driver to
            # 300 s instead of 1800 s. Here we simply patch the namelist
            # before running TEB CMake as a workaround.
            patch_nml = {'parameters': {'XTSTEP_SURF': 300.}}
            patch_nml = [patch_nml, patch_nml]
            make = [True, False]

    elif case == 'psychrometrics':
        trial = 'a9c4c46889dca998ea9ca40e5e0edc62ecc75715'
        ref = '2aa67d0e961890704516567571a84e3e74776b59'
        case_name = 'CAPITOUL'
        patch_nml_trial = {'parameters': {'CCOOL_COIL': 'DXCOIL', 'LSOLAR_PANEL': True}}
        patch_nml_ref = {'parameters': {'CCOOL_COIL': 'DXCOIL', 'LSOLAR_PANEL': True}}
        patch_nml = [patch_nml_trial, patch_nml_ref]
        make = [False, False]

    elif case == 'minimal_dx':
        trial = '__THIS__minimal_dx'
        ref = 'master'
        case_name = 'CAPITOUL'
        patch_nml_trial = {'parameters': {'CCOOL_COIL': 'MinimalDX', 'LSOLAR_PANEL': True}}
        patch_nml_ref = {'parameters': {'CCOOL_COIL': 'DXCOIL', 'LSOLAR_PANEL': True}}
        patch_nml = [patch_nml_trial, patch_nml_ref]
        make = [False, False]

    else:
        raise RuntimeError('The correct case was not selected')

    print(f'Using configuration: {trial}, {ref}, {case_name}, {build_type}, {patch_nml}, {make}, {allow_failure}')

    compare(trial, ref, case_name, build_type, patch_nml, make, allow_failure)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--build_type', required=True, choices=['Debug', 'Release'],
                        help='CMAKE_BUILD_TYPE')
    parser.add_argument('--case', required=True, choices=['make_cmake', 'integration',
                                                          'psychrometrics', 'minimal_dx'],
                        help='The test case to run')
    parser.add_argument('--allow_failure', required=False, type=bool, nargs='?', const=True, default=False)
    args = parser.parse_args()

    print(f'CMAKE_BUILD_TYPE: {args.build_type}')
    print(f'case: {args.case}')
    print(f'allow-failure: {args.allow_failure}')
    main(args.build_type, args.case, args.allow_failure)
