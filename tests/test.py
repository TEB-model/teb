import sys
import os
from pathlib import Path
import shutil
import argparse
import datetime
import pandas as pd
import f90nml


THIS_DIR = Path(__file__).parent.resolve()
TEMP_DIR = THIS_DIR / 'temp'
PROJ_DIR = THIS_DIR.parent
sys.path.append(str(PROJ_DIR / 'tools'))

from helpers import compare


def test_make_cmake(build_type: str) -> None:
    # This test only works on Linux as TEB make only supported gcc.
    if (not os.name == 'posix') or (build_type != 'Debug'):
        print('test_make_cmake: This test only works on Linux -- skipping')
        return

    old_id = 'https://github.com/teb-model/teb/archive/3_sfx8.1.zip'
    new_id = 'test-make-cmake'
    case_name = 'CAPITOUL'
    # In the old TEB, XTSTEP_SURF was hardcoded in the driver to
    # 300 s instead of 1800 s. Here we simply patch the namelist
    # before running TEB CMake as a workaround.
    patch_nml = {'parameters': {'XTSTEP_SURF': 300.}}
    compare(old_id, new_id, case_name, build_type, [patch_nml, patch_nml], [True, False])
    return None

def test_integration(build_type: str) -> None:
    old_id = 'test-make-cmake'
    new_id = 'fix-psychrometrics'
    case_name = 'CAPITOUL'
    compare(old_id, new_id, case_name, build_type, [None, None], [False, False])
    return None

def test_minimal_dx(build_type: str) -> None:
    old_id = 'test-make-cmake'
    new_id = 'fix-psychrometrics'
    case_name = 'CAPITOUL'
    patch_nml = {'parameters': {'CCOOL_COIL': 'MinimalDX'}}
    compare(old_id, new_id, case_name, build_type, [None, patch_nml], [False, False])
    return None


if __name__ == "__main__":
    # Make sure there are no previous versions installed
    shutil.rmtree(PROJ_DIR / 'build', ignore_errors=True)
    shutil.rmtree(PROJ_DIR / 'temp', ignore_errors=True)

    parser = argparse.ArgumentParser()
    parser.add_argument('--build_type', required=True, choices=['Debug', 'Release'],
                        help='CMAKE_BUILD_TYPE')
    args = parser.parse_args()

    print(f'CMAKE_BUILD_TYPE: {args.build_type}')
    test_make_cmake(args.build_type)
    #test_integration(args.build_type)
    #test_minimal_dx(args.build_type)
