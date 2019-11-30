import sys
from pathlib import Path
import shutil

import datetime
import pandas as pd
import f90nml

THIS_DIR = Path(__file__).parent.resolve()
TEMP_DIR = THIS_DIR / 'temp'
PROJ_DIR = THIS_DIR.parents[1]
print(PROJ_DIR)
sys.path.append(str(PROJ_DIR / 'tools'))

from helpers import run_teb, load_txt, get_date_params

def main():
    path_to_data = TEMP_DIR / 'CAPITOUL'
    path_to_make_exe = TEMP_DIR / 'teb-3_sfx8.1' / 'driver.exe'
    path_to_cmake_exe = TEMP_DIR / 'teb_cmake' / 'driver'

    path_to_namelist = path_to_data / 'input.nml'
    case_nml = f90nml.read(path_to_namelist)
    start, freq = get_date_params(case_nml)
    tz = 'UTC'

    # TEB Make
    run_teb(path_to_data, path_to_make_exe)
    shutil.move(path_to_data / 'output', path_to_data / 'output_make')
    df_make = load_txt(path_to_data / 'output_make', start, freq, tz=tz)

    # TEB CMake
    # In the old TEB, XTSTEP_SURF was hardcoded in the driver to
    # 300 s instead of 1800 s. Here we simply patch the namelist
    # before running TEB CMake as a workaround.
    path_to_unpatched_namelist = path_to_namelist.parent / 'input.nml_unpatched'
    shutil.move(path_to_namelist, path_to_unpatched_namelist )
    patch_nml = {'parameters': {'XTSTEP_SURF': 300.}}
    f90nml.patch(str(path_to_unpatched_namelist), patch_nml, str(path_to_namelist))

    run_teb(path_to_data, path_to_cmake_exe)
    df_cmake = load_txt(path_to_data / 'output', start, freq, tz=tz)
    # We have more outputs in the newer versions
    # drop any quantity not present in both versions.
    df_cmake = df_cmake[df_make.columns]

    diff = (df_make - df_cmake)
    num_unequal_samples = len(diff[diff.values > 0])
    if num_unequal_samples != 0:
        raise RuntimeError(f"{num_unequal_samples} samples are not equal")
    else:
        print("All output samples are equals")

if __name__ == "__main__":
    main()