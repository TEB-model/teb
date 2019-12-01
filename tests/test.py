import sys
import os
from pathlib import Path
import shutil

import datetime
import pandas as pd
import f90nml


THIS_DIR = Path(__file__).parent.resolve()
TEMP_DIR = THIS_DIR / 'temp'
PROJ_DIR = THIS_DIR.parent
print(PROJ_DIR)
sys.path.append(str(PROJ_DIR / 'tools'))

from helpers import run_teb, load_txt, get_date_params, build_teb, prepare_case, git_checkout, build_teb_make


def checkout_run_load(commit_id: str, case_name: str, patch_nml, download_zip=False) -> pd.DataFrame:
    if download_zip:
        path_to_exe, commit_id = build_teb_make(commit_id)
    else:
        git_checkout('test-make-cmake') # FIXME: this should be the ref id
        path_to_exe = build_teb(commit_id)

    path_to_case_dir = prepare_case(commit_id, case_name)
    run_teb(path_to_case_dir, path_to_exe, patch_nml)
    start, freq = get_date_params(path_to_case_dir / 'input.nml')
    df = load_txt(path_to_case_dir / 'output', start, freq, tz='UTC')
    return df

def compare(ref_id: str, trial_id: str, case_name: str, patch_nml) -> pd.DataFrame:
    df_ref = checkout_run_load(ref_id, case_name, patch_nml, True)
    df_trial = checkout_run_load(trial_id, case_name, patch_nml, False)
    # We have more outputs in the newer versions
    # drop any quantity not present in both versions.
    print(df_trial.columns)
    df_trial = df_trial[df_ref.columns]
    df_diff = df_ref - df_trial
    print(df_diff)
    num_unequal_samples = len(df_diff[df_diff.values > 0])
    if num_unequal_samples != 0:
        raise RuntimeError(f"{num_unequal_samples} samples are not equal")
    else:
        print("All output samples are equals")
    return df_diff



if __name__ == "__main__":
    shutil.rmtree(PROJ_DIR / 'build', ignore_errors=True)
    ref_id = 'https://github.com/teb-model/teb/archive/3_sfx8.1.zip'
    trial_id = 'test_2'
    case_name = 'CAPITOUL'
    patch_nml = {'parameters': {'XTSTEP_SURF': 300.}}
    diff = compare(ref_id, trial_id, case_name, patch_nml)
    print(diff)
