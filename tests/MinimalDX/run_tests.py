# TEB version x.x.x (https://github.com/teb-model/teb).
# Copyright 2018 D. Meyer. Licensed under CeCILL-C version 1.

import os
import shutil
import subprocess
from pathlib import Path

from datetime import datetime, timedelta

from typing import List

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import f90nml


THIS_DIR = Path(__file__).resolve().parent
CAPITOUL_DIR = THIS_DIR.parent / 'CAPITOUL'
PROJ_DIR = THIS_DIR.parents[1]
# TEB does not output cols of datetime, we need to know when the sim started.
# The start date in `input.nml` for the 'CAPITOUL' data is '2004-02-20 00:00:00'
# however the first data point in the output files is for 30 mins later -- i.e. '2004-02-20 00:30:00'.
TEB_OUT_START_DATETIME = '2004-02-20 00:30:00'


def run_teb_coil_case(case_name: str, path_to_cases: Path) -> None:
    """ Run TEB with different CCOOL_COIL configuration.
    """

    # Check environment to determine folder structure
    if os.name == 'nt':
        path_to_exe = PROJ_DIR / 'build' / 'Release' / 'driver.exe'
    elif os.name == 'posix' or 'mac':
        path_to_exe = PROJ_DIR / 'build' / 'driver'
    else:
        raise Exception('Error: this program cannot run on your system')

    case_dir = path_to_cases / case_name
    Path(case_dir).mkdir(parents=True, exist_ok=True)
    # TEB needs the output folder to be present.
    Path(case_dir / 'output').mkdir(parents=True, exist_ok=True)

    # Use the CAPITOUL data as input for the experiment.
    shutil.copytree(CAPITOUL_DIR / 'input', case_dir / 'input')

    # We only want to change the coil to use in the system.
    patch_nml = {'parameters': {'CCOOL_COIL': case_name}}
    f90nml.patch(CAPITOUL_DIR / 'input.nml', patch_nml, case_dir / 'input.nml')

    with open(case_dir / str('driver_' + case_name + '.log'), 'w') as f:
        print('Running case: ' + case_name)
        subprocess.check_call([path_to_exe], cwd=case_dir, stdout=f)

def get_datetime_end(df, base, freq):
    base_datetime = datetime.strptime(base, '%Y-%m-%d %H:%M:%S')
    df_length = len(df.index)
    end = base_datetime + freq * (df_length - 1)
    return end


def load_outputs_to_df(case_name: str, path_to_cases: Path) -> pd.DataFrame:
    path_to_outdir = path_to_cases / case_name / 'output'
    df_list = []
    for filename in path_to_outdir.glob('*.txt'):
        df = pd.read_csv(path_to_outdir / filename, delim_whitespace=True,
                                       skip_blank_lines=True, header=None, names=[str(case_name + ':' + filename.stem)])
        end = get_datetime_end(df, TEB_OUT_START_DATETIME, timedelta(minutes=30))
        date_rng = pd.date_range(start=TEB_OUT_START_DATETIME, end=end, freq='30T')
        df = df.set_index(date_rng)
        df = df.tz_localize(tz='UTC')
        df_list.append(df)
    df = pd.concat(df_list, axis=1)
    return df


def plot_results(df_list: List[pd.DataFrame], path_to_outdir: Path) -> None:
    assert len(df_list) == 2
    assert df_list[0].size == df_list[1].size
    for index, col_name in enumerate(df_list[0]):
        print('Plotting: ' + col_name)
        styles = ['r-', 'k--']
        for idx, df in enumerate(df_list):
            df.iloc[:, index].plot(style=styles[idx])
        plt.legend()
        plt.ylabel('See var description')
        plt.savefig(path_to_outdir / str(col_name.split(':')[1] + '.png'))
        plt.close()


if __name__ == "__main__":
    # Clean up from previous runs.
    shutil.rmtree(THIS_DIR / 'artefacts', ignore_errors=True)
    shutil.rmtree(THIS_DIR / 'plots', ignore_errors=True)
    # Create plot dir.
    Path(THIS_DIR / 'plots').mkdir()
    # Run the two cases ('MinimalDX' and 'DXCOIL') separately and compare the results at the end.
    CASES = ['MinimalDX', 'DXCOIL']
    df_list = []
    for case in CASES:
        run_teb_coil_case(case, THIS_DIR / 'artefacts')
        df = load_outputs_to_df(case, THIS_DIR / 'artefacts')
        df_list.append(df)
    plot_results(df_list, THIS_DIR / 'plots')