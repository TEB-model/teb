# TEB version x.x.x (https://github.com/teb-model/teb).
# Copyright 2018 D. Meyer. Licensed under CeCILL-C version 1.

import os
import subprocess
from pathlib import Path

import datetime

from typing import List

import pandas as pd
import f90nml


def run_teb(path_to_case_dir: Path, path_to_exe: Path) -> None:
    """Helpers to run teb from python.

    Parameters
    ----------
    path_to_case_dir : Path
        [description]
    path_to_exe : Path
        [description]

    Raises
    ------
    RuntimeError
        [description]
    RuntimeError
        [description]
    RuntimeError
        [description]
    """
    # Sanity check
    if not list(path_to_case_dir.glob('*.nml')):
        raise RuntimeError(f'Namelist not found in {path_to_case_dir}')
    if not (path_to_case_dir / 'input').exists():
        raise RuntimeError(f'Input folder not found in {path_to_case_dir}')
    if not (path_to_exe).exists():
        raise RuntimeError(f'{path_to_exe} not found.')

    # TEB needs the output folder to be present.
    Path(path_to_case_dir / 'output').mkdir(parents=True, exist_ok=True)

    case_name = path_to_case_dir.stem
    with open(path_to_case_dir / str(case_name + '.log'), 'w') as f:
        print('Running case: ' + case_name)
        subprocess.check_call([path_to_exe], cwd=path_to_case_dir, stdout=f)


def get_date_params(path_to_nml: Path):
    """Returns start and freq from a teb namelist file.
    
    Parameters
    ----------
    path_to_nml : Path
        [description]
    
    Returns
    -------
    [type]
        [description]
    """

    params = path_to_nml['parameters']
    freq = datetime.timedelta(seconds=params['xtstep_surf'])
    sec_from_start = datetime.timedelta(seconds=params['ztime_start'])
    start = datetime.datetime(year=params['iyear'], month=params['imonth'], day=params['iday']) + sec_from_start
    return start, freq


def load_txt(path_to_files: Path, start: datetime.datetime,
                            freq: datetime.timedelta, tz: str) -> pd.DataFrame:
    """Returns a Pandas DataFrame given the path contaning teb inputs/outputs.

    Parameters
    ----------
    path_to_files : Path
        Path to teb input/output folder
    start : datetime.datetime
        Start date
    freq : datetime.timedelta
        Frequency of data/time step
    tz : str
        Timezone

    Returns
    -------
    pd.DataFrame
        Timezone localized Pandas DataFrame contaning teb inputs/outputs.
    """
    ts_list = []
    for filename in path_to_files.glob('*.txt'):
        ts = pd.read_csv(filename, delim_whitespace=True, skip_blank_lines=True,
                         header=None, names=[filename.stem])
        end = start + freq * (len(ts.index) - 1)
        rng = pd.date_range(start=start, end=end, freq=freq)
        ts = ts.set_index(rng)
        ts = ts.tz_localize(tz=tz)
        ts_list.append(ts)
    df = pd.concat(ts_list, axis=1)
    return df

if __name__ == "__main__":
    pass
