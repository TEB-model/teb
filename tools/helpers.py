# TEB version x.x.x (https://github.com/teb-model/teb).
# Copyright 2018 D. Meyer. Licensed under CeCILL-C version 1.

import os
import subprocess
from pathlib import Path
import shutil

import datetime

from typing import List

import pandas as pd
import f90nml

from urllib.parse import urlparse
import urllib.request

import zipfile

THIS_DIR = Path(__file__).parent.resolve()
PROJ_DIR = THIS_DIR.parent


def run_teb(path_to_case_dir: Path, path_to_exe: Path, patch_nml) -> None:
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

    if patch_nml:
        path_to_namelist = path_to_case_dir / 'input.nml'
        path_to_unpatched_namelist = path_to_case_dir / 'input.nml_unpatched'
        shutil.move(path_to_namelist, path_to_unpatched_namelist)
        f90nml.patch(str(path_to_unpatched_namelist), patch_nml, str(path_to_namelist))

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
    case_nml = f90nml.read(path_to_nml)
    params = case_nml['parameters']
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

def build_teb(commit_id: str, build_type: str):
    path_to_exe_dir = PROJ_DIR / 'build' / commit_id
    path_to_exe_dir.mkdir(parents=True, exist_ok=True)
    config_command = ['cmake', f'-DCMAKE_BUILD_TYPE={build_type}', str(PROJ_DIR), '-LA']
    build_command = ['make']
    with open(path_to_exe_dir / str(commit_id + '_config.log'), 'w') as f:
        print('Configuring TEB for case: ' + commit_id)
        subprocess.check_call(config_command, cwd=path_to_exe_dir, stdout=f)
    with open(path_to_exe_dir / str(commit_id + '_build.log'), 'w') as f:
        print('Building TEB for case: ' + commit_id)
        subprocess.check_call(build_command, cwd=path_to_exe_dir, stdout=f)
    path_to_exe = path_to_exe_dir / 'driver'
    return path_to_exe

def prepare_case(commit_id, case_name):
    path_to_case_dir = PROJ_DIR / 'temp' / commit_id / case_name
    if path_to_case_dir.exists():
        shutil.rmtree(path_to_case_dir)
    shutil.copytree(PROJ_DIR / 'examples' / case_name, path_to_case_dir)
    return path_to_case_dir

def git_checkout(commit_id: str):
    build_folder = PROJ_DIR / 'build' / commit_id
    build_folder.mkdir(parents=True, exist_ok=True)
    command = ['git', 'checkout', str(commit_id)]
    with open(build_folder / str(commit_id + '_git.log'), 'w') as f:
        print('checking out case' + commit_id)
        subprocess.check_call(command, cwd=build_folder, stdout=f)
    return None

def build_teb_make(url: str):
    url_parsed = urlparse(url)
    f = url_parsed.path.split('/')[-1]
    path_to_f = PROJ_DIR / 'temp' / f
    path_to_f_dir = path_to_f.parent
    path_to_f_dir.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, path_to_f)
    path_to_f_dir_temp = path_to_f_dir / 'temp'
    with zipfile.ZipFile(path_to_f, 'r') as zip_ref:
        zip_ref.extractall(path_to_f_dir_temp)
    # look inside the temp to get the folder name
    path_to_temp_unzipped_dir = list(path_to_f_dir_temp.glob('*'))[0]
    path_to_unzipped_dir = PROJ_DIR / 'build' / path_to_temp_unzipped_dir.name
    shutil.move(str(path_to_temp_unzipped_dir), path_to_unzipped_dir)
    commit_id = str(path_to_unzipped_dir.name)
    path_to_exe_dir = path_to_unzipped_dir
    config_command = ['./mkmf.pl', '-t', 'gfortran_args',  '-p', 'driver.exe', 'src_driver', 'src_struct', 'src_proxi_SVAT', 'src_solar', 'src_teb']
    build_command = ['make']
    with open(path_to_exe_dir / str(commit_id + '_config.log'), 'w') as f:
        print('Configuring TEB for case: ' + commit_id)
        subprocess.check_call(config_command, cwd=path_to_exe_dir, stdout=f)
    with open(path_to_exe_dir / str(commit_id + '_build.log'), 'w') as f:
        print('Building TEB for case: ' + commit_id)
        subprocess.check_call(build_command, cwd=path_to_exe_dir, stdout=f)
    path_to_exe = path_to_exe_dir / 'driver.exe'
    return path_to_exe, commit_id

def checkout_run_load(commit_id: str, case_name: str, build_type: str, patch_nml: dict, download_zip) -> pd.DataFrame:
    if download_zip:
        path_to_exe, commit_id = build_teb_make(commit_id)
    else:
        git_checkout(commit_id)
        path_to_exe = build_teb(commit_id, build_type)

    path_to_case_dir = prepare_case(commit_id, case_name)
    run_teb(path_to_case_dir, path_to_exe, patch_nml)
    start, freq = get_date_params(path_to_case_dir / 'input.nml')
    df = load_txt(path_to_case_dir / 'output', start, freq, tz='UTC')
    return df

def compare(ref_id: str, trial_id: str, case_name: str, build_type: str, patch_nml=None, make=[False, False]) -> pd.DataFrame:
    df_ref = checkout_run_load(ref_id, case_name, build_type, patch_nml, make[0])
    df_trial = checkout_run_load(trial_id, case_name, build_type, patch_nml, make[1])
    # We have more outputs in the newer versions
    # drop any quantity not present in both versions.
    df_trial = df_trial[df_ref.columns]
    df_diff = df_ref - df_trial
    num_unequal_samples = len(df_diff[df_diff.values > 0])
    if num_unequal_samples != 0:
        raise RuntimeError(f"{num_unequal_samples} samples are not equal")
    else:
        print("All output samples are equals")
    return df_diff


if __name__ == "__main__":
    pass
