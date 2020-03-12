# MinimalDX version 0.1.1 (https://www.github.com/dmey/minimal-dx).
# Copyright 2018 D. Meyer and R. Raustad. Licensed under MIT.

import subprocess
import os
from collections import namedtuple
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_style('ticks')

def main(mode, with_eplus_psychro=False):
    '''
    Compares the results from the simplified coil implementation and the original version in EnergyPlus
    '''

    # Check environment to determine folder structure
    if os.name == 'nt':
        path_to_exe = os.path.join('..', 'build', 'tests', mode, 'Debug')
    elif os.name == 'posix' or 'mac':
        path_to_exe = os.path.join('..', 'build', 'tests', mode)
    else:
        raise Exception('Error: this program cannot run on your system')

    # Define the path to the 2 programs called by run_tests.py
    path_to_make_test_data_program = os.path.join(path_to_exe, 'make_test_data_' + mode)
    if with_eplus_psychro:
        path_to_tests_program = os.path.join(path_to_exe, 'test_' + mode + '_eplus')
        eplus_suffix = ' with EPlus Psychro'
        eplus_suffix_folder = '-with-eplus-psychro'
    else:
        path_to_tests_program = os.path.join(path_to_exe, 'test_' + mode + '')
        eplus_suffix = ''
        eplus_suffix_folder = ''

    # Define the path to the output folder and check if exists else create it
    path_to_output_folder = os.path.join('outputs', mode + eplus_suffix_folder)
    if not os.path.exists(path_to_output_folder):
        os.makedirs(path_to_output_folder)

    # Define the path to the figure folder and check if exists else create it
    path_to_figures = os.path.join(path_to_output_folder, 'figures')
    if not os.path.exists(path_to_figures):
        os.makedirs(path_to_figures)

    # Create the test dataset
    path_to_test_dataset = os.path.join(path_to_output_folder, 'test_data_' + mode + '.csv')
    subprocess.check_call([path_to_make_test_data_program, path_to_test_dataset])

    # Read input dataset to use later
    input_dataset = pd.read_csv(path_to_test_dataset, delim_whitespace=True, na_values='Infinity')

    # Call implementations and read outputs
    Impl = namedtuple('Implementation', ['name', 'exe_path'])
    impls = [
        Impl('MinimalDX', os.path.abspath(path_to_tests_program))
        #Impl('SimDXCoolingCoil', os.path.abspath(path_to_tests_program))
    ]

    ref_impl = Impl('EnergyPlus', os.path.abspath(path_to_tests_program))

    outputs = {}
    for impl in [ref_impl] + impls:
        out_path = os.path.join(path_to_output_folder, impl.name + '.csv')
        with open('tests.log', 'w') as f:
            subprocess.check_call([impl.exe_path, path_to_test_dataset, out_path, impl.name], stdout=f)
            outputs[impl.name] = pd.read_csv(out_path, delim_whitespace=True, na_values='Infinity')

    # determine output variables and row count from ref impl
    Var = namedtuple('Variable', ['key', 'name', 'unit'])

    out_vars = [
        Var(name, *name.split('|'))
        for name in outputs[ref_impl.name].columns.values
    ]

    for out_var in out_vars:
        df = pd.DataFrame({
            impl.name: outputs[impl.name][out_var.key]
            for impl in [ref_impl] + impls
        })

        fig, ax = plt.subplots()
        sns.boxplot(data=df, orient="v", color='white', width=.25, ax=ax)
        plt.setp(ax.artists, edgecolor='k', facecolor='w')
        plt.setp(ax.lines, color='k')
        ax.set(xlabel='Model Name', ylabel=out_var.name + ' in ' + out_var.unit)
        sns.despine(offset=5, trim=True, ax=ax)
        fig.tight_layout()
        fig.savefig(os.path.join(path_to_figures, 'BoxPlot' + impl.name + out_var.name + '.png'))

        # Plot boxplot of relative error against reference implementation
        # Use absolute error metric to compare the new implementation with ASHRAE psychrometric
        # routines as we expect the differences to be large
        # Switch to relative error (expressed as a percentage of the reference implementation)
        # when comparing the new implementation that uses EnergyPlus psychrometric routines as
        # the two implementation should yield almost identical results
        # https://en.wikipedia.org/wiki/Approximation_error

        if impl.name == 'MinimalDX-ABS':
            error_type = ' absolute error '
            error_unit = out_var.unit
            error = pd.DataFrame({
                impl.name + ' - ' + ref_impl.name: (outputs[impl.name][out_var.key] - outputs[ref_impl.name][out_var.key])
                        for impl in impls

        })
        elif impl.name == 'MinimalDX':
            error_type = ' relative error '
            error_unit =  '%'
            error = pd.DataFrame({
                impl.name + ' - ' + ref_impl.name: ((outputs[impl.name][out_var.key] - outputs[ref_impl.name][out_var.key])  \
                                                    / outputs[ref_impl.name][out_var.key])*100
                        for impl in impls
        })

        max_diff = error.max()[0]
        if max_diff != 0.0:
            print('-------------------------- Max Diff----------------------')
            print('max_diff: ', error.max())
            print('\n')
            idxmax = error.idxmax()
            print('Input data that produce max diff', input_dataset.iloc[idxmax])
            print('\n\n')

        min_diff = error.min()[0]
        if min_diff != 0.0:
            print('-------------------------- Min Diff----------------------')
            print('min_diff: ', error.min())
            print('\n')
            idxmin = error.idxmin()
            print('Input data that produce min diff: ', input_dataset.iloc[idxmin])
            print('\n\n')

        fig, ax = plt.subplots()
        sns.boxplot(data=error, orient="v", color='white', width=.25, ax=ax)
        plt.setp(ax.artists, edgecolor='k', facecolor='w')
        plt.setp(ax.lines, color='k')
        ax.set(xlabel='Model Name', ylabel=out_var.name + error_type + 'in ' + error_unit)
        sns.despine(offset=5, trim=True, ax=ax)
        fig.tight_layout()
        fig.savefig(os.path.join(path_to_figures, 'BoxPlotError' + impl.name + out_var.name + '.png'))

        fig, ax = plt.subplots()
        ax.scatter(error.index.values, error, color='k', s=0.01)
        ax.set_xlabel('Run Number in 1')
        ax.set_ylabel(out_var.name + error_type + 'in '+ error_unit)
        ax.set_ylim([ min_diff, max_diff])
        sns.despine(offset=5, trim=True, ax=ax)
        fig.tight_layout()
        fig.savefig(os.path.join(path_to_figures, 'ScatterPlotError' + impl.name + out_var.name + '.png'))

        # Free memory at the end of each iteration
        plt.close()

if __name__ == "__main__":
    main('cooling')
    main('cooling', with_eplus_psychro=True)