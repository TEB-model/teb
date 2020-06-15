<!-- omit in toc -->
# TEB Software documentation

- [Structure overview](#structure-overview)
  - [Source files](#source-files)
  - [Experiment files](#experiment-files)
    - [`input.nml`](#inputnml)
    - [`input`](#input)
    - [`output`](#output)
- [How to run TEB](#how-to-run-teb)

## Structure overview

Source and experiment files are organized in the `source` and `experiment` directories respectively. Please see below for more information. For a quick overview of the source file structure you can visit the code browser at [https://teb-model.github.io/teb](https://teb-model.github.io/teb).

### Source files

`src/teb`: directory contains the physical sources of TEB: If you wish to implement TEB into your own software (e.g. atmospheric model). You only need to use these routines.

`src/driver`: directory contains the routines of the driver environment. Among these routines, you should only look at the driver.F90 program. This program is the main driver of this driver version of TEB. You can modify the values of the urban parameters using the `input.nml` file.

`src/driver`: directory contains the modules defining the several fortran structures used to group the TEB variables (and hence limit the number of arguments in all teb subroutines). It also contains the `teb_garden_struct.F90`, `bem_morpho_struct.F90` and `window_data_struct.F90` called from the main program `driver.F90`. The `teb_garden_struct.F90` routine stores all the variables defined in the driver into the variables of the fortran structures, before calling the physical `routineteb_garden.F90`.

`src/solar`: directory contains a simplified sky model, that is in SURFEX, to estimate the part of the diffuse solar radiation that is near solar direction.This code is here only in order to retrieve exactly the same forcing data for TEB that what is in SURFEX (to insure that we have the same results).

`src/proxi-SVAT`: directory contains extremely simple vegetation models in case you use the 'garden' option or the 'greenroof' option. These vegetation models are very crude, and you are invited to use your own vegetation models instead (or ISBA using SURFEX if you wish). Note that the physical variables that need to be exchanged between TEB and these vegetation model (either on roofs or gardens) are only those that are in arguments of the routines in `src/proxi-SVAT` directory.

### Experiment files

The general structure to run an experiment in TEB is as follows: the user provides a configuration file named `input.nml` and a folder named `input` containing a set of input files. TEB simulates and writes output files in the `output` folder. Each input/output file represent a quantity with rows representing  individual time steps. The folder `examples/CAPITOUL` contains data used to set-up the CAPITOUL experiment (Masson et al 2008) and can be used as a template to start new experiments.


#### `input.nml`

This file contains the configuration for parameters used in the model run. For all the options currently available in TEB, please refer to [namelist options](namelist-options.md).


#### `input`

The following input files are required to run TEB. All quantities must be adjusted at forcing level.

| Name            | Description                                     | Unit       |
| --------------- | ----------------------------------------------- | ---------- |
| Forc_TA.txt     | Dry-bulb temperature                            | K          |
| Forc_QA.txt     | Air specific humidity                           | kg/m³      |
| Forc_PS.txt     | Atmospheric pressure                            | Pa         |
| Forc_WIND.txt   | Wind speed                                      | m/s        |
| Forc_DIR.txt    | Wind direction                                  | arc degree |
| Forc_DIR_SW.txt | Direct downwelling shortwave solar radiation    | W/m²       |
| Forc_SCA_SW.txt | Scattered downwelling shortwave solar radiation | W/m²       |
| Forc_LW.txt     | Downwelling longwave solar radiation            | W/m²       |
| Forc_RAIN.txt   | Liquid precipitation                            | kg/m²/s    |
| Forc_SNOW.txt   | Snow precipitation                              | kg/m²/s    |
| Forc_CO2.txt    | CO2 concentration                               | kg/m³      |


#### `output`

The following files are produced by the TEB `driver`.

| Name         | Description                                      | Unit  |
| ------------ | ------------------------------------------------ | ----- |
| RN_TOWN.txt  | Net all-wave radiation flux density              | W/m²  |
| H_TOWN.txt   | Turbulent sensible heat flux density             | W/m²  |
| LE_TOWN.txt  | Turbulent latent heat flux density               | W/m²  |
| P_CANYON.txt | Pressure at road level                           | Pa    |
| Q_CANYON.txt | Specific humidity at half building height        | kg/kg |
| T_CANYON.txt | Dry-bulb air temperature at half building height | K     |
| U_CANYON.txt | Wind speed at half building height               | m/s   |
| TI_BLD.txt   | Internal building temperature                    | K     |
| T_ROAD1.txt  | Road surface temperature                         | K     |
| T_ROOF1.txt  | Roof surface temperature                         | K     |
| T_WALLA1.txt | Wall surface temperature (wall A)                | K     |
| T_WALLB1.txt | Wall surface temperature (wall B)                | K     |


If building energy model option is activated (`CBEM = "BEM"`) in `input.nml`, the following outputs will also be present in W/m² of buildings ground surface  (not by m² of floor).

| Name          | Description          | Unit          |
| ------------- | -------------------- | ------------- |
| HVAC_COOL.txt | Cooling power demand | W/m² building |
| HVAC_HEAT.txt | Heating power demand | W/m² building |

If the solar panel option is activated (`LSOLAR_PANEL = .TRUE.`) in `input.nml`, the following outputs will also be present in W/m² panel.

| Name                | Description                                 | Unit       |
| ------------------- | ------------------------------------------- | ---------- |
| THER_PROD_PANEL.txt | Thermal power production of solar panels    | W/m² panel |
| PHOT_PROD_PANEL.txt | Electrical power production of solar panels | W/m² panel |


## How to run TEB

After building TEB as described in [README.md#installing](../README.md#installing) you can manually run the teb `driver` program from your experiment directory (e.g. `CAPITOUL`) or use the helper functions as described in the [tutorial](tutorial.ipynb). To manually run the `driver` program open you command prompt and type:

```
cd examples/CAPITOUL
../../build/driver
```

This should print the following to standard output:

```
17994/17999
17995/17999
17996/17999
17997/17999
17998/17999
17999/17999

     --------------------------
     |  DRIVER ENDS CORRECTLY |
     --------------------------
```

Results are available in the `output` directory and contain the temporal evolution of several TEB variables such as:

