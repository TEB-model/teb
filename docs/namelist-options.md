# Namelist options

List of options currently available in TEB.

## `&dimensions`

| Name           | Unit | Description            | Example |
| -------------- | ---- | ---------------------- | ------- |
| `NROOF_LAYER`  | -    | Number of roof layers  | `5`     |
| `NROAD_LAYER`  | -    | Number of road layers  | `5`     |
| `NWALL_LAYER`  | -    | Number of wall layers  | `5`     |
| `NFLOOR_LAYER` | -    | Number of floor layers | `5`     |


## `&parameters`

| Name            | Option          | Unit            | Description                                                                                          | Example    |
| --------------- | --------------- | --------------- | ---------------------------------------------------------------------------------------------------- | ---------- |
| `XTSTEP_SURF`   | -               | s               | Time step used in TEB (s)                                                                            | `300.`     |
| `IYEAR`         | -               | Year            | Current year (UTC)                                                                                   | `2004`     |
| `IMONTH`        | -               | Month           | Current month (UTC)                                                                                  | `2`        |
| `IDAY`          | -               | Day             | Current day (UTC)                                                                                    | `20`       |
| `ZTIME_START`   | -               | s               | Time at start of the run (UTC) (s)                                                                   | `82800.`   |
| `ZLON`          | -               | deg             | Longitude                                                                                            | `1.3`      |
| `ZLAT`          | -               | deg             | Latitude                                                                                             | `43.484`   |
| `INB_STEP_ATM`  | -               | s               | Number of time steps in input data                                                                   | `17999`    |
| `INB_ATM`       | -               | -               | Number of TEB calls within one atmospheric forcing time step. Forcing time step=XTSTEP_SURF*INB_ATM. | `6`        |
| `KSW`           | -               | -               | Number of spectral bands in SW forcing                                                               | `1`        |
| `ZZREF`         | -               | m               | Atmospheric forcing height above roof level                                                          | `28.`      |
| `NO_OUTPUT`     | `{TRUE, FALSE}` | -               | Whether to output results to file                                                                    | `.FALSE.`  |
| `CBEM`          |                 | -               | Building energy Model                                                                                | `DEF`      |
|                 | `DEF`           | -               | No Building Energy Model (BEM)                                                                       | -          |
|                 | `BEM`           | -               | Building Energy Model (BEM)                                                                          | -          |
| `HROAD_DIR`     |                 | -               | Road direction                                                                                       | `UNIF`     |
|                 | `UNIF`          | -               | Uniform roads                                                                                        | -          |
|                 | `ORIE`          | -               | Specified road orientation                                                                           | -          |
| `HWALL_OPT`     |                 | -               | Building energy Model                                                                                | `UNIF`     |
|                 | `UNIF`          | -               | Uniform walls                                                                                        | -          |
|                 | `TWO`           | -               | Two opposite  wall                                                                                   | -          |
| `LGARDEN`       | `{TRUE, FALSE}` | -               | Whether to activate garden module                                                                    | `.TRUE.`   |
| `LGREENROOF`    | `{TRUE, FALSE}` | -               | Whether to activate green roof module                                                                | `.TRUE.`   |
| `VEG_ALB`       | -               | -               | Garden and greenroof surface albedo                                                                  | `0.15`     |
| `VEG_EMIS`      | -               | -               | Garden and greenroof surface emissivity                                                              | `0.98`     |
| `GARDEN_BR`     | -               | -               | Bowen ratio for gardens                                                                              | `0.25`     |
| `GREENROOF_BR`  | -               | -               | Bowen ratio for greenroof                                                                            | `0.50`     |
| `ZFRAC_GR`      | -               | -               | Green roof fraction                                                                                  | `0.3`      |
| `LSOLAR_PANEL`  | `{TRUE, FALSE}` | -               | Whether to activate solar panel module                                                               | `.TRUE.`   |
| `ZFRAC_PANEL`   | -               | -               | Solar panels fraction on roofs                                                                       | `0.5`      |
| `LPAR_RD_IRRIG` | `{TRUE, FALSE}` | -               | Whether to activate road watering                                                                    | `.FALSE.`  |
| `LSOLAR_PANEL`  | `{TRUE, FALSE}` | -               | Whether to activate solar panel module                                                               | `.TRUE.`   |
| `CCOOL_COIL`    |                 | -               | Cooling system                                                                                       | `MINIDX`   |
|                 | `IDEAL`         | -               | Idealized                                                                                            | -          |
|                 | `DXCOIL`        | -               | Bueno et al., 2014                                                                                   | -          |
|                 | `MINIDX`        | -               | Meyer and Raustad, 2019                                                                              | -          |
| `CHEAT_COIL`    |                 | -               | Heating system                                                                                       | `IDEAL`    |
|                 | `IDEAL`         | -               | Idealized                                                                                            | -          |
|                 | `FINCAP`        | -               | Fincap                                                                                               | -          |
| `HZ0H`          |                 | -               | `z0h` roof & road:                                                                                   | `KAND07`   |
|                 | `MASC95`        | -               | Mascart et al, 1995                                                                                  | -          |
|                 | `BRUT82`        | -               | Brustaert, 1982                                                                                      | -          |
|                 | `KAND07`        | -               | Kanda, 2007                                                                                          | -          |
| `ZZ0`           | -               | m               | Roughness length                                                                                     | `2.`       |
| `ZBLD`          | -               | -               | Horizontal building area density                                                                     | `0.62`     |
| `ZGARDEN`       | -               | -               | Fraction of garden areas                                                                             | `0.11`     |
| `ZBLD_HEIGHT`   | -               | m               | Canyon height                                                                                        | `20.`      |
| `ZWALL_O_HOR`   | -               | -               | Vertical to horizontal surf ratio                                                                    | `1.05`     |
| `ZROAD_DIR`     | -               | deg             | N-S road  (deg from North, clockwise)                                                                | `0.`       |
| `ZALB_ROOF`     | -               | -               | Roof albedo                                                                                          | `0.40`     |
| `ZEMIS_ROOF`    | -               | -               | Roof emissivity                                                                                      | `0.97`     |
| `ZHC_ROOF`      | -               | J m-3 K-1       | Roof volumetric heat capacity (`NROOF_LAYER`)                                                        | `1580000.` |
| `ZTC_ROOF`      | -               | W m-1 K-1       | Roof thermal conductivity (`NROOF_LAYER`)                                                            | `1.15`     |
| `ZD_ROOF`       | -               | m               | Roof thickness  (`NROOF_LAYER`)                                                                      | `0.001`    |
| `ZALB_ROAD`     | -               | -               | Road albedo                                                                                          | `0.08`     |
| `ZEMIS_ROAD`    | -               | -               | Road emissivity                                                                                      | `0.97`     |
| `ZHC_ROAD`      | -               | J m-3 K-1       | Road volumetric heat capacity (`NROAD_LAYER`)                                                        | `1580000.` |
| `ZTC_ROAD`      | -               | W m-1 K-1       | Road thermal conductivity (`NROAD_LAYER`)                                                            | `1.15`     |
| `ZD_ROAD`       | -               | m               | Road thickness (`NROAD_LAYER`)                                                                       | `0.001`    |
| `ZALB_WALL`     | -               | -               | Wall albedo                                                                                          | `0.40`     |
| `ZEMIS_WALL`    | -               | -               | Wall emissivity                                                                                      | `0.97`     |
| `ZHC_WALL`      | -               | J m-3 K-1       | Wall volumetric heat capacity (`NWALL_LAYER`)                                                        | `1580000.` |
| `ZTC_WALL`      | -               | W m-1 K-1       | Wall thermal conductivity (`NWALL_LAYER`)                                                            | `1.15`     |
| `ZD_WALL`       | -               | m               | Wall thickness (`NWALL_LAYER`)                                                                       | `0.001`    |
| `ZHC_FLOOR`     | -               | J m-3 K-1       | Floor volumetric heat capacity (`NFLOOR_LAYER`)                                                      | `1580000.` |
| `ZTC_FLOOR`     | -               | W m-1 K-1       | Floor thermal conductivity (`NFLOOR_LAYER`)                                                          | `1.15`     |
| `ZD_FLOOR`      | -               | m               | Floor thickness (`NFLOOR_LAYER`)                                                                     | `0.001`    |
| `ZH_TRAFFIC`    | -               | W m-2           | Sensible heat fluxes due to traffic                                                                  | `0.`       |
| `ZLE_TRAFFIC`   | -               | W m-2           | Latent heat fluxes due to traffic                                                                    | `0.`       |
| `ZH_INDUSTRY`   | -               | W m-2           | Sensible heat fluxes due to factories                                                                | `0.`       |
| `ZLE_INDUSTRY`  | -               | W m-2           | Latent heat fluxes due to factories                                                                  | `0.`       |
| `ZEMIS_PANEL`   | -               | W m-2           | Emissivity of solar panel                                                                            | `0.9`      |
| `ZALB_PANEL`    | -               | -               | Albedo of solar panel                                                                                | `0.1`      |
| `ZEFF_PANEL`    | -               | -               | Efficiency of solar panel                                                                            | `0.14`     |
| `ZRESIDENTIAL`  | -               | -               | Fraction of residential use in buildings                                                             | `0.3`      |
| `ZDT_RES`       | -               | K               | Target temperature change when unoccupied (residential)                                              | `3.`       |
| `ZDT_OFF`       | -               | K               | Target temperature change when unoccupied (office)                                                   | `3.`       |
| `PFLOOR_HEIGHT` | -               | m               | Floor height                                                                                         | `3.`       |
| `ZINF`          | -               | AC/H            | Infiltration flow rate                                                                               | `0.5`      |
| `ZQIN`          | -               | W m-2(floor)    | Internal heat gains                                                                                  | `5.8`      |
| `ZQIN_FRAD`     | -               | -               | Radiant fraction of internal heat gains                                                              | `0.2`      |
| `ZQIN_FLAT`     | -               | -               | Latent fraction of internal heat gains                                                               | `0.2`      |
| `ZGR`           | -               | -               | Glazing ratio                                                                                        | `0.1`      |
| `ZSHGC`         | -               | -               | Window solar transmittance                                                                           | `0.763`    |
| `PU_WIN`        | -               | W m-2 K-1       | Window glass-to-glass U-factor                                                                       | `2.716`    |
| `LSHADE`        | `{TRUE, FALSE}` | -               | Whether to turn on shading devices                                                                   | `.FALSE.`  |
| `ZSHADE`        | `{0., 1.}`      | -               | Flag to activate shading devices                                                                     | `0.`       |
| `ZSHGC_SH`      | -               | -               | Window + shading solar heat gain coefficient                                                         | `0.025`    |
| `ZNATVENT`      | `{0., 1., 2.}`  | -               | Natural ventilation; none, manual, auto.                                                             | `0.`       |
| `ZV_VENT`       | -               | AC/H            | Ventilation flow rate                                                                                | `0.`       |
| `ZV_VENT`       | -               | -               | Fraction of evaporation for the condensers                                                           | `0.`       |
| `ZV_VENT`       | -               | -               | Fraction of waste heat released into the canyon                                                      | `1.`       |
| `ZTCOOL_TARGET` | -               | K               | Cooling setpoint of HVAC system                                                                      | `297.16`   |
| `ZTHEAT_TARGET` | -               | K               | Heating setpoint of HVAC system                                                                      | `297.16`   |
| `ZHR_TARGET`    | -               | -               | Relative humidity setpoint                                                                           | `0.5`      |
| `ZEFF_HEAT`     | -               | -               | Efficiency of the heating system                                                                     | `0.9`      |
| `ZCAP_SYS_HEAT` | -               | W m-2(bld)      | Capacity of the heating system                                                                       | `90.`      |
| `ZCAP_SYS_RAT`  | -               | W m-2(bld)      | Rated capacity of the cooling system                                                                 | `100.`     |
| `ZT_ADP`        | -               | K               | Apparatus dewpoint temperature of the cooling coil                                                   | `285.66`   |
| `ZM_SYS_RAT`    | -               | kg s-1 m-2(bld) | Rated HVAC mass flow rate                                                                            | `0.0067`   |
| `ZCOP_RAT`      | -               | -               | Rated COP of the cooling system                                                                      | `2.5`      |
| `CCH_BEM`       | -               | -               | TEB option for building outside convection. coeff.                                                   | `DOE-2`    |
|                 | `ROW30`         | -               | Default                                                                                              | -          |
|                 | `DOE-2`         | -               | DOE                                                                                                  | -          |
| `ZROUGH_ROOF`   | -               | -               | Roof roughness coeff. (`DOE-2`)                                                                      | `1.52`     |
| `ZROUGH_WALL`   | -               | -               | Wall roughness coeff. (`DOE-2`)                                                                      | `1.52`     |
| `ZT_ROAD`       | -               | K               | Initial road layers temperatures (`NROAD_LAYER`)                                                     | `274.`     |
| `ZT_ROOF`       | -               | K               | Initial roof layers temperatures (`NROOF_LAYER`)                                                     | `274.`     |
| `ZT_WALL`     | -               | K               | Initial wall layers temperatures (`NWALL_LAYER`)                                                     | `274.`     |
| `ZT_FLOOR`      | -               | K               | Initial building floor layers temperatures (`NFLOOR_LAYER`)                                          | `274.`     |
| `ZT_MASS`       | -               | K               | Initial building mass layers temperatures (`NFLOOR_LAYER`)                                           | `274.`     |
| `ZTI_BLD`       | -               | K               | Initial indoor air temperature                                                                       | `292.15`   |
| `ZT_CANYON`     | -               | K               | Initial canyon air temperature                                                                       | `274.`     |
| `ZT_WIN1`       | -               | K               | Initial external window temperature                                                                  | `275.`     |
| `ZT_WIN2`       | -               | K               | Initial internal window temperature                                                                  | `292.`     |
| `ZQ_CANYON`     | -               | kg kg-1         | Initial outdoor specific humidity                                                                    | `0.002`    |
| `ZQI_BLD`       | -               | kg kg-1         | Initial indoor air specific humidity                                                                 | `0.006`    |

