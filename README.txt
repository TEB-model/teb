                     *******************************
                     ** TEB simple driver version **
                     *******************************

                           METEO-FRANCE / CNRS

                            Version March 2018


contact : Valéry Masson (valery.masson@meteo.fr)


This simple driver for the TEB model (Building Energy Balance; 
Masson 2000 and subsequent papers), is intended to be used by scientists 
whishing to implement TEB in their own software. This simple driver 
provides all the source code of TEB alone, making it easier to analyse 
and to integrate in its own atmospheric model for example.

If you wish to use TEB for physical simulations, wihtout intention to extract 
all the TEB routines into another software environment, you could use the 
SURFEX platform. It contains TEB and much more, especially the ISBA scheme 
for vegetation, or several I/O formats, including Netcdf. 
(http://www.cnrm.meteo.fr/surfex/)

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
Documentation:
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

Please refer to the SURFEX website for the complete documentation on TEB:
http://www.umr-cnrm.fr/surfex/IMG/pdf/surfex_scidoc_v8.1.pdf

code browser:

There is a code browser available for all SURFEX code available from the SURFEX
page. The code browser corresponding to the physical routines of teb scheme
can be found starting from the teb_garden.F90 routine:

http://www.umr-cnrm.fr/surfex/data/BROWSER/doc_surf81/teb__garden_8F90.html

Note that, this browser is in the SURFEX environment. Because the garden and greenroofs
are replaced by proxy schemes in this offline version, you should not take attention
to the routines called by greenroof.F90, garden.F90 and teb_veg_properties.F90.

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
Physical description of TEB
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

The present simple "driver" version of TEB is based on :

------------------------------------------------------------
version8_1 of SURFEX code (r 5f2f00d8)
------------------------------------------------------------


The following physical features are included :

1) Historical version of TEB (Masson 2000) :
  - canyon shape geometry
  - isotropic direction of roads
  - 3 energy balances for roofs, roads and walls layers
  - simplified internal temperature evolution (force-restore)
  - water reservoirs on roofs and roads
  - snow mantel on roofs and roads

2) Building Energy Module (BEM) (Bueno et al 2012 doi:10.5194/gmd-5-433-2012)
  - adds Internal Building Energy Balance
  - adds floor and building mass energy balances
  - windows
  - Heat-Ventilation-air-Conditionning (HVAC) systems
  - inflitration / ventilation
  - solar protections
  - natural ventilation (opening of windows)

3) Roads orientation (Lemonsu et al 2012, DOI 10.1007/s10584-012-0521-6)
  - possibility to use a specified road orientation 
    (instead of a genericisotropic roads)
  - possibility to separate the 2 facing walls energy balances
    (then one have energy balances for walls 'A' and 'wall 'B')

4) Garden (Lemonsu et al 2012, DOI 10.1007/s10584-012-0521-6)
  - possibility to add a fraction of gardens inside the canyon itself
    (modifies all energetic exchanges within the canyon)
    ==> The garden needs its own vegetation energy balance. In this driver version,
        a simple proxi is used (based on fixed albedo and Bowen ratio)

5) Greenroofs (De Munck et al 2013)
  - possibility to add greenroofs on the roofs
    ==> The greenroof needs its own vegetation energy balance. In this driver version,
        a simple proxi is used (based on fixed albedo and Bowen ratio)

6) Irrigation and solar panels (DeMunck 2013 and Masson 2013)
  - possibility to have irrigation of greenroofs, gardens and watering of roads
  - possibility to have solar panels (hot water and/or photovoltaic)


Note: the only change of the TEB physical routines compared to the surfex v8.1 
version (open_surfex_v8_1.tar.gz (r 5f2f00d8) version ) are:

- a simpler ini_csts.F90 routine (less physical constants need to be initialized)
- the snow constants are initialized directly in module modd_snow_par.F90
- a correction of bug in teb.F90 routine for the extrapolation of air temperature 
  and humidity from forcing level to roof level.
- in teb_garden.F90, all references to garden and greenroofs (list of variables
  and arguments) are simplified, because of the use of proxys schemes for these.

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
Technical description of the driver structure
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

'src_teb' directory contains the physical sources of TEB
  ==> If you wish to implement TEB into your own software (e.g. atmospheric model)
      You only need to use these routines

'src_driver' directory contains the routines of the driver environment
  ==> Among these routines, you should only look at the driver.F90 program.
      This program is the main driver of this driver version of TEB.
      You can modify the values of the urban parameters (where indicated,
      from line 530 to line 865).

'src_struct' directory contains the modules defining the several fortran structures
      used to group the TEB variables (and hence limit the number of arguments in all
      teb subroutines). It also contains the teb_garden_struct.F90, bem_morpho_struct.F90 
      and window_data_struct.F90 called from the main program driver.F90. The  
      teb_garden_struct.F90 routine stores all the variables defined in the driver
      into the variables of the fortran structures, before calling the physical routine
      teb_garden.F90.   

'src_solar' directory contains a simplified sky model, that is in SURFEX, to
      estimate the part of the diffuse solar radiation that is near solar direction.
      This code is here only in order to retrieve exactly the same forcings data
      for TEB that what is in SURFEX (to insure that we have the same results).

'src_proxi-SVAT' directory contains extremely simple vegetation models in case you
     use the 'garden' option or the 'greenroof' option. 
  ==> These vegetation models are very crude, and you are invited to use your own
      vegetation models instead (or ISBA using SURFEX if you wish). Note that the
      physical varaibles that need to be exchanged between TEB and these vegetation
      model (either on roofs or gradens) are only those that are in arguments of
      the routines in 'src_proxi-SVAT'directory.

'input' directory contains data from the CAPITOUL experiment (Masson et al 2008),
      that is used as a test case.

'output' directory contains the results of your run once you have executed driver.exe

'output_ref' contains the results you should obtain, in order to check that everthing
       went well when you first ran driver.exe (do not change the urban parameters 
       and options in driver.F90 if you want to do this).
 
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
How to use the TEB driver program ?
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

1) You first need to create the makefile
   -------------------------------------

To do this, you have to:

a) define your fortran compiler options in file gfortran_args
b) run the script this way :

mkmf.pl -t gfortran_args -p driver.exe src_driver src_struct src_proxi_SVAT src_solar src_teb

c) You now should have a 'Makefile' file in your directory

Documentation on this mkmf.pl makefile creation tool is provided in mkmf.html


2) Compiling of the code
   ---------------------

Just use the command :
make

This will compile all the routines and create the executable file driver.exe

3) Running TEB !
   -----------

Just use the command :
driver.exe

The driver program should be running and show prints on your window like this:

17994/17999
17995/17999
17996/17999
17997/17999
17998/17999
17999/17999
  
     --------------------------
     |  DRIVER ENDS CORRECTLY |
     --------------------------

4) Looking at the results
   ----------------------

Just go to the 'output' directory. It contains the temporal evolution of
a few of TEB variables :

RN_TOWN.txt     : Net radiation      (W/m²)
H_TOWN.txt      : sensible heat flux (W/m²)
LE_TOWN.txt     : latent heat flux   (W/m²)

P_CANYON.txt    : pressure at road level      (Pa)
Q_CANYON.txt    : specific humidity in canyon (kg/kg)
T_CANYON.txt    : air temperature in canyon   (K)
U_CANYON.txt    : Wind strength in canyon     (m/s)

TI_BLD.txt      : Internal building temperature     (K)
T_ROAD1.txt     : road surface temperature          (K)
T_ROOF1.txt     : roof surface temperature          (K)
T_WALLA1.txt    : wall surface temperature (wall A) (K)
T_WALLB1.txt    : wall surface temperature (wall B) (K)


If BEM option is activated, you also have these outputs:
W / m² of buildings ground surface  (not by m² of floor)

HVAC_COOL.txt   : Cooling Energy Demand
HVAC_HEAT.txt   : Heating Energy Demand


5) Enjoy
   -----

You can now modify the urban options and parameters in driver.F90 in order
to see their impacts.


6) If you want to change the atmospheric forcing data
   --------------------------------------------------

a) change the characteristics of your forcing in driver.F90 (lines 97 to 107)
b) change the forcing data in the 'input' directory, but you should use the same
   type of text files
c) run driver.exe


7) If you want to implement TEB in your own software
   -------------------------------------------------

All the TEB sources you need are in the 'src_teb' directory.


--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
End. Have fun with TEB.
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
