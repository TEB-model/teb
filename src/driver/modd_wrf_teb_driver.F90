! TEB model version 4.0.0 (https://github.com/dmey/teb-model).
! Copyright 2018 D. Meyer. Licensed under CeCILL-C version 1.

MODULE MODD_WRF_TEB_DRIVER

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: TEB_DRIVER

  CONTAINS

    SUBROUTINE TEB_DRIVER( &
                     ! inputs
                     num_roof_layers, num_wall_layers, num_road_layers,       &
                     num_snow_layers, num_floor_layers,                       &
                     ZLON, ZLAT,                                              &
                     IYEAR, IMONTH, IDAY, ZTIME_START, XTSTEP_SURF,           &
                     XPS, XPA,                                                &
                     XTA, ZQA_KGKG, XRHOA, XCO2,                              &
                     XLW,                                                     &
                     XRAIN, XSNOW,                                            &
                     ZZREF,                                                   &
                     ZDIR_SW, ZSCA_SW, ZWIND, ZDIR,                           &
                     ! inputs / outputs
                     ZT_CANYON, ZQ_CANYON,                                    &
                     ZTI_BLD,                                                 &
                     ZT_ROOF, ZT_ROAD, ZT_WALL_A, ZT_WALL_B,                  &
                     ZWS_ROOF, ZWS_ROAD,                                      &
                     ZWSNOW_ROOF, ZTSNOW_ROOF, ZRSNOW_ROOF, ZASNOW_ROOF,      &
                     ZTSSNOW_ROOF, ZESNOW_ROOF,                               &
                     ZWSNOW_ROAD, ZTSNOW_ROAD, ZRSNOW_ROAD, ZASNOW_ROAD,      &
                     ZTSSNOW_ROAD, ZESNOW_ROAD,                               &
                     ZAUX_MAX, ZT_FLOOR,                                      &
                     ZT_MASS,                                                 &
                     ZT_WIN1, ZT_WIN2, ZQI_BLD,                               &
                     ZTHER_PRODC_DAY,                                         &
                     ! parameters
                     LGARDEN, LGREENROOF, LSOLAR_PANEL,                       &
                     HZ0H, HROAD_DIR, HWALL_OPT,                              &
                     HSNOW_ROAD,                                              &
                     HSNOW_ROOF,                                              &
                     ZH_TRAFFIC, ZLE_TRAFFIC, ZH_INDUSTRY, ZLE_INDUSTRY,      &
                     ZZ0,                                                     &
                     ZBLD, ZGARDEN, ZROAD_DIR, ZFRAC_GR,                      &
                     ZBLD_HEIGHT, ZWALL_O_HOR,                                &
                     ZALB_ROOF, ZEMIS_ROOF,                                   &
                     ZHC_ROOF, ZTC_ROOF, ZD_ROOF,                             &
                     ZALB_ROAD, ZEMIS_ROAD,                                   &
                     ZHC_ROAD, ZTC_ROAD,ZD_ROAD,                              &
                     ZALB_WALL, ZEMIS_WALL,                                   &
                     ZHC_WALL, ZTC_WALL, ZD_WALL,                             &
                     CCOOL_COIL, ZF_WATER_COND, CHEAT_COIL,                   &
                     HNATVENT, ZNATVENT,                                      &
                     XF_WASTE_CAN, ZQIN, ZQIN_FRAD,                           &
                     ZQIN_FLAT, ZGR, ZEFF_HEAT, ZINF,                         &
                     ZTCOOL_TARGET, ZTHEAT_TARGET, ZHR_TARGET,                &
                     ZV_VENT, ZCAP_SYS_HEAT, XCAP_SYS_RAT, ZT_ADP,            &
                     XM_SYS_RAT, ZCOP_RAT, ZHC_FLOOR, ZTC_FLOOR,              &
                     ZD_FLOOR,  ZSHGC,                                        &
                     ZSHGC_SH,                                                &
                     LSHADE, ZSHADE,                                          &
                     CBEM,                                                    &
                     CCH_BEM, ZROUGH_ROOF, ZROUGH_WALL,                       &
                     LPAR_RD_IRRIG, ZRD_START_MONTH, ZRD_END_MONTH,           &
                     ZRD_START_HOUR, ZRD_END_HOUR, ZRD_24H_IRRIG,             &
                     ZEMIS_PANEL, ZALB_PANEL, ZEFF_PANEL, ZFRAC_PANEL,        &
                     ZRESIDENTIAL,                                            &
                     ZDT_RES, ZDT_OFF,                                        &
                     ZFLOOR_HEIGHT, PU_WIN,                                   &
                     ! outputs
                     ZHVAC_COOL, ZHVAC_HEAT,                                  &
                     ZTHER_PROD_PANEL, ZPHOT_PROD_PANEL,                      &
                     XU_CANYON,                                               &
                     ZRN_TOWN, ZH_TOWN, ZLE_TOWN, ZGFLUX_TOWN, ZEVAP_TOWN,    &
                     ZUSTAR_TOWN,                                             &
                     ZTS_TOWN, ZEMIS_TOWN,                                    &
                     ZALB_TOWN, ZDIR_CANYON, ZQ_TOWN                          &
                     )

      ! TEB imports
      USE MODD_CSTS,     ONLY : XCPD, XPI, &
                                XRD, XP00
      USE MODD_SURF_ATM, ONLY: XCISMIN, XVMODMIN, LALDTHRES, XRIMAX
      USE MODD_SURF_PAR, ONLY: XUNDEF
      USE MODD_TYPE_DATE_SURF
      USE MODE_THERMOS
      USE MODD_REPROD_OPER, ONLY : CQSAT
      !
      USE MODI_SUNPOS
      USE MODI_OL_READ_ATM
      USE MODI_OL_ALLOC_ATM
      USE MODI_OL_TIME_INTERP_ATM
      USE MODI_TEB_GARDEN_STRUCT
      USE MODI_WINDOW_DATA_STRUCT
      USE MODI_BEM_MORPHO_STRUCT
      USE MODI_CIRCUMSOLAR_RAD

      IMPLICIT NONE

      INTEGER,      INTENT(IN)    :: num_roof_layers, &
                                     num_wall_layers, &
                                     num_road_layers, &
                                     num_snow_layers, &
                                     num_floor_layers

      ! Coordinates
      REAL, DIMENSION(1),                   INTENT(IN) :: ZLON                ! Longitude (deg)
      REAL, DIMENSION(1),                   INTENT(IN) :: ZLAT                ! Latitude (deg)
      REAL,                                 INTENT(IN) :: XTSTEP_SURF         ! length of timestep in seconds
      INTEGER,                              INTENT(IN) :: IYEAR, IMONTH, IDAY ! date at current time step
      REAL,                                 INTENT(IN) :: ZTIME_START         ! time at beginning of time step, since midnight (UTC, s)
      ! End coordinates

      ! Parameters
      LOGICAL,                              INTENT(IN) :: LSOLAR_PANEL    ! Flag to use a solar panels on roofs
      LOGICAL,                              INTENT(IN) :: LGARDEN         ! Flag to use a vegetation scheme in gardens
      LOGICAL,                              INTENT(IN) :: LGREENROOF      ! Flag to use a green roofs scheme
      CHARACTER(LEN=4),                     INTENT(IN) :: HROAD_DIR       ! road direction option :
                                                                          ! 'UNIF' : uniform roads
                                                                          ! 'ORIE' : specified road orientation
      CHARACTER(LEN=4),                     INTENT(IN) :: HWALL_OPT       ! Wall option
                                                                          ! 'UNIF' : uniform walls
                                                                          ! 'TWO ' : 2 opposite  walls
      REAL, DIMENSION(1),                   INTENT(IN) :: ZGARDEN         ! fraction of GARDEN areas
      REAL, DIMENSION(1),                   INTENT(IN) :: ZROAD_DIR       ! road direction (° from North, clockwise)
      REAL, DIMENSION(1),                   INTENT(IN) :: ZFRAC_GR        ! fraction of greenroofs on roofs
      CHARACTER(LEN=12),                    INTENT(IN) :: CCOOL_COIL      ! option for cooling device type
      CHARACTER(LEN=6),                     INTENT(IN) :: CHEAT_COIL      ! option for heating device type
      REAL, DIMENSION(1),                   INTENT(IN) :: ZF_WATER_COND   ! fraction of evaporation for the condensers
      CHARACTER(LEN=4), DIMENSION(1),       INTENT(IN) :: HNATVENT
      REAL, DIMENSION(1),                   INTENT(IN) :: ZNATVENT        ! flag to describe surventilation system for
                                                                          ! i/o 0 for NONE, 1 for MANU and 2 for AUTO
      REAL, DIMENSION(1),                   INTENT(IN) :: XF_WASTE_CAN    ! fraction of waste heat released into the canyon
      REAL, DIMENSION(1),                   INTENT(IN) :: ZQIN            ! Internal heat gains [W m-2(floor)]
      REAL, DIMENSION(1),                   INTENT(IN) :: ZQIN_FRAD       ! Radiant fraction of internal heat gains
      REAL, DIMENSION(1),                   INTENT(IN) :: ZQIN_FLAT       ! Latent franction of internal heat gains
      REAL, DIMENSION(1),                   INTENT(IN) :: ZGR             ! Glazing ratio
      REAL, DIMENSION(1),                   INTENT(IN) :: ZEFF_HEAT       ! Efficiency of the heating system
      REAL, DIMENSION(1),                   INTENT(IN) :: ZINF            ! Infiltration flow rate [AC/H]
      REAL, DIMENSION(1),                   INTENT(IN) :: ZTCOOL_TARGET   ! Cooling setpoint of HVAC system [K]
      REAL, DIMENSION(1),                   INTENT(IN) :: ZTHEAT_TARGET   ! Heating setpoint of HVAC system [K]
      REAL, DIMENSION(1),                   INTENT(IN) :: ZHR_TARGET      ! Relative humidity setpoint
      REAL, DIMENSION(1),                   INTENT(IN) :: ZV_VENT         ! Ventilation flow rate [AC/H]
      REAL, DIMENSION(1),                   INTENT(IN) :: ZCAP_SYS_HEAT   ! Capacity of the heating system [W m-2(bld)]
      REAL, DIMENSION(1),                   INTENT(IN) :: XCAP_SYS_RAT    ! Rated capacity of the cooling system [W m-2(bld)]
      REAL, DIMENSION(1),                   INTENT(IN) :: ZT_ADP          ! Apparatus dewpoint temperature of the cooling coil [K]
      REAL, DIMENSION(1),                   INTENT(IN) :: XM_SYS_RAT      ! Rated HVAC mass flow rate [kg s-1 m-2(bld)]
      REAL, DIMENSION(1),                   INTENT(IN) :: ZCOP_RAT        ! Rated COP of the cooling system
      REAL, DIMENSION(1, num_floor_layers), INTENT(IN) :: ZHC_FLOOR       ! heat capacity for floor layers
      REAL, DIMENSION(1, num_floor_layers), INTENT(IN) :: ZTC_FLOOR       ! thermal conductivity for floor layers
      REAL, DIMENSION(1, num_floor_layers), INTENT(IN) :: ZD_FLOOR        ! depth of floor layers
      REAL, DIMENSION(1),                   INTENT(IN) :: ZSHGC           ! window solar transmittance
      REAL, DIMENSION(1),                   INTENT(IN) :: ZSHGC_SH        ! window + shading solar heat gain coef.
      LOGICAL, DIMENSION(1),                INTENT(IN) :: LSHADE          ! Flag to use shading devices
      REAL, DIMENSION(1),                   INTENT(IN) :: PU_WIN          ! window U-factor [K m W-2]
      REAL, DIMENSION(1),                   INTENT(IN) :: ZSHADE          ! flag to activate shading devices -> REAL for i/o 0. or 1
      REAL, DIMENSION(1),                   INTENT(IN) :: ZFLOOR_HEIGHT   ! Floor height (m)
      CHARACTER(LEN=5),                     INTENT(IN) :: CCH_BEM         ! TEB option for building outside conv. coef
      CHARACTER(LEN=3),                     INTENT(IN) :: CBEM            ! Building Energy model 'DEF' or 'BEM'
      REAL, DIMENSION(1),                   INTENT(IN) :: ZROUGH_ROOF     ! roof roughness coef.
      REAL, DIMENSION(1),                   INTENT(IN) :: ZROUGH_WALL     ! wall roughness coef.
      LOGICAL,                              INTENT(IN) :: LPAR_RD_IRRIG   ! Flag for road watering
      REAL, DIMENSION(1),                   INTENT(IN) :: ZRD_START_MONTH ! start month for watering of roads(included)
      REAL, DIMENSION(1),                   INTENT(IN) :: ZRD_END_MONTH   ! end   month for watering of roads(included)
      REAL, DIMENSION(1),                   INTENT(IN) :: ZRD_START_HOUR  ! start hour  for watering of roads(included)
      REAL, DIMENSION(1),                   INTENT(IN) :: ZRD_END_HOUR    ! end   hour  for watering of roads(excluded)
      REAL, DIMENSION(1),                   INTENT(IN) :: ZRD_24H_IRRIG   ! 24h quantity of water used for road watering (liter/m2)
      REAL, DIMENSION(1),                   INTENT(IN) :: ZEMIS_PANEL     ! Emissivity of solar panel [-]
      REAL, DIMENSION(1),                   INTENT(IN) :: ZALB_PANEL      ! albedo of solar panel  [-]
      REAL, DIMENSION(1),                   INTENT(IN) :: ZEFF_PANEL      ! Efficiency of solar panel [-]
      REAL, DIMENSION(1),                   INTENT(IN) :: ZFRAC_PANEL     ! fraction of solar panels on roofs
      REAL, DIMENSION(1),                   INTENT(IN) :: ZRESIDENTIAL    ! Fraction of residential use in buildings
      REAL,                                 INTENT(IN) :: ZDT_RES         ! target temperature change when unoccupied
                                                                          ! (K) (residential buildings)
      REAL,                                 INTENT(IN) :: ZDT_OFF         ! target temperature change when unoccupied
                                                                          ! (K) (office buildings)
      ! End parameters

      ! Variable inputs
      REAL, DIMENSION(1),                  INTENT(IN) :: XCO2         ! CO2 concentration in the air    (kg/m3)
      CHARACTER(LEN=6),                    INTENT(IN) :: HZ0H         ! TEB option for z0h roof & road
                                                                      ! 'MASC95' : Mascart et al 1995
                                                                      ! 'BRUT82' : Brustaert     1982
                                                                      ! 'KAND07' : Kanda         2007
      REAL, DIMENSION(1),                  INTENT(IN) :: ZWIND        ! module of the horizontal wind
      CHARACTER(LEN=4),                    INTENT(IN) :: HSNOW_ROOF   ! snow roof scheme
                                                                      ! 'NONE'
                                                                      ! 'D95 '
                                                                      ! '1-L '
      CHARACTER(LEN=4),                    INTENT(IN) :: HSNOW_ROAD   ! snow road scheme
                                                                      ! 'NONE'
                                                                      ! 'D95 '
                                                                      ! '1-L '
      REAL, DIMENSION(1),                  INTENT(IN) :: XTA          ! air temperature forcing (K)
      REAL, DIMENSION(1),                  INTENT(IN) :: XPS          ! pressure at the surface
      REAL, DIMENSION(1),                  INTENT(IN) :: XPA          ! pressure at the first atmospheric level
      REAL, DIMENSION(1),                  INTENT(IN) :: XRHOA        ! air density at the lowest level
      REAL, DIMENSION(1),                  INTENT(IN) :: ZDIR         ! wind direction
      REAL, DIMENSION(1),                  INTENT(IN) :: XLW          ! atmospheric infrared radiation
      REAL, DIMENSION(1,1),                INTENT(IN) :: ZDIR_SW      ! incoming direct solar radiation on an horizontal surface
      REAL, DIMENSION(1,1),                INTENT(IN) :: ZSCA_SW      ! scattered incoming solar rad.
      REAL, DIMENSION(1),                  INTENT(IN) :: XRAIN        ! rain rate (kg/m2/s)
      REAL, DIMENSION(1),                  INTENT(IN) :: XSNOW        ! snow rate (kg/s/m2 of snow)
      REAL, DIMENSION(1),                  INTENT(IN) :: ZH_TRAFFIC   ! anthropogenic sensible heat fluxes due to traffic
      REAL, DIMENSION(1),                  INTENT(IN) :: ZLE_TRAFFIC  ! anthropogenic latent heat fluxes due to traffic
      REAL, DIMENSION(1),                  INTENT(IN) :: ZH_INDUSTRY  ! anthropogenic sensible heat fluxes due to factories
      REAL, DIMENSION(1),                  INTENT(IN) :: ZLE_INDUSTRY ! anthropogenic latent heat fluxes due to factories
      REAL, DIMENSION(1),                  INTENT(IN) :: ZZREF        ! reference height of the first atmospheric level (temperature)
      REAL, DIMENSION(1),                  INTENT(IN) :: ZZ0          ! town roughness length for momentum
      REAL, DIMENSION(1),                  INTENT(IN) :: ZBLD         ! fraction of buildings
      REAL, DIMENSION(1),                  INTENT(IN) :: ZBLD_HEIGHT  ! buildings h
      REAL, DIMENSION(1),                  INTENT(IN) :: ZWALL_O_HOR  ! wall surf. / hor. surf.
      REAL, DIMENSION(1),                  INTENT(IN) :: ZALB_ROOF    ! roof albedo
      REAL, DIMENSION(1),                  INTENT(IN) :: ZEMIS_ROOF   ! roof emissivity
      REAL, DIMENSION(1, num_roof_layers), INTENT(IN) :: ZHC_ROOF     ! heat capacity for roof layers (J/m3/K)
      REAL, DIMENSION(1, num_roof_layers), INTENT(IN) :: ZTC_ROOF     ! thermal conductivity for roof layers (w/m/K)
      REAL, DIMENSION(1, num_roof_layers), INTENT(IN) :: ZD_ROOF      ! depth of roof layers
      REAL, DIMENSION(1),                  INTENT(IN) :: ZALB_ROAD    ! road albedo
      REAL, DIMENSION(1),                  INTENT(IN) :: ZEMIS_ROAD   ! road emissivity
      REAL, DIMENSION(1, num_road_layers), INTENT(IN) :: ZHC_ROAD     ! heat capacity for road layers (J/m3/K)
      REAL, DIMENSION(1, num_road_layers), INTENT(IN) :: ZTC_ROAD     ! thermal conductivity for road layers (w/m/K)
      REAL, DIMENSION(1, num_road_layers), INTENT(IN) :: ZD_ROAD      ! depth of road layers
      REAL, DIMENSION(1),                  INTENT(IN) :: ZALB_WALL    ! wall albedo
      REAL, DIMENSION(1),                  INTENT(IN) :: ZEMIS_WALL   ! wall emissivity
      REAL, DIMENSION(1, num_wall_layers), INTENT(IN) :: ZHC_WALL     ! heat capacity for wall layers (J/m3/K)
      REAL, DIMENSION(1, num_wall_layers), INTENT(IN) :: ZTC_WALL     ! thermal conductivity for wall layers
      REAL, DIMENSION(1, num_wall_layers), INTENT(IN) :: ZD_WALL      ! depth of wall layers
      REAL, DIMENSION(1),                  INTENT(IN) :: ZQA_KGKG     ! air humidity at forcing level (kg/kg)
      ! End variable inputs

      ! State
      REAL, DIMENSION(1),                       INTENT(INOUT)   :: ZT_WIN1         ! outdoor window temperature [K]
      REAL, DIMENSION(1),                       INTENT(INOUT)   :: ZT_WIN2         ! Indoor window temperature [K]
      REAL, DIMENSION(1),                       INTENT(INOUT)   :: ZAUX_MAX        ! Auxiliar variable for autosize calcs (not used)
      REAL, DIMENSION(1),                       INTENT(INOUT)   :: ZTHER_PRODC_DAY ! Present day integrated thermal production of energy (J/m2 panel)
      REAL, DIMENSION(1),                       INTENT(INOUT)   :: ZQI_BLD         ! Indoor air specific humidity [kg kg-1]
      REAL, DIMENSION(1, num_floor_layers),     INTENT(INOUT)   :: ZT_FLOOR        ! Floor layers temperatures [K]
      REAL, DIMENSION(1, num_floor_layers),     INTENT(INOUT)   :: ZT_MASS         ! Internal mass layers temperatures [K]
      REAL, DIMENSION(1, num_road_layers),      INTENT(INOUT)   :: ZT_ROAD         ! road layers temperatures
      REAL, DIMENSION(1, num_roof_layers),      INTENT(INOUT)   :: ZT_ROOF         ! roof layers temperatures
      REAL, DIMENSION(1, num_wall_layers),      INTENT(INOUT)   :: ZT_WALL_A       ! wall layers temperatures (wall 'A')
      REAL, DIMENSION(1, num_wall_layers),      INTENT(INOUT)   :: ZT_WALL_B       ! wall layers temperatures (wall 'B')
      REAL, DIMENSION(1),                       INTENT(INOUT)   :: ZT_CANYON       ! canyon air temperature
      REAL, DIMENSION(1),                       INTENT(INOUT)   :: ZQ_CANYON       ! canyon air specific humidity
      REAL, DIMENSION(1),                       INTENT(INOUT)   :: ZTI_BLD         ! inside building temperature
      REAL, DIMENSION(1),                       INTENT(INOUT)   :: ZWS_ROOF        ! roof water reservoir(kg/m2)
      REAL, DIMENSION(1),                       INTENT(INOUT)   :: ZWS_ROAD        ! road water reservoir(kg/m2)
      REAL, DIMENSION(1, num_snow_layers),      INTENT(INOUT)   :: ZWSNOW_ROOF     ! snow layers reservoir
      REAL, DIMENSION(1, num_snow_layers),      INTENT(INOUT)   :: ZTSNOW_ROOF     ! snow layers temperature
      REAL, DIMENSION(1, num_snow_layers),      INTENT(INOUT)   :: ZRSNOW_ROOF     ! snow layers density
      REAL, DIMENSION(1),                       INTENT(INOUT)   :: ZASNOW_ROOF     ! snow albedo
      REAL, DIMENSION(1),                       INTENT(INOUT)   :: ZESNOW_ROOF     ! snow emissivity
      REAL, DIMENSION(1),                       INTENT(INOUT)   :: ZTSSNOW_ROOF    ! snow surface temperature
      REAL, DIMENSION(1, num_snow_layers),      INTENT(INOUT)   :: ZWSNOW_ROAD     ! snow layers reservoir
      REAL, DIMENSION(1, num_snow_layers),      INTENT(INOUT)   :: ZTSNOW_ROAD     ! snow layers temperature
      REAL, DIMENSION(1, num_snow_layers),      INTENT(INOUT)   :: ZRSNOW_ROAD     ! snow layers density
      REAL, DIMENSION(1),                       INTENT(INOUT)   :: ZASNOW_ROAD     ! snow albedo
      REAL, DIMENSION(1),                       INTENT(INOUT)   :: ZESNOW_ROAD     ! snow emissivity
      REAL, DIMENSION(1),                       INTENT(INOUT)   :: ZTSSNOW_ROAD    ! snow surface temperature
      ! End state

      ! Diagnostics
      REAL, DIMENSION(1),                   INTENT(OUT) :: ZHVAC_COOL    ! Energy consumption of the cooling system [W m-2(bld)]
      REAL, DIMENSION(1),                   INTENT(OUT) :: ZHVAC_HEAT    ! Energy consumption of the heating system [W m-2(bld)]
      REAL, DIMENSION(1),                   INTENT(OUT) :: ZTHER_PROD_PANEL ! Thermal energy production of solar panel on roofs (W/m2 panel)
      REAL, DIMENSION(1),                   INTENT(OUT) :: ZPHOT_PROD_PANEL ! Photovoltaic Energy production of solar panel on roofs (W/m2 panel)
      REAL, DIMENSION(1),                   INTENT(OUT) :: ZRN_TOWN      ! net radiation over town
      REAL, DIMENSION(1),                   INTENT(OUT) :: ZH_TOWN       ! sensible heat flux over town
      REAL, DIMENSION(1),                   INTENT(OUT) :: ZLE_TOWN      ! latent heat flux over town
      REAL, DIMENSION(1),                   INTENT(OUT) :: ZGFLUX_TOWN   ! flux through the ground
      REAL, DIMENSION(1),                   INTENT(OUT) :: ZEVAP_TOWN    ! evaporation flux (kg/m2/s)
      REAL, DIMENSION(1),                   INTENT(OUT) :: ZUSTAR_TOWN   ! friction velocity over town
      REAL, DIMENSION(1),                   INTENT(OUT) :: ZTS_TOWN      ! town surface temperature
      REAL, DIMENSION(1),                   INTENT(OUT) :: ZEMIS_TOWN    ! town equivalent emissivity
      REAL, DIMENSION(1),                   INTENT(OUT) :: ZALB_TOWN     ! town equivalent albedo
      REAL, DIMENSION(1),                   INTENT(OUT) :: XU_CANYON     ! canyon hor. wind
      REAL, DIMENSION(1),                   INTENT(OUT) :: ZDIR_CANYON   ! canyon wind direction
      REAL, DIMENSION(1),                   INTENT(OUT) :: ZQ_TOWN       ! Town humidity
      ! End diagnostics

      ! Constants
      CHARACTER(LEN=3)   :: HIMPLICIT_WIND = 'NEW'   ! Implicitation option for wind fluxes
      INTEGER            :: KSW            = 1       ! number of spectral bands in SW forcing
      LOGICAL            :: LCANOPY        = .FALSE. ! is canopy active ?
      REAL, DIMENSION(1) :: XSW_BANDS      = XUNDEF  ! mean wavelength of each shortwave band (m)
                                                     ! (unused in TEB, but needed as argument)
      ! End constants

      ! Local variables
      REAL                  :: ZTIME_BEG       ! Time at beginning of time step
      REAL                  :: ZTIME           ! Time at end       of time step
      TYPE(DATE_TIME)       :: TPTIME
      REAL, DIMENSION(1)    :: ZU_CANYON       ! canyon hor. wind
      REAL, DIMENSION(1)    :: ZU_LOWCAN       ! wind near the road
      REAL, DIMENSION(1)    :: ZT_LOWCAN       ! temp. near the road
      REAL, DIMENSION(1)    :: ZQ_LOWCAN       ! hum. near the road
      REAL, DIMENSION(1)    :: ZZ_LOWCAN       ! height of atm. var. near the road
      REAL, DIMENSION(1)    :: ZVMOD           ! module of the horizontal wind
      REAL, DIMENSION(1)    :: ZEXNA           ! exner function at the lowest level
      REAL, DIMENSION(1)    :: ZEXNS           ! surface exner function
      REAL, DIMENSION(1)    :: ZF1_o_B         ! Coefficient for sky model
      REAL, DIMENSION(1)    :: ZCOEF           ! work array
      REAL, DIMENSION(1)    :: XZENITH         ! solar zenithal angle (radian)
      REAL, DIMENSION(1)    :: XZENITH2        ! solar zenithal angle (radian)
      REAL, DIMENSION(1)    :: XTSUN           ! solar time (s from midnight)
      REAL, DIMENSION(1)    :: XAZIM           ! azimuthal angle (radian from North, clockwise)
      REAL, DIMENSION(1)    :: ZPEW_A_COEF_LOWCAN
      REAL, DIMENSION(1)    :: ZPEW_B_COEF_LOWCAN
      REAL, DIMENSION(1)    :: ZPEW_A_COEF     ! implicit coefficients
      REAL, DIMENSION(1)    :: ZPEW_B_COEF     ! for wind coupling
      ! V = ZPEW_A_COEF * v'w'  +   ZPEW_B_COEF
      ! IF explicit: ZPEW_A_COEF = 0 and  ZPEW_B_COEF
      REAL, DIMENSION(1)    :: ZWAKE           ! reduction of average wind speed
      REAL, DIMENSION(1)    :: ZSVF_WALL       ! wall sky view factor
      REAL, DIMENSION(1)    :: ZCAN_HW_RATIO   ! canyon    h/W
      REAL, DIMENSION(1)    :: ZSVF_ROAD       ! road sky view factor
      REAL, DIMENSION(1)    :: ZROAD           ! fraction of roads
      REAL, DIMENSION(1)    :: ZROAD_O_GRND    ! Road to ground surface ratio
      REAL, DIMENSION(1)    :: ZGARDEN_O_GRND  ! Garden to ground surface ratio
      REAL, DIMENSION(1)    :: ZWALL_O_GRND    ! Wall to ground surface ratio
      REAL, DIMENSION(1)    :: ZSVF_GARDEN     ! garden sky view factor
      REAL, DIMENSION(1)    :: ZTOTS_O_HORS    ! total canyon+roof surf over horizontal surf
      REAL, DIMENSION(1)    :: ZROOF_FRAC      ! roof area fraction
      REAL, DIMENSION(1)    :: ZWALL_FRAC      ! wall area fraction
      REAL, DIMENSION(1)    :: ZROAD_FRAC      ! road area fraction
      REAL, DIMENSION(1)    :: ZGARDEN_FRAC    ! garden area fraction
      REAL, DIMENSION(1)    :: ZUGG_WIN        ! window glass-to-glass U-factro [W m-2 K-1]
      REAL, DIMENSION(1)    :: ZALB_WIN        ! window albedo
      REAL, DIMENSION(1)    :: ZABS_WIN        ! window absortance
      REAL, DIMENSION(1)    :: ZTRAN_WIN       ! window transmittance
      REAL, DIMENSION(1)    :: ZN_FLOOR        ! Number of floors
      REAL, DIMENSION(1)    :: ZWALL_O_BLD     ! Wall area [m2_wall/m2_bld]
      REAL, DIMENSION(1)    :: ZGLAZ_O_BLD     ! Window area [m2_win/m2_bld]
      REAL, DIMENSION(1)    :: ZMASS_O_BLD     ! Mass area [m2_mass/m2_bld]
      REAL, DIMENSION(1)    :: ZFLOOR_HW_RATIO ! H/W ratio of 1 floor level
      REAL, DIMENSION(1)    :: ZF_FLOOR_MASS   ! View factor floor-mass
      REAL, DIMENSION(1)    :: ZF_FLOOR_WALL   ! View factor floor-wall
      REAL, DIMENSION(1)    :: ZF_FLOOR_WIN    ! View factor floor-window
      REAL, DIMENSION(1)    :: ZF_FLOOR_ROOF   ! View factor floor-roof
      REAL, DIMENSION(1)    :: ZF_WALL_FLOOR   ! View factor wall-floor
      REAL, DIMENSION(1)    :: ZF_WALL_MASS    ! View factor wall-mass
      REAL, DIMENSION(1)    :: ZF_WALL_WIN     ! View factor wall-win
      REAL, DIMENSION(1)    :: ZF_WIN_FLOOR    ! View factor win-floor
      REAL, DIMENSION(1)    :: ZF_WIN_MASS     ! View factor win-mass
      REAL, DIMENSION(1)    :: ZF_WIN_WALL     ! View factor win-wall
      REAL, DIMENSION(1)    :: ZF_MASS_FLOOR   ! View factor mass-floor
      REAL, DIMENSION(1)    :: ZF_MASS_WALL    ! View factor mass-wall
      REAL, DIMENSION(1)    :: ZF_MASS_WIN     ! View factor mass-window
      REAL, DIMENSION(1)    :: ZF_WIN_WIN      ! indoor win to win view factor
      REAL, DIMENSION(1)    :: ZDIR_ALB_TOWN   ! town equivalent direct albedo
      REAL, DIMENSION(1)    :: ZSCA_ALB_TOWN   ! town equivalent diffuse albedo
      REAL, DIMENSION(1,1)  :: ZTDIR_SW        ! total direct SW
      REAL, DIMENSION(1,1)  :: ZTSCA_SW        ! total diffuse SW
      ! TODO: what about these? make into state variables?
      ! For now set to off state as they are not considered in the eval.
      LOGICAL, DIMENSION(1) :: GSHAD_DAY       ! has shading been necessary this day ?
      LOGICAL, DIMENSION(1) :: GNATVENT_NIGHT  ! has natural ventilation been necessary/possible this night ?
      REAL, DIMENSION(1,1)  :: XDIR_SW         ! incoming direct solar radiation on an horizontal surface
      REAL, DIMENSION(1,1)  :: XSCA_SW         ! scattered incoming solar rad.
      REAL, DIMENSION(1)    :: PFLOOR_HEIGHT   ! Floor height (m)
      REAL, DIMENSION(1)    :: ZF_WASTE_CAN    ! fraction of waste heat released into the canyon
      REAL, DIMENSION(1)    :: ZCAP_SYS_RAT    ! Rated capacity of the cooling system [W m-2(bld)]
      REAL, DIMENSION(1)    :: ZM_SYS_RAT      ! Rated HVAC mass flow rate [kg s-1 m-2(bld)]
      ! End local variables

      ! Unused model outputs
      REAL, DIMENSION(1) :: ZRN_STRLROOF       ! net radiation over structural roof
      REAL, DIMENSION(1) :: ZH_STRLROOF        ! sensible heat flux over structural roof
      REAL, DIMENSION(1) :: ZLE_STRLROOF       ! latent heat flux over structural roof
      REAL, DIMENSION(1) :: ZGFLUX_STRLROOF    ! flux through the structural roof
      REAL, DIMENSION(1) :: ZRUNOFF_STRLROOF   ! water runoff on the structural roof
      REAL, DIMENSION(1) :: ZRN_BLT            ! net radiation over built covers
      REAL, DIMENSION(1) :: ZH_BLT             ! sensible heat flux over built covers
      REAL, DIMENSION(1) :: ZLE_BLT            ! latent heat flux over built covers
      REAL, DIMENSION(1) :: ZGFLUX_BLT         ! flux through the built covers
      REAL, DIMENSION(1) :: ZSFCO2             ! Surface Flux of CO2
      REAL, DIMENSION(1) :: ZUW_GRND           ! friction flux over ground
      REAL, DIMENSION(1) :: ZUW_ROOF           ! Momentum flux for roofs
      REAL, DIMENSION(1) :: ZDUWDU_GRND        ! d(u'w')/du for ground
      REAL, DIMENSION(1) :: ZDUWDU_ROOF        ! d(u'w')/du for roof
      REAL, DIMENSION(1) :: ZAC_GARDEN         ! garden aerodynamical conductance
      REAL, DIMENSION(1) :: ZAC_GREENROOF      ! green roofs aerodynamical conductance
      REAL, DIMENSION(1) :: ZAC_ROAD_WAT       ! road aerodynamical conductance (for water)
      REAL, DIMENSION(1) :: ZAC_GARDEN_WAT     ! garden aerodynamical conductance for vapor
      REAL, DIMENSION(1) :: ZAC_GREENROOF_WAT  ! green roofs aerodynamical conductance for vapor
      REAL, DIMENSION(1) :: ZABS_SW_ROOF       ! Shortwave radiation absorbed by roofs
      REAL, DIMENSION(1) :: ZABS_LW_ROOF       ! Longwave  radiation absorbed by roofs
      REAL, DIMENSION(1) :: ZABS_SW_SNOW_ROOF  ! Shortwave radiation absorbed by roof snow
      REAL, DIMENSION(1) :: ZABS_LW_SNOW_ROOF  ! abs. LW rad. by snow
      REAL, DIMENSION(1) :: ZABS_SW_ROAD       ! Shortwave radiation absorbed by roads
      REAL, DIMENSION(1) :: ZABS_LW_ROAD       ! Longwave  radiation absorbed by roads
      REAL, DIMENSION(1) :: ZABS_SW_SNOW_ROAD  ! Shortwave radiation absorbed by road snow
      REAL, DIMENSION(1) :: ZABS_LW_SNOW_ROAD  ! abs. LW rad. by snow
      REAL, DIMENSION(1) :: ZABS_SW_WALL_A     ! Shortwave radiation absorbed by wall A
      REAL, DIMENSION(1) :: ZABS_LW_WALL_A     ! Longwave  radiation absorbed by wall A
      REAL, DIMENSION(1) :: ZABS_SW_WALL_B     ! Shortwave radiation absorbed by wall B
      REAL, DIMENSION(1) :: ZABS_LW_WALL_B     ! Longwave  radiation absorbed by wall B
      REAL, DIMENSION(1) :: ZABS_SW_PANEL      ! Shortwave radiation absorbed by solar panels
      REAL, DIMENSION(1) :: ZABS_LW_PANEL      ! Longwave  radiation absorbed by solar panels
      REAL, DIMENSION(1) :: ZABS_SW_GARDEN     ! Shortwave radiation absorbed by gardens
      REAL, DIMENSION(1) :: ZABS_LW_GARDEN     ! Longwave  radiation absorbed by gardens
      REAL, DIMENSION(1) :: ZABS_SW_GREENROOF  ! Shortwave radiation absorbed by greenroofs
      REAL, DIMENSION(1) :: ZABS_LW_GREENROOF  ! Longwave  radiation absorbed by greenroofs
      REAL, DIMENSION(1) :: ZG_GREENROOF_ROOF  ! heat flux between base of greenroof and structural roof
      REAL, DIMENSION(1) :: ZH_BLD_COOL        ! Sensible cooling energy demand of the building [W m-2(bld)]
      REAL, DIMENSION(1) :: ZT_BLD_COOL        ! Total cooling energy demand of the building [W m-2(bld)]
      REAL, DIMENSION(1) :: ZH_BLD_HEAT        ! Heating energy demand of the building [W m-2(bld)]
      REAL, DIMENSION(1) :: ZLE_BLD_COOL       ! Latent cooling energy demand of the building [W m-2(bld)]
      REAL, DIMENSION(1) :: ZLE_BLD_HEAT       ! Latent heating energy demand of the building [W m-2(bld)]
      REAL, DIMENSION(1) :: ZH_WASTE           ! Sensible waste heat from HVAC system [W m-2(tot)]
      REAL, DIMENSION(1) :: ZLE_WASTE          ! Latent waste heat from HVAC system [W m-2(tot)]
      REAL, DIMENSION(1) :: ZCAP_SYS           ! Actual capacity of the cooling system [W m-2(bld)]
      REAL, DIMENSION(1) :: ZM_SYS             ! Actual HVAC mass flow rate [kg s-1 m-2(bld)]
      REAL, DIMENSION(1) :: ZCOP               ! COP of the cooling system
      REAL, DIMENSION(1) :: ZQ_SYS             ! Supply air specific humidity [kg kg-1]
      REAL, DIMENSION(1) :: ZT_SYS             ! Supply air temperature [K]
      REAL, DIMENSION(1) :: ZTR_SW_WIN         ! Solar radiation transmitted throught windows [W m-2(bld)]
      REAL, DIMENSION(1) :: ZFAN_POWER         ! HVAC fan power
      REAL, DIMENSION(1) :: ZABS_SW_WIN        ! window absorbed shortwave radiation [W m-2]
      REAL, DIMENSION(1) :: ZABS_LW_WIN        ! absorbed infrared rad. [W m-2]
      REAL, DIMENSION(1) :: ZEMIT_LW_FAC       ! LW flux emitted by the facade (W/m² ground)
      REAL, DIMENSION(1) :: ZEMIT_LW_GRND      ! LW flux emitted by the ground (W/m² ground)
      REAL, DIMENSION(1) :: ZT_RAD_IND         ! Indoor mean radiant temperature [K]
      REAL, DIMENSION(1) :: ZREF_SW_GRND       ! total solar rad reflected from ground
      REAL, DIMENSION(1) :: ZREF_SW_FAC        ! total solar rad reflected from facade
      REAL, DIMENSION(1) :: ZHU_BLD            ! Indoor relative humidity 0 < (-) < 1
      REAL, DIMENSION(1) :: ZIRRIG_ROAD        ! road irrigation during current time-step
      REAL, DIMENSION(1) :: ZPROD_PANEL        ! Averaged energy production of solar panel on roofs (W/m2 panel)
      REAL, DIMENSION(1) :: ZTHER_PROD_BLD     ! Thermal energy production of solar panel on roofs (W/m2 bld)
      REAL, DIMENSION(1) :: ZPHOT_PROD_BLD     ! Photovoltaic energy production of solar panel on roofs (W/m2 bld)
      REAL, DIMENSION(1) :: ZPROD_BLD          ! Averaged energy production of solar panel on roofs (W/m2 bld)
      REAL, DIMENSION(1) :: ZH_PANEL           ! Sensible heat flux from solar panels (W/m2 panel)
      REAL, DIMENSION(1) :: ZRN_PANEL          ! Net radiation of solar panel (W/m2 panel)
      REAL, DIMENSION(1) :: ZCUR_QIN           ! Internal heat gains [W m-2(floor)]
      REAL, DIMENSION(1) :: ZCUR_TCOOL_TARGET  ! Cooling setpoint of HVAC system [K]
      REAL, DIMENSION(1) :: ZCUR_THEAT_TARGET  ! Heating setpoint of HVAC system [K]
      REAL, DIMENSION(1) :: ZRN_ROOF           ! net radiation over roof
      REAL, DIMENSION(1) :: ZH_ROOF            ! sensible heat flux over roof
      REAL, DIMENSION(1) :: ZLE_ROOF           ! latent heat flux over roof
      REAL, DIMENSION(1) :: ZLEW_ROOF          ! latent heat flux over roof (snow)
      REAL, DIMENSION(1) :: ZGFLUX_ROOF        ! flux through the roof
      REAL, DIMENSION(1) :: ZRUNOFF_ROOF       ! runoff over the ground
      REAL, DIMENSION(1) :: ZRN_ROAD           ! net radiation over road
      REAL, DIMENSION(1) :: ZH_ROAD            ! sensible heat flux over road
      REAL, DIMENSION(1) :: ZLE_ROAD           ! latent heat flux over road
      REAL, DIMENSION(1) :: ZLEW_ROAD          ! latent heat flux over road (snow)  (sublimation)
      REAL, DIMENSION(1) :: ZGFLUX_ROAD        ! flux through the road
      REAL, DIMENSION(1) :: ZRUNOFF_ROAD       ! runoff over the ground
      REAL, DIMENSION(1) :: ZRN_WALL_A         ! net radiation over wall
      REAL, DIMENSION(1) :: ZH_WALL_A          ! sensible heat flux over wall
      REAL, DIMENSION(1) :: ZLE_WALL_A         ! latent heat flux over wall
      REAL, DIMENSION(1) :: ZGFLUX_WALL_A      ! flux through the wall
      REAL, DIMENSION(1) :: ZRN_WALL_B         ! net radiation over wall
      REAL, DIMENSION(1) :: ZH_WALL_B          ! sensible heat flux over wall
      REAL, DIMENSION(1) :: ZLE_WALL_B         ! latent heat flux over wall
      REAL, DIMENSION(1) :: ZGFLUX_WALL_B      ! flux through the wall
      REAL, DIMENSION(1) :: ZRNSNOW_ROOF       ! net radiation over snow
      REAL, DIMENSION(1) :: ZHSNOW_ROOF        ! sensible heat flux over snow
      REAL, DIMENSION(1) :: ZLESNOW_ROOF       ! latent heat flux over snow
      REAL, DIMENSION(1) :: ZGSNOW_ROOF        ! flux under the snow
      REAL, DIMENSION(1) :: ZMELT_ROOF         ! snow melt
      REAL, DIMENSION(1) :: ZRNSNOW_ROAD       ! net radiation over snow
      REAL, DIMENSION(1) :: ZHSNOW_ROAD        ! sensible heat flux over snow
      REAL, DIMENSION(1) :: ZLESNOW_ROAD       ! latent heat flux over snow
      REAL, DIMENSION(1) :: ZGSNOW_ROAD        ! flux under the snow
      REAL, DIMENSION(1) :: ZMELT_ROAD         ! snow melt
      REAL, DIMENSION(1) :: ZRN_GRND           ! net radiation over ground
      REAL, DIMENSION(1) :: ZH_GRND            ! sensible heat flux over ground
      REAL, DIMENSION(1) :: ZLE_GRND           ! latent heat flux over ground
      REAL, DIMENSION(1) :: ZGFLUX_GRND        ! flux through the ground
      REAL, DIMENSION(1) :: ZRUNOFF_TOWN       ! runoff over the ground
      REAL, DIMENSION(1) :: ZCD                ! town averaged drag coefficient
      REAL, DIMENSION(1) :: ZCDN               ! town averaged neutral drag coefficient
      REAL, DIMENSION(1) :: ZCH_TOWN           ! town averaged heat transfer coefficient
      REAL, DIMENSION(1) :: ZRI_TOWN           ! town averaged Richardson number
      REAL, DIMENSION(1) :: ZRESA_TOWN         ! town aerodynamical resistance
      REAL, DIMENSION(1) :: ZDQS_TOWN          ! heat storage inside town
      REAL, DIMENSION(1) :: ZQF_TOWN           ! total anthropogenic heat
      REAL, DIMENSION(1) :: ZQF_BLD            ! anthropogenic heat flux of domestic heating
      REAL, DIMENSION(1) :: ZFLX_BLD           ! heat flux between inside of the bld and its structure
      REAL, DIMENSION(1) :: ZAC_ROAD           ! road conductance
      ! End unused model outputs

      INTEGER :: JLOOP  ! loop counter

            ! The following are IN arguments to this subroutine but are further modified,
            ! hence need to be copied to a local variable first.
            ZTIME = ZTIME_START
            XDIR_SW = ZDIR_SW
            XSCA_SW = ZSCA_SW
            PFLOOR_HEIGHT = ZFLOOR_HEIGHT
            ZF_WASTE_CAN = XF_WASTE_CAN
            ZCAP_SYS_RAT = XCAP_SYS_RAT
            ZM_SYS_RAT = XM_SYS_RAT

            ! =========================================================================
            !
            ! We try to keep the following module as close as possible to the TEB driver program
            ! to allow to easily update between future versions of TEB.
            ! Any changes to the following code should be requested and made directly into the
            ! TEB driver first and moved here.
            ! -----------------------------------------------------------
            ! Initializations (SYSTEM)
            ! -----------------------------------------------------------
            !
            GSHAD_DAY = .FALSE. ! has shading been necessary this day ?
            GNATVENT_NIGHT =.FALSE. ! has natural ventilation been necessary/possible this night ?
            !
            ! coherence check
            IF ( (.NOT. LGREENROOF) .AND. ZFRAC_GR(1)>0.) THEN
              print*, 'Greenroofs option   is not activated but a non-zero greenroof fraction is given'
              STOP
            END IF
            IF ( (.NOT. LGARDEN) .AND. ZGARDEN(1)>0.) THEN
              print*, 'Garden     option   is not activated but a non-zero garden    fraction is given'
              STOP
            END IF
            IF ( ZBLD(1)+ZGARDEN(1)>=1.) THEN
              print*, 'The sum of garden and building fraction is larger than one, so road fraction is negative. Please check their values.'
              STOP
            END IF
            IF ( (.NOT. LSOLAR_PANEL) .AND. ZFRAC_PANEL(1)>0.) THEN
              print*, 'Solar panels option is not activated but a non-zero solar panels fraction is given'
              STOP
            END IF
            IF ( (.NOT. CBEM=='BEM') .AND. ZGR(1)>0.) THEN
              print*, 'Building Energy Module (BEM) is not activated but a non-zero glazing ratio is given'
              STOP
            END IF
            !
            ! Multi layer option
            LCANOPY= .FALSE.  ! DO NOT CHANGE: whether multi-layer canopy is active.
            !

            CQSAT='OLD' ! saturation is computed relative to water above 0°C, and relative to ice below 0°C
            !
            ! Thresholds
            XCISMIN = 0.5           ! Minimum wind shear
            XVMODMIN = 0.5          ! Minimum wind speed
            LALDTHRES = .FALSE.     ! DO NOT CHANGE: whether to activate Aladin threshold for wind
            XRIMAX = 0.2            ! Maximum Richardson number for exchange coefficients computations
            !
            !
            ZTIME = ZTIME_START
            !
            ! Geometric parameters
            ZROAD         = (1. - ZBLD - ZGARDEN)
            ZCAN_HW_RATIO = 0.5 * ZWALL_O_HOR / (1.-ZBLD)
            ZSVF_ROAD     = (SQRT(ZCAN_HW_RATIO**2+1.) - ZCAN_HW_RATIO)
            ZSVF_GARDEN   = ZSVF_ROAD
            ZSVF_WALL     =  0.5*(ZCAN_HW_RATIO+1.-SQRT(ZCAN_HW_RATIO**2+1.))/ZCAN_HW_RATIO
            ZZ_LOWCAN     = ZBLD_HEIGHT / 2
            ZTOTS_O_HORS = 1. + ZWALL_O_HOR
            ZROOF_FRAC   = ZBLD        / ZTOTS_O_HORS
            ZWALL_FRAC   = ZWALL_O_HOR / ZTOTS_O_HORS
            ZROAD_FRAC   = ZROAD       / ZTOTS_O_HORS
            ZGARDEN_FRAC = ZGARDEN     / ZTOTS_O_HORS
            ZWALL_O_GRND = ZWALL_FRAC  / (ZROAD_FRAC+ZGARDEN_FRAC)
            ZROAD_O_GRND   = ZROAD / (ZROAD + ZGARDEN)
            ZGARDEN_O_GRND = ZGARDEN / (ZROAD + ZGARDEN)
            !
            CALL BEM_MORPHO_STRUCT(ZBLD, ZWALL_O_HOR, ZBLD_HEIGHT, PFLOOR_HEIGHT,                  &
                                  ZGR, ZN_FLOOR, ZWALL_O_BLD, ZGLAZ_O_BLD, ZMASS_O_BLD,     &
                                  ZFLOOR_HW_RATIO,                                          &
                                  ZF_FLOOR_MASS, ZF_FLOOR_WALL, ZF_FLOOR_WIN,               &
                                  ZF_FLOOR_ROOF, ZF_WALL_FLOOR, ZF_WALL_MASS,               &
                                  ZF_WALL_WIN, ZF_WIN_FLOOR, ZF_WIN_MASS, ZF_WIN_WALL,      &
                                  ZF_MASS_FLOOR, ZF_MASS_WALL, ZF_MASS_WIN, ZF_WASTE_CAN,   &
                                  ZF_WIN_WIN)
            !
            ! Window characteristics
            CALL WINDOW_DATA_STRUCT(1, ZSHGC, PU_WIN, ZALB_WIN, ZABS_WIN, ZUGG_WIN, ZTRAN_WIN)
            !
            !
                ! Solar time and position
                  CALL SUNPOS(IYEAR, IMONTH, IDAY, ZTIME, ZLON, ZLAT, XTSUN, XZENITH, XAZIM)
                  CALL SUNPOS(IYEAR, IMONTH, IDAY, ZTIME+XTSTEP_SURF, ZLON, ZLAT, XTSUN, &
                   XZENITH2, XAZIM)
                  !
                  ! Update time
                  ZTIME_BEG = ZTIME               ! Time at beginning of time step
                  ZTIME = ZTIME + XTSTEP_SURF     ! Time at end of time step
                  CALL ADD_FORECAST_TO_DATE_SURF(IYEAR, IMONTH, IDAY, ZTIME)
                  !
                  TPTIME%TIME= ZTIME
                  TPTIME%TDATE%YEAR =IYEAR
                  TPTIME%TDATE%MONTH=IMONTH
                  TPTIME%TDATE%DAY  =IDAY
                  !
                  ! Exner functions
                  ZEXNS = (XPS/XP00)**(XRD/XCPD)
                  ZEXNA = (XPA/XP00)**(XRD/XCPD)
                  !
                  !
                  ! Coherence between solar zenithal angle and radiation
                  ! when solar beam close to horizontal -> reduction of direct radiation to
                  ! the benefit of scattered radiation
                  ! when pi/2 - 0.1 < ZENITH < pi/2 - 0.05 => weight of direct to scattered radiation decreases linearly with zenith
                  ! when pi/2 - 0.05 < ZENITH => all the direct radiation is converted to scattered radiation
                  ! coherence between solar zenithal angle and radiation
                  !
                  ZCOEF(:) = (XPI/2. - XZENITH(:) - 0.05) / 0.05
                  ZCOEF(:) = MAX(MIN(ZCOEF,1.),0.)
                  DO JLOOP=1,SIZE(XDIR_SW,2)
                    XSCA_SW(:,JLOOP) = XSCA_SW(:,JLOOP) + XDIR_SW(:,JLOOP) * (1 - ZCOEF)
                    XDIR_SW(:,JLOOP) = XDIR_SW(:,JLOOP) * ZCOEF(:)
                  ENDDO
                  !
                  ZTDIR_SW = XDIR_SW(1,1)
                  ZTSCA_SW = XSCA_SW(1,1)
                  KSW = 1  ! Only one spectral band here
                  !
                  ! Sky model for diffuse radiation
                  ! Add directionnal contrib from scattered radiation
                  CALL CIRCUMSOLAR_RAD(XDIR_SW(:,1), XSCA_SW(:,1), XZENITH, ZF1_o_B)
                  ZTDIR_SW(:,1) = XDIR_SW(:,1) + XSCA_SW(:,1) * ZF1_o_B
                  ZTSCA_SW(:,1) = XSCA_SW(:,1) * (1. - ZF1_o_B)
                  !
                  ! Wind speed
                  ZVMOD = ZWIND(1)
                  ZWAKE = 1. + (2./XPI-1.) * 2. * (ZCAN_HW_RATIO-0.5)
                  ZWAKE = MAX(MIN(ZWAKE,1.),2./XPI)
                  ZU_CANYON = ZWAKE * EXP(-ZCAN_HW_RATIO/4.) * ZVMOD     &
                          * LOG( (           2.* ZBLD_HEIGHT/3.) / ZZ0)   &
                          / LOG( (ZZREF + 2.* ZBLD_HEIGHT/3.) / ZZ0)
                  ZU_LOWCAN = ZU_CANYON
                  ZT_LOWCAN = ZT_CANYON
                  ZQ_LOWCAN = ZQ_CANYON
                  !
                  ZPEW_A_COEF        = 0.
                  ZPEW_A_COEF_LOWCAN = 0.
                  ZPEW_B_COEF        = ZVMOD
                  ZPEW_B_COEF_LOWCAN = ZU_LOWCAN
                  !
    CALL TEB_GARDEN_STRUCT (LGARDEN, LGREENROOF, LSOLAR_PANEL,                &
                     HZ0H, HIMPLICIT_WIND, HROAD_DIR, HWALL_OPT, TPTIME,      &
                     XTSUN, ZT_CANYON, ZQ_CANYON, ZU_CANYON,                  &
                     ZT_LOWCAN, ZQ_LOWCAN, ZU_LOWCAN, ZZ_LOWCAN,              &
                     ZTI_BLD,                                                 &
                     ZT_ROOF, ZT_ROAD, ZT_WALL_A, ZT_WALL_B,                  &
                     ZWS_ROOF,ZWS_ROAD,                                       &
                     HSNOW_ROOF,                                              &
                     ZWSNOW_ROOF, ZTSNOW_ROOF, ZRSNOW_ROOF, ZASNOW_ROOF,      &
                     ZTSSNOW_ROOF, ZESNOW_ROOF,                               &
                     HSNOW_ROAD,                                              &
                     ZWSNOW_ROAD, ZTSNOW_ROAD, ZRSNOW_ROAD, ZASNOW_ROAD,      &
                     ZTSSNOW_ROAD, ZESNOW_ROAD,                               &
                     ZPEW_A_COEF, ZPEW_B_COEF,                                &
                     ZPEW_A_COEF_LOWCAN, ZPEW_B_COEF_LOWCAN,                  &
                     XPS, XPA, ZEXNS, ZEXNA,                                  &
                     XTA, ZQA_KGKG, XRHOA, XCO2,                              &
                     XLW, ZTDIR_SW, ZTSCA_SW, XSW_BANDS, KSW,                 &
                     XZENITH, XAZIM,                                          &
                     XRAIN, XSNOW,                                            &
                     ZZREF, ZZREF, ZVMOD,                                     &
                     ZH_TRAFFIC, ZLE_TRAFFIC, ZH_INDUSTRY, ZLE_INDUSTRY,      &
                     XTSTEP_SURF,                                             &
                     ZZ0,                                                     &
                     ZBLD,ZGARDEN,ZROAD_DIR,ZROAD,ZFRAC_GR,                   &
                     ZBLD_HEIGHT,ZWALL_O_HOR,ZCAN_HW_RATIO,                   &
                     ZROAD_O_GRND, ZGARDEN_O_GRND, ZWALL_O_GRND,              &
                     ZALB_ROOF, ZEMIS_ROOF,                                   &
                     ZHC_ROOF,ZTC_ROOF,ZD_ROOF,                               &
                     ZALB_ROAD, ZEMIS_ROAD, ZSVF_ROAD,                        &
                     ZHC_ROAD,ZTC_ROAD,ZD_ROAD,                               &
                     ZALB_WALL, ZEMIS_WALL, ZSVF_WALL,                        &
                     ZSVF_GARDEN,                                             &
                     ZHC_WALL,ZTC_WALL,ZD_WALL,                               &
                     ZRN_ROOF, ZH_ROOF, ZLE_ROOF, ZLEW_ROOF, ZGFLUX_ROOF,     &
                     ZRUNOFF_ROOF,                                            &
                     ZRN_ROAD, ZH_ROAD, ZLE_ROAD, ZLEW_ROAD, ZGFLUX_ROAD,     &
                     ZRUNOFF_ROAD,                                            &
                     ZRN_WALL_A, ZH_WALL_A, ZLE_WALL_A, ZGFLUX_WALL_A,        &
                     ZRN_WALL_B, ZH_WALL_B, ZLE_WALL_B, ZGFLUX_WALL_B,        &
                     ZRN_STRLROOF,ZH_STRLROOF,ZLE_STRLROOF, ZGFLUX_STRLROOF,  &
                     ZRUNOFF_STRLROOF,                                        &
                     ZRN_BLT,ZH_BLT,ZLE_BLT, ZGFLUX_BLT,                      &
                     ZRNSNOW_ROOF, ZHSNOW_ROOF, ZLESNOW_ROOF, ZGSNOW_ROOF,    &
                     ZMELT_ROOF,                                              &
                     ZRNSNOW_ROAD, ZHSNOW_ROAD, ZLESNOW_ROAD, ZGSNOW_ROAD,    &
                     ZMELT_ROAD,                                              &
                     ZRN_GRND, ZH_GRND, ZLE_GRND, ZGFLUX_GRND,                &
                     ZRN_TOWN, ZH_TOWN, ZLE_TOWN, ZGFLUX_TOWN, ZEVAP_TOWN,    &
                     ZRUNOFF_TOWN, ZSFCO2,                                    &
                     ZUW_GRND, ZUW_ROOF, ZDUWDU_GRND, ZDUWDU_ROOF,            &
                     ZUSTAR_TOWN, ZCD, ZCDN, ZCH_TOWN, ZRI_TOWN,              &
                     ZTS_TOWN, ZEMIS_TOWN, ZDIR_ALB_TOWN, ZSCA_ALB_TOWN,      &
                     ZRESA_TOWN, ZDQS_TOWN, ZQF_TOWN, ZQF_BLD,                &
                     ZFLX_BLD, ZAC_ROAD, ZAC_GARDEN, ZAC_GREENROOF,           &
                     ZAC_ROAD_WAT, ZAC_GARDEN_WAT, ZAC_GREENROOF_WAT,         &
                     ZABS_SW_ROOF,ZABS_LW_ROOF,                               &
                     ZABS_SW_SNOW_ROOF,ZABS_LW_SNOW_ROOF,                     &
                     ZABS_SW_ROAD,ZABS_LW_ROAD,                               &
                     ZABS_SW_SNOW_ROAD,ZABS_LW_SNOW_ROAD,                     &
                     ZABS_SW_WALL_A,ZABS_LW_WALL_A,                           &
                     ZABS_SW_WALL_B,ZABS_LW_WALL_B,                           &
                     ZABS_SW_PANEL,ZABS_LW_PANEL,                             &
                     ZABS_SW_GARDEN,ZABS_LW_GARDEN,                           &
                     ZABS_SW_GREENROOF,ZABS_LW_GREENROOF,                     &
                     ZG_GREENROOF_ROOF,    &
                     CCOOL_COIL, ZF_WATER_COND, CHEAT_COIL,  &
                     HNATVENT, ZNATVENT, IDAY, ZAUX_MAX, ZT_FLOOR,            &
                     ZT_MASS, ZH_BLD_COOL, ZT_BLD_COOL, ZH_BLD_HEAT,          &
                     ZLE_BLD_COOL, ZLE_BLD_HEAT, ZH_WASTE, ZLE_WASTE,         &
                     ZF_WASTE_CAN, ZHVAC_COOL, ZHVAC_HEAT, ZQIN, ZQIN_FRAD,   &
                     ZQIN_FLAT, ZGR, ZEFF_HEAT, ZINF,                         &
                     ZTCOOL_TARGET, ZTHEAT_TARGET, ZHR_TARGET, ZT_WIN2,       &
                     ZQI_BLD, ZV_VENT, ZCAP_SYS_HEAT, ZCAP_SYS_RAT, ZT_ADP,   &
                     ZM_SYS_RAT, ZCOP_RAT, ZCAP_SYS, ZM_SYS, ZCOP, ZQ_SYS,    &
                     ZT_SYS, ZTR_SW_WIN, ZFAN_POWER, ZHC_FLOOR, ZTC_FLOOR,    &
                     ZD_FLOOR, ZT_WIN1, ZABS_SW_WIN, ZABS_LW_WIN, ZSHGC,      &
                     ZSHGC_SH, ZUGG_WIN, ZALB_WIN, ZABS_WIN, ZEMIT_LW_FAC,    &
                     ZEMIT_LW_GRND, ZT_RAD_IND, ZREF_SW_GRND, ZREF_SW_FAC,    &
                     ZHU_BLD, ZTIME_BEG, LSHADE, ZSHADE, GSHAD_DAY,           &
                     GNATVENT_NIGHT,                                          &
                     CBEM,                                                    &
                     ZN_FLOOR, ZWALL_O_BLD, ZGLAZ_O_BLD, ZMASS_O_BLD,         &
                     ZFLOOR_HW_RATIO, ZF_FLOOR_MASS, ZF_FLOOR_WALL,           &
                     ZF_FLOOR_WIN, ZF_FLOOR_ROOF, ZF_WALL_FLOOR, ZF_WALL_MASS,&
                     ZF_WALL_WIN, ZF_WIN_FLOOR, ZF_WIN_MASS, ZF_WIN_WALL,     &
                     ZF_MASS_FLOOR, ZF_MASS_WALL, ZF_MASS_WIN, LCANOPY,       &
                     ZTRAN_WIN, CCH_BEM, ZROUGH_ROOF, ZROUGH_WALL, ZF_WIN_WIN,&
                     LPAR_RD_IRRIG, ZRD_START_MONTH, ZRD_END_MONTH,           &
                     ZRD_START_HOUR, ZRD_END_HOUR, ZRD_24H_IRRIG, ZIRRIG_ROAD,&
                     ZEMIS_PANEL, ZALB_PANEL, ZEFF_PANEL, ZFRAC_PANEL,        &
                     ZRESIDENTIAL,                                            &
                     ZTHER_PROD_PANEL, ZPHOT_PROD_PANEL, ZPROD_PANEL,         &
                     ZTHER_PROD_BLD  , ZPHOT_PROD_BLD  , ZPROD_BLD  ,         &
                     ZTHER_PRODC_DAY, ZH_PANEL, ZRN_PANEL,                    &
                     ZDT_RES, ZDT_OFF,                                        &
                     ZCUR_TCOOL_TARGET, ZCUR_THEAT_TARGET, ZCUR_QIN           )

      ! Calculate the total (direct + diffuse radiations) town albedo
      IF (ZTDIR_SW(1,1) == 0.) THEN
        ZALB_TOWN = 0.
      ELSE
        ZALB_TOWN = (ZDIR_ALB_TOWN * ZTDIR_SW(1,1) + ZSCA_ALB_TOWN * ZTSCA_SW(1,1)) / (ZTDIR_SW(1,1) + ZTSCA_SW(1,1))
      END IF
      ! TEB does not modify the wind direction
      ZDIR_CANYON = ZDIR(1)
      ! Town specific humidity assumed to be same as canyon humidity
      ZQ_TOWN = ZQ_CANYON

      XU_CANYON = ZU_CANYON

    END SUBROUTINE TEB_DRIVER
END MODULE MODD_WRF_TEB_DRIVER