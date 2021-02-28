!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL licence version 2.1.
! See the following links for details:
! https://cecill.info/licences/Licence_CeCILL_V2.1-en.txt
! https://cecill.info/licences/Licence_CeCILL_V2.1-fr.txt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ======================================================================
PROGRAM DRIVER
! ======================================================================
!
! ......................................................................
!  METHOD
! ......................................................................
!
! Program designed to create, modify, and test TEB routines before being
! actually implemented into another driver or model (e.g. SURFEX or in an
! atmospheric model)
!
! Using just TEB physical routines, new arguments can be created or
! eliminated as local variables in Driver.
!
!
! Notes:
! The current version of Driver does not include the option TEB_CANOPY.
! Outputs have to be declared in Driver.
!
!    AUTHOR
!
!	B. Bueno, Meteo-France
!
!    MODIFICATIONS
!
!      Original    08/12/10
!      Modification   04/13 (V. Masson) adds garden     (with a proxi SVAT)
!                                            greenroofs (with a proxi SVAT)
!                                            road orientation option
!                                            separated walls option
!      Modification   10/13 (V. Masson) adds irrigation and solar panels
! ---------------------------------------------------------------
! Modules
! ---------------------------------------------------------------
!
USE MODD_CSTS,     ONLY : XCPD, XSTEFAN, XPI, XDAY, XKARMAN,   &
                          XLVTT, XLSTT, XLMTT, XRV, XRD, XG, XP00
USE MODD_SURF_ATM, ONLY: XCISMIN, XVMODMIN, LALDTHRES, XRIMAX
USE MODD_SURF_PAR, ONLY: XUNDEF
USE MODD_TYPE_DATE_SURF
USE MODE_THERMOS
USE MODD_REPROD_OPER, ONLY : CQSAT
!
!USE MODI_INIT_SURFCONSPHY
USE MODI_INI_CSTS
USE MODI_ADD_FORECAST_TO_DATE_SURF
USE MODI_SUNPOS
USE MODI_OL_READ_ATM
USE MODI_OL_ALLOC_ATM
USE MODI_OL_TIME_INTERP_ATM
USE MODI_TEB_GARDEN_STRUCT
USE MODI_WINDOW_DATA_STRUCT
USE MODI_BEM_MORPHO_STRUCT
USE MODI_CIRCUMSOLAR_RAD
USE MODI_OPEN_CLOSE_BIN_ASC_FORC
!
USE MODD_FORC_ATM, ONLY: CSV         ,&! name of all scalar variables
                         XDIR_ALB    ,&! direct albedo for each band
                         XSCA_ALB    ,&! diffuse albedo for each band
                         XEMIS       ,&! emissivity
                         XTSRAD      ,&! radiative temperature
                         XTSUN       ,&! solar time (s from midnight)
                         XZS         ,&! orography               (m)
                         XZREF       ,&! height of T,q forcing   (m)
                         XUREF       ,&! height of wind forcing  (m)
                         XTA         ,&! air temperature forcing (K)
                         XQA         ,&! air humidity forcing    (kg/m3)
                         XSV         ,&! scalar variables
                         XU          ,&! zonal wind              (m/s)
                         XV          ,&! meridian wind           (m/s)
                         XDIR_SW     ,&! direct  solar radiation (on horizontal surf.)
                         XSCA_SW     ,&! diffuse solar radiation (on horizontal surf.)
                         XSW_BANDS   ,&! mean wavelength of each shortwave band (m)
                         XZENITH     ,&! zenithal angle  (radian from the vertical)
                         XZENITH2    ,&! zenithal angle  (radian from the vertical)
                         XAZIM       ,&! azimuthal angle (radian from North, clockwise)
                         XLW         ,&! longwave radiation (on horizontal surf.)
                         XPS         ,&! pressure at atmospheric model surface (Pa)
                         XPA         ,&! pressure at forcing level      (Pa)
                         XRHOA       ,&! density at forcing level       (kg/m3)
                         XCO2        ,&! CO2 concentration in the air   (kg/m3)
                         XSNOW       ,&! snow precipitation             (kg/m2/s)
                         XRAIN       ,&! liquid precipitation           (kg/m2/s)
                         XSFTH       ,&! flux of heat                   (W/m2)
                         XSFTQ       ,&! flux of water vapor            (kg/m2/s)
                         XSFU        ,&! zonal momentum flux            (m/s)
                         XSFV        ,&! meridian momentum flux         (m/s)
                         XSFCO2      ,&! flux of CO2                    (kg/m2/s)
                         XSFTS       ,&! flux of scalar var.            (kg/m2/s)
                         XPEW_A_COEF ,&! implicit coefficients
                         XPEW_B_COEF ,&!
                         XPET_A_COEF ,&
                         XPEQ_A_COEF ,&
                         XPET_B_COEF ,&
                         XPEQ_B_COEF

USE PROXI_SVAT_CONSTS, ONLY: VEG_ALB, VEG_EMIS, GARDEN_BR, GREENROOF_BR
!
IMPLICIT NONE
! ---------------------------------------------------------------
! Namelist: Add to allow parameters to be changed dynamically using
!           namelist file.
! ---------------------------------------------------------------
!
namelist /dimensions/   NROOF_LAYER, NROAD_LAYER, NWALL_LAYER, NFLOOR_LAYER

namelist /parameters/   XTSTEP_SURF, IYEAR, IMONTH, IDAY, ZTIME_START, ZLON, ZLAT ,&
                        INB_STEP_ATM, INB_ATM, KSW, ZZREF, CBEM, HROAD_DIR ,&
                        HWALL_OPT, LGARDEN, LGREENROOF, ZFRAC_GR, LSOLAR_PANEL ,&
                        ZFRAC_PANEL, LPAR_RD_IRRIG, CCOOL_COIL ,&
                        CHEAT_COIL, HZ0H, ZZ0, ZBLD, ZGARDEN, ZBLD_HEIGHT ,&
                        ZWALL_O_HOR, ZROAD_DIR, ZALB_ROOF, ZEMIS_ROOF ,&
                        ZHC_ROOF, ZTC_ROOF, ZD_ROOF, ZALB_ROAD, ZEMIS_ROAD ,&
                        ZHC_ROAD, ZTC_ROAD, ZD_ROAD, ZALB_WALL, ZEMIS_WALL, ZHC_WALL ,&
                        ZTC_WALL, ZD_WALL, ZHC_FLOOR, ZTC_FLOOR, ZD_FLOOR, ZH_TRAFFIC ,&
                        ZLE_TRAFFIC, ZH_INDUSTRY, ZLE_INDUSTRY, ZRD_START_MONTH, ZRD_END_MONTH ,&
                        ZRD_START_HOUR, ZRD_END_HOUR, ZRD_24H_IRRIG, ZEMIS_PANEL, ZALB_PANEL ,&
                        ZEFF_PANEL, ZRESIDENTIAL, ZDT_RES, ZDT_OFF, PFLOOR_HEIGHT, ZINF, ZQIN ,&
                        ZQIN_FRAD, ZQIN_FLAT, ZGR, ZSHGC, PU_WIN, ZSHADE, ZSHGC_SH ,&
                        ZNATVENT, ZV_VENT, ZF_WATER_COND, ZF_WASTE_CAN, ZTCOOL_TARGET, ZTHEAT_TARGET ,&
                        ZHR_TARGET, ZEFF_HEAT, ZCAP_SYS_HEAT, ZCAP_SYS_RAT, ZT_ADP, ZM_SYS_RAT ,&
                        ZCOP_RAT, CCH_BEM, ZROUGH_ROOF, ZROUGH_WALL, ZT_ROAD, ZT_ROOF, ZT_WALL ,&
                        ZT_FLOOR, ZT_MASS, ZTI_BLD, ZT_CANYON, ZT_WIN1, ZT_WIN2, ZQ_CANYON, ZQI_BLD, &
                        VEG_ALB, VEG_EMIS, GARDEN_BR, GREENROOF_BR

! ---------------------------------------------------------------
! Declarations of local variables (INPUTS)
! ---------------------------------------------------------------
!
! Simulation - prepared for CAPITOUL atmospheric forcing
REAL               :: XTSTEP_SURF
INTEGER            :: IYEAR      ! Current year (UTC)
INTEGER            :: IMONTH     ! Current month (UTC)
INTEGER            :: IDAY       ! Current day (UTC)
REAL               :: ZTIME_START! Time at start of the run (s)
REAL,DIMENSION(1)  :: ZLON       ! Longitude (deg)
REAL,DIMENSION(1)  :: ZLAT       ! Latitude (deg)
INTEGER           :: INB_STEP_ATM! Forcing time-steps
INTEGER            :: INB_ATM    ! number time the driver calls the TEB
!                                ! routines during a forcing time-step
!                                ! --> it defines the time-step for TEB
INTEGER            :: KSW        ! number of spectral bands in SW forcing
REAL,DIMENSION(1)  :: ZZREF      ! Atm. Forcing height above roof level
!
! GO TO LINE 550 FOR MORE INPUTS
!
! ---------------------------------------------------------------------               ! ||   ||
! declarations of local variables (PARAMETERS)                                        ! ||   ||
! ---------------------------------------------------------------------               ! ||   ||
!                                                                                     ! ||   ||
REAL,DIMENSION(1) :: ZTI_BLD        ! indoor air temperature                          ! ||   ||
REAL,DIMENSION(1) :: ZT_CANYON      ! air canyon temperature                          ! ||   ||
REAL,DIMENSION(1) :: ZQA_KGKG       ! air humidity at forcing level (kg/kg)           ! || G ||
REAL              :: ZTIMEC   = 0.0 ! Current duration since start of the run(s)      ! || O ||
INTEGER           :: JFORC_STEP     ! atmospheric loop index                          ! ||   ||
INTEGER           :: JSURF_STEP     ! Driver loop index                               ! || B ||
INTEGER           :: JLOOP          ! loop counter                                    ! || E ||
INTEGER           :: NWALL_LAYER    ! number of wall layers                           ! || L ||
INTEGER           :: NROOF_LAYER    ! number of roof layers                           ! || O ||
INTEGER           :: NROAD_LAYER    ! number of road layers                           ! || W ||
INTEGER           :: NFLOOR_LAYER   ! number of floor layers                          ! ||   ||
TYPE(DATE_TIME)   :: TPTIME                                                           ! ||   ||
CHARACTER(LEN=3)  :: HIMPLICIT_WIND='NEW' ! Implicitation option for wind fluxes      ! ||   ||
REAL,DIMENSION(:,:), ALLOCATABLE :: ZHC_WALL ! Heat capacity        of wall layers    ! ||   ||
REAL,DIMENSION(:,:), ALLOCATABLE :: ZTC_WALL ! Thermal conductivity of wall layers    ! ||   ||
REAL,DIMENSION(:,:), ALLOCATABLE :: ZD_WALL  ! Thickness            of wall layers    !\\     //
REAL,DIMENSION(:,:), ALLOCATABLE :: ZHC_ROOF ! Heat capacity        of roof layers    ! \\   //
REAL,DIMENSION(:,:), ALLOCATABLE :: ZTC_ROOF ! Thermal conductivity of roof layers    !  \\ //
REAL,DIMENSION(:,:), ALLOCATABLE :: ZD_ROOF  ! Thickness            of roof layers    !   \\/
REAL,DIMENSION(:,:), ALLOCATABLE :: ZHC_ROAD ! Heat capacity        of road layers    !
REAL,DIMENSION(:,:), ALLOCATABLE :: ZTC_ROAD ! Thermal conductivity of road layers    ! ||   ||
REAL,DIMENSION(:,:), ALLOCATABLE :: ZD_ROAD  ! Thickness            of road layers    ! ||   ||
REAL,DIMENSION(:,:), ALLOCATABLE :: ZT_ROAD  ! road layers temperatures               ! || F ||
REAL,DIMENSION(:,:), ALLOCATABLE :: ZT_ROOF  ! roof layers temperatures               ! || O ||
REAL,DIMENSION(:,:), ALLOCATABLE :: ZT_WALL! wall layers temperatures (wall 'A')    ! || R ||
REAL,DIMENSION(:,:), ALLOCATABLE :: ZT_WALL_B! wall layers temperatures (wall 'B')    ! ||   ||
REAL,DIMENSION(1)  :: ZROAD_DIR         ! road direction (° from North, clockwise)    ! || U ||
REAL,DIMENSION(1)  :: ZFRAC_GR          ! fraction of greenroofs on roofs             ! || R ||
REAL,DIMENSION(1)  :: ZFRAC_PANEL       ! fraction of solar panels on roofs           ! || B ||
                                                                                      ! || A ||
!* snow variables                                                                     ! || N ||
CHARACTER(LEN=3)   :: HSNOW_ROOF = "1-L"! Option for roof snow                        ! ||   ||
REAL,DIMENSION(1,1):: ZWSNOW_ROOF= 0.   ! Initial Amount      of roof snow reservoir  ! || P ||
REAL,DIMENSION(1,1):: ZTSNOW_ROOF       ! layer temperature   of roof snow            ! || A ||
REAL,DIMENSION(1,1):: ZRSNOW_ROOF       ! density             of roof snow            ! || R ||
REAL,DIMENSION(1)  :: ZTSSNOW_ROOF      ! surface temperature of roof snow            ! || A ||
CHARACTER(LEN=3)   :: HSNOW_ROAD = "1-L"! Option for road snow                        ! || M ||
REAL,DIMENSION(1,1):: ZWSNOW_ROAD= 0.   ! Initial amount      of road snow reservoir  ! || E ||
REAL,DIMENSION(1,1):: ZTSNOW_ROAD       ! layer temperature   of road snow            ! || T ||
REAL,DIMENSION(1,1):: ZRSNOW_ROAD       ! density             of road snow            ! || E ||
REAL,DIMENSION(1)  :: ZTSSNOW_ROAD      ! surface temperature of road snow            ! || R ||
!* radiative variables                                                                ! || S ||
REAL,DIMENSION(1)  :: ZALB_WALL         ! albedo of walls                             ! ||   ||
REAL,DIMENSION(1)  :: ZEMIS_WALL        ! emissivity of walls                         ! ||   ||
REAL,DIMENSION(1)  :: ZALB_ROOF         ! albedo of roofs                             !\\     //
REAL,DIMENSION(1)  :: ZEMIS_ROOF        ! emissivity of roofs                         ! \\   //
REAL,DIMENSION(1)  :: ZESNOW_ROAD       ! road snow emissivity                        !  \\ //
REAL,DIMENSION(1)  :: ZALB_ROAD         ! albedo of roads                             !   \\/
REAL,DIMENSION(1)  :: ZEMIS_ROAD        ! emissivity of roads                         !
REAL,DIMENSION(1)  :: ZASNOW_ROOF       ! roof snow albedo                            ! ||   ||
REAL,DIMENSION(1)  :: ZASNOW_ROAD       ! road snow albedo                            ! || A ||
REAL,DIMENSION(1)  :: ZABS_SW_ROOF      ! Shortwave radiation absorbed by roofs       ! || N ||
REAL,DIMENSION(1)  :: ZABS_SW_ROAD      ! Shortwave radiation absorbed by roads       ! || D ||
REAL,DIMENSION(1)  :: ZABS_SW_WALL_A    ! Shortwave radiation absorbed by wall A      ! ||   ||
REAL,DIMENSION(1)  :: ZABS_SW_WALL_B    ! Shortwave radiation absorbed by wall B      ! || T ||
REAL,DIMENSION(1)  :: ZABS_SW_GARDEN    ! Shortwave radiation absorbed by gardens     ! || E ||
REAL,DIMENSION(1)  :: ZABS_SW_GREENROOF ! Shortwave radiation absorbed by greenroofs  ! || B ||
REAL,DIMENSION(1)  :: ZABS_LW_ROOF      ! Longwave  radiation absorbed by roofs       ! ||   ||
REAL,DIMENSION(1)  :: ZABS_LW_ROAD      ! Longwave  radiation absorbed by roads       ! || O ||
REAL,DIMENSION(1)  :: ZABS_LW_WALL_A    ! Longwave  radiation absorbed by wall A      ! || P ||
REAL,DIMENSION(1)  :: ZABS_LW_WALL_B    ! Longwave  radiation absorbed by wall B      ! || T ||
REAL,DIMENSION(1)  :: ZABS_LW_GARDEN    ! Longwave  radiation absorbed by gardens     ! || I ||
REAL,DIMENSION(1)  :: ZABS_LW_GREENROOF ! Longwave  radiation absorbed by greenroofs  ! || O ||
REAL,DIMENSION(1)  :: ZABS_SW_SNOW_ROOF ! Shortwave radiation absorbed by roof snow   ! || N ||
REAL,DIMENSION(1)  :: ZABS_SW_SNOW_ROAD ! Shortwave radiation absorbed by road snow   ! || S ||
REAL,DIMENSION(1)  :: ZABS_SW_PANEL  ! Shortwave radiation absorbed by solar panels   ! ||   ||
REAL,DIMENSION(1)  :: ZABS_LW_PANEL  ! Longwave  radiation absorbed by solar panels   ! ||   ||
REAL,DIMENSION(1)  :: ZDIR_ALB_TOWN  ! town direct albedo                             ! ||   ||
REAL,DIMENSION(1)  :: ZSCA_ALB_TOWN  ! town scaterred albedo                          ! ||   ||
REAL,DIMENSION(1)  :: ZCAN_HW_RATIO  ! canyon    h/W                                  ! ||   ||
REAL,DIMENSION(1)  :: ZSVF_ROAD      ! road sky view factor                           ! ||   ||
REAL,DIMENSION(1)  :: ZSVF_GARDEN    ! garden sky view factor                         ! ||   ||
REAL,DIMENSION(1)  :: ZSVF_WALL      ! wall sky view factor                           ! ||   ||
REAL,DIMENSION(1)  :: ZWAKE          ! reduction of average wind speed                ! ||   ||
REAL,DIMENSION(1)  :: ZGSNOW_ROOF  = 0.0 ! roof snow conduction                       !\\     //
REAL,DIMENSION(1)  :: ZHSNOW_ROOF                                                     ! \\   //
! Anthropogenic heat fluxes                                                           !  \\ //
REAL,DIMENSION(1)  :: ZH_TRAFFIC       ! heat fluxes due to traffic                   !   \\/
REAL,DIMENSION(1)  :: ZLE_TRAFFIC      ! heat fluxes due to traffic                   !
REAL,DIMENSION(1)  :: XH_TRAFFIC       ! heat fluxes due to traffic                   !   \\/
REAL,DIMENSION(1)  :: XLE_TRAFFIC      ! heat fluxes due to traffic
REAL,DIMENSION(1)  :: ZH_INDUSTRY      ! heat fluxes due to factories                 ! ||   ||
REAL,DIMENSION(1)  :: ZLE_INDUSTRY     ! heat fluxes due to factories                 ! ||   ||
! Urban options                                                                       ! ||   ||
LOGICAL            :: LPAR_RD_IRRIG    ! Flag for road watering                       ! || S ||
LOGICAL            :: LSOLAR_PANEL     ! Flag to use a solar panels on roofs          ! || P ||
LOGICAL            :: LGARDEN          ! Flag to use a vegetation scheme in gardens   ! || E ||
LOGICAL            :: LGREENROOF       ! Flag to use a green roofs scheme             ! || C ||
CHARACTER(LEN=6)   :: HZ0H             ! TEB option for z0h roof & road               ! || I ||
CHARACTER(LEN=4)   :: HROAD_DIR        ! road direction option :                      ! || F ||
!                                      ! 'UNIF' : uniform roads                       ! || I ||
!                                      ! 'ORIE' : specified road orientation          ! || C ||
CHARACTER(LEN=4)   :: HWALL_OPT        ! Wall option                                  ! || A ||
!                                      ! 'UNIF' : uniform walls                       ! || T ||
!                                      ! 'TWO ' : 2 opposite  walls                   ! || I ||
! Canyon geometry                                                                     ! || O ||
REAL,DIMENSION(1)  :: ZZ0              ! Roughness length (m)                         ! || N ||
REAL,DIMENSION(1)  :: ZBLD             ! Horizontal building area density             ! || S ||
REAL,DIMENSION(1)  :: ZGARDEN          ! fraction of GARDEN areas                     ! ||   ||
REAL,DIMENSION(1)  :: ZBLD_HEIGHT      ! Canyon height (m)                            ! ||   ||
REAL,DIMENSION(1)  :: ZWALL_O_HOR      ! Vertical to horizonal surf ratio             ! ||   ||
! Urban variables                                                                     ! ||   ||
REAL,DIMENSION(1)  :: ZU_CANYON = 1.0  ! hor. wind in canyon                          ! ||   ||
REAL,DIMENSION(1)  :: ZQ_CANYON        ! canyon air specific humidity     (kg/kg)     ! ||   ||
REAL,DIMENSION(1)  :: ZDIR_CANYON      ! canyon wind direction                        ! ||   ||
REAL,DIMENSION(1)  :: ZTS_ROOF         ! surface temperature                          ! ||   ||
REAL,DIMENSION(1)  :: ZTS_ROAD         ! surface temperature                          !\\     //
REAL,DIMENSION(1)  :: ZTS_WALL_A       ! surface temperature                          ! \\   //
REAL,DIMENSION(1)  :: ZTS_WALL_B       ! surface temperature                          !  \\ //
!                                                                                     !   \\/
REAL,DIMENSION(1)  :: ZZ_LOWCAN        ! altitude of air layer above road             !
REAL,DIMENSION(1)  :: ZT_LOWCAN        ! temperature of air above road                ! ||   ||
REAL,DIMENSION(1)  :: ZU_LOWCAN        ! wind above road                              ! ||   ||
REAL,DIMENSION(1)  :: ZQ_LOWCAN        ! humidity above road                          ! ||   ||
REAL,DIMENSION(1)  :: ZEXNS            ! surface exner function                       ! ||   ||
REAL,DIMENSION(1)  :: ZEXNA            ! exner function                               ! ||   ||
REAL,DIMENSION(1)  :: ZROAD            ! fraction of roads                            ! ||   ||
REAL,DIMENSION(1)  :: ZWS_ROOF     = 0.! roof water content (kg/m2)                   ! ||   ||
REAL,DIMENSION(1)  :: ZWS_ROAD     = 0.! road water content (kg/m2)                   ! ||   ||
REAL,DIMENSION(1)  :: ZPEW_A_COEF      ! implicit coefficients                        ! ||   ||
REAL,DIMENSION(1)  :: ZPEW_B_COEF      ! for wind coupling                            ! ||   ||
REAL,DIMENSION(1)  :: ZPEW_A_COEF_LOWCAN!                                             ! ||   ||
REAL,DIMENSION(1)  :: ZPEW_B_COEF_LOWCAN!                                             ! ||   ||
REAL,DIMENSION(1)  :: ZWALL_O_GRND     ! Wall to ground surface ratio                 ! ||   ||
REAL,DIMENSION(1)  :: ZROAD_O_GRND     ! Road to ground surface ratio                 ! ||   ||
REAL,DIMENSION(1)  :: ZGARDEN_O_GRND   ! Garden to ground surface ratio               ! ||   ||
REAL,DIMENSION(1)  :: ZVMOD            ! Module of wind speed at the top of the canyon! ||   ||
REAL,DIMENSION(1)  :: ZCD              ! drag coefficient                             ! ||   ||
REAL,DIMENSION(1)  :: ZCDN             ! neutral drag coefficient                     ! ||   ||
REAL,DIMENSION(1)  :: ZAC_ROAD         ! road aerodynamical conductance               ! ||   ||
REAL,DIMENSION(1)  :: ZAC_ROAD_WAT     ! road aerodynamical conductance (for water)   ! ||   ||
REAL,DIMENSION(1)  :: ZAC_GARDEN       ! garden aerodynamical conductance             !\\     //
REAL,DIMENSION(1)  :: ZAC_GARDEN_WAT   ! garden aerodynamical conductance for vapor   ! \\   //
REAL,DIMENSION(1)  :: ZAC_GREENROOF    ! green roofs aerodynamical conductance        !  \\ //
REAL,DIMENSION(1)  :: ZAC_GREENROOF_WAT! green roofs aerodynamical conductance        !   \\/
!                                      ! for vapor                                    !
REAL,DIMENSION(1)  :: ZUW_ROOF         ! Momentum flux for roofs                      ! ||   ||
REAL,DIMENSION(1)  :: ZDUWDU_GRND      ! d(u'w')/du for ground                        ! ||   ||
REAL,DIMENSION(1)  :: ZDUWDU_ROOF      ! d(u'w')/du for roof                          ! ||   ||
REAL,DIMENSION(1)  :: ZUSTAR_TOWN      ! Fraction velocity for town                   ! ||   ||
REAL,DIMENSION(1)  :: ZRESA_TOWN       ! Aerodynamical resistance                     ! ||   ||
REAL,DIMENSION(1)  :: ZCH_TOWN         ! Heat exchange coefficient                    ! ||   ||
REAL,DIMENSION(1)  :: ZRI_TOWN         ! Richardson number                            ! ||   ||
REAL,DIMENSION(1)  :: ZRUNOFF_TOWN     ! runoff for town                              ! ||   ||
REAL,DIMENSION(1)  :: ZRUNOFF_ROOF     ! runoff for roof                              ! ||   ||
REAL,DIMENSION(1)  :: ZRUNOFF_ROAD     ! runoff for road                              ! ||   ||
! Road-wall-budget routines                                                           ! ||   ||
REAL,DIMENSION(1)  :: ZGSNOW_ROAD = 0.0! road snow conduction                         ! ||   ||
! Urban fluxes variables                                                              ! ||   ||
REAL,DIMENSION(1)  :: ZROOF_FRAC       ! roof, wall,                                  ! ||   ||
REAL,DIMENSION(1)  :: ZWALL_FRAC       ! road, and green area                         ! ||   ||
REAL,DIMENSION(1)  :: ZROAD_FRAC       ! fractions                                    ! ||   ||
REAL,DIMENSION(1)  :: ZGARDEN_FRAC     ! fractions                                    ! ||   ||
REAL,DIMENSION(1)  :: ZTOTS_O_HORS     ! total canyon+roof surf over horizontal surf  ! ||   ||
REAL,DIMENSION(1)  :: ZUW_GRND         ! friction flux over ground                    ! ||   ||
REAL,DIMENSION(1)  :: ZRN_ROOF         ! net radiation over roof                      ! ||   ||
REAL,DIMENSION(1)  :: ZH_ROOF          ! sensible heat flux over roof                 ! ||   ||
REAL,DIMENSION(1)  :: ZLE_ROOF         ! latent heat flux over roof                   !\\     //
REAL,DIMENSION(1)  :: ZGFLUX_ROOF      ! flux through the roof                        ! \\   //
REAL,DIMENSION(1)  :: ZRN_ROAD         ! net radiation over road                      !  \\ //
REAL,DIMENSION(1)  :: ZH_ROAD          ! sensible heat flux over road                 !   \\/
REAL,DIMENSION(1)  :: ZLE_ROAD         ! latent heat flux over road                   !
REAL,DIMENSION(1)  :: ZGFLUX_ROAD      ! flux through the road                        ! ||   ||
REAL,DIMENSION(1)  :: ZRN_WALL_A       ! net radiation over wall                      ! ||   ||
REAL,DIMENSION(1)  :: ZRN_WALL_B       ! net radiation over wall                      ! ||   ||
REAL,DIMENSION(1)  :: ZH_WALL_A = 0.0  ! sensible heat flux over wall                 ! ||   ||
REAL,DIMENSION(1)  :: ZH_WALL_B = 0.0  ! sensible heat flux over wall                 ! ||   ||
REAL,DIMENSION(1)  :: ZLE_WALL_A       ! latent heat flux over wall                   ! ||   ||
REAL,DIMENSION(1)  :: ZLE_WALL_B       ! latent heat flux over wall                   ! ||   ||
REAL,DIMENSION(1)  :: ZGFLUX_WALL_A    ! flux through the wall                        ! ||   ||
REAL,DIMENSION(1)  :: ZGFLUX_WALL_B    ! flux through the wall                        ! ||   ||
                                                                                      ! ||   ||
REAL,DIMENSION(1)  :: ZG_GREENROOF_ROOF! heat flux between base of greenroof          ! ||   ||
!                                      ! and structural roof                          ! ||   ||
                                                                                      ! ||   ||
REAL,DIMENSION(1)  :: ZRNSNOW_ROOF = 0.0 ! net radiation over snow                    ! ||   ||
REAL,DIMENSION(1)  :: ZRNSNOW_ROAD = 0.0 ! net radiation over snow                    ! ||   ||
REAL,DIMENSION(1)  :: ZHSNOW_ROAD  = 0.0 ! sensible heat flux over snow               ! ||   ||
REAL,DIMENSION(1)  :: ZLESNOW_ROAD = 0.0 ! latent heat flux over snow                 ! ||   ||
REAL,DIMENSION(1)  :: ZRN_BLT            ! net radiation over built covers            ! ||   ||
REAL,DIMENSION(1)  :: ZH_BLT             ! sensible heat flux over built covers       ! ||   ||
REAL,DIMENSION(1)  :: ZLE_BLT            ! latent heat flux over built covers         !\\     //
REAL,DIMENSION(1)  :: ZGFLUX_BLT         ! flux through the built covers              ! \\   //
REAL,DIMENSION(1)  :: ZFLX_BLD           ! heat flx from bld to its structure         !  \\ //
REAL,DIMENSION(1)  :: ZDQS_TOWN          ! storage inside town materials              !   \\/
REAL,DIMENSION(1)  :: ZQF_TOWN           ! total anthropogenic heat                   !
REAL,DIMENSION(1)  :: ZMELT_ROOF = 0.0   ! snow melting on roof                       ! ||   ||
REAL,DIMENSION(1)  :: ZMELT_ROAD  = 0.0  ! snow melting on road                       ! ||   ||
REAL,DIMENSION(1)  :: ZQF_BLD                                                         ! ||   ||
REAL,DIMENSION(1)  :: ZSFCO2             ! Surface Flux of CO2                        ! ||   ||
!                                                                                     ! ||   ||
! new arguments created after BEM                                                     ! ||   ||
!                                                                                     ! ||   ||
CHARACTER(LEN=3)    :: CBEM              ! Building Energy model 'DEF' or 'BEM'       ! ||   ||
CHARACTER(LEN=12)   :: CCOOL_COIL        ! option for cooling device type             ! ||   ||
CHARACTER(LEN=6)    :: CHEAT_COIL        ! option for heating device type             ! ||   ||
REAL, DIMENSION(1)  :: ZF_WATER_COND     ! fraction of evaporation for the condensers ! ||   ||
CHARACTER(LEN=4), DIMENSION(1) :: HNATVENT                                            ! ||   ||
REAL, DIMENSION(1)  :: ZNATVENT          ! flag to describe surventilation system for ! ||   ||
!                                        ! i/o 0 for NONE, 1 for MANU and 2 for AUTO  ! ||   ||
REAL, DIMENSION(1)  :: ZAUX_MAX = 5. ! Auxiliar variable for autosize calcs (not used)! ||   ||
REAL, DIMENSION(:,:), ALLOCATABLE :: ZT_FLOOR ! Floor layers temperatures [K]         ! ||   ||
REAL, DIMENSION(:,:), ALLOCATABLE :: ZT_MASS  ! Internal mass layers temperatures [K] ! ||   ||
REAL, DIMENSION(1) :: PFLOOR_HEIGHT      ! Floor height (m)                           ! ||   ||
REAL, DIMENSION(1) :: ZH_BLD_COOL        ! Sensible cooling energy demand             ! ||   ||
                                         ! of the building [W m-2(bld)]               ! ||   ||
REAL, DIMENSION(1) :: ZT_BLD_COOL        ! Total cooling energy demand                ! ||   ||
                                         ! of the building [W m-2(bld)]               ! ||   ||
REAL, DIMENSION(1) :: ZH_BLD_HEAT        ! Heating energy demand                      ! ||   ||
                                         ! of the building [W m-2(bld)]               ! ||   ||
REAL, DIMENSION(1) :: ZLE_BLD_COOL       ! Latent cooling energy demand               !\\     //
                                         ! of the building [W m-2(bld)]               ! \\   //
REAL, DIMENSION(1) :: ZLE_BLD_HEAT       ! Latent heating energy demand               !  \\ //
                                         ! of the building [W m-2(bld)]               !   \\/
REAL, DIMENSION(1) :: ZH_WASTE           ! Sensible waste heat from HVAC system       !
                                         ! [W m-2(tot)]                               ! ||   ||
REAL, DIMENSION(1) :: ZLE_WASTE          ! Latent waste heat from HVAC system         ! ||   ||
                                         ! [W m-2(tot)]                               ! ||   ||
REAL, DIMENSION(1) :: ZF_WASTE_CAN       ! fraction of waste heat released            ! ||   ||
!                                        ! into the canyon                            ! ||   ||
REAL, DIMENSION(1) :: ZHVAC_COOL         ! Energy consumption of the cooling system   ! ||   ||
                                         ! [W m-2(bld)]                               ! ||   ||
REAL, DIMENSION(1) :: ZHVAC_HEAT         ! Energy consumption of the heating system   ! ||   ||
                                         ! [W m-2(bld)]                               ! ||   ||
REAL, DIMENSION(1) :: ZHVAC_COOL_TOT     ! Energy consumption of the cooling system   ! ||   ||
                                         ! [W m-2(tot)]                               ! ||   ||
REAL, DIMENSION(1) :: ZHVAC_HEAT_TOT     ! Energy consumption of the heating system   ! ||   ||
                                         ! [W m-2(tot)]                               ! ||   ||
REAL, DIMENSION(1) :: ZQIN               ! Internal heat gains [W m-2(floor)]         ! ||   ||
REAL, DIMENSION(1) :: ZQIN_FRAD          ! Radiant fraction of internal heat gains    ! ||   ||
REAL, DIMENSION(1) :: ZQIN_FLAT          ! Latent franction of internal heat gains    ! ||   ||
REAL, DIMENSION(1) :: ZGR                ! Glazing ratio                              ! ||   ||
REAL, DIMENSION(1) :: ZEFF_HEAT          ! Efficiency of the heating system           ! ||   ||
REAL, DIMENSION(1) :: ZINF               ! Infiltration flow rate [AC/H]              ! ||   ||
REAL, DIMENSION(1) :: ZTCOOL_TARGET      ! Cooling setpoint of HVAC system [K]        ! ||   ||
REAL, DIMENSION(1) :: ZTHEAT_TARGET      ! Heating setpoint of HVAC system [K]        ! ||   ||
REAL, DIMENSION(1) :: ZCUR_QIN           ! Internal heat gains [W m-2(floor)]         ! ||   ||
REAL, DIMENSION(1) :: ZCUR_TCOOL_TARGET  ! Cooling setpoint of HVAC system [K]        ! ||   ||
REAL, DIMENSION(1) :: ZCUR_THEAT_TARGET  ! Heating setpoint of HVAC system [K]        ! ||   ||
REAL, DIMENSION(1) :: ZHR_TARGET         ! Relative humidity setpoint                 ! ||   ||
REAL, DIMENSION(1) :: ZT_WIN2            ! Indoor window temperature [K]              !\\     //
REAL, DIMENSION(1) :: ZQI_BLD            ! Indoor air specific humidity [kg kg-1]     ! \\   //
REAL, DIMENSION(1) :: ZV_VENT            ! Ventilation flow rate [AC/H]               !  \\ //
REAL, DIMENSION(1) :: ZCAP_SYS_HEAT      ! Capacity of the heating system             !   \\/
                                         ! [W m-2(bld)]                               !
REAL, DIMENSION(1) :: ZCAP_SYS_RAT       ! Rated capacity of the cooling system       ! ||   ||
                                         ! [W m-2(bld)]                               ! ||   ||
REAL, DIMENSION(1) :: ZT_ADP             ! Apparatus dewpoint temperature of the      ! ||   ||
                                         ! cooling coil [K]                           ! ||   ||
REAL, DIMENSION(1) :: ZM_SYS_RAT         ! Rated HVAC mass flow rate                  ! ||   ||
                                         ! [kg s-1 m-2(bld)]                          ! ||   ||
REAL, DIMENSION(1) :: ZCOP_RAT           ! Rated COP of the cooling system            ! ||   ||
REAL, DIMENSION(1) :: ZCAP_SYS           ! Actual capacity of the cooling system      ! ||   ||
                                         ! [W m-2(bld)]                               ! ||   ||
REAL, DIMENSION(1) :: ZM_SYS             ! Actual HVAC mass flow rate                 ! ||   ||
                                         ! [kg s-1 m-2(bld)]                          ! ||   ||
REAL, DIMENSION(1) :: ZCOP               ! COP of the cooling system                  ! ||   ||
REAL, DIMENSION(1) :: ZQ_SYS             ! Supply air specific humidity [kg kg-1]     ! ||   ||
REAL, DIMENSION(1) :: ZT_SYS             ! Supply air temperature [K]                 ! ||   ||
REAL, DIMENSION(1) :: ZTR_SW_WIN         ! Solar radiation transmitted throught       ! ||   ||
                                         ! windows [W m-2(bld)]                       ! ||   ||
REAL, DIMENSION(1) :: ZFAN_POWER         ! HVAC fan power                             ! ||   ||
REAL, DIMENSION(:,:), ALLOCATABLE :: ZHC_FLOOR ! heat capacity for floor layers       ! ||   ||
REAL, DIMENSION(:,:), ALLOCATABLE :: ZTC_FLOOR ! thermal conductivity for             ! ||   ||
                                               ! floor layers                         ! ||   ||
REAL, DIMENSION(:,:), ALLOCATABLE :: ZD_FLOOR  ! depth of floor layers                ! ||   ||
                                                                                      ! ||   ||
REAL, DIMENSION(1) :: ZT_WIN1            ! outdoor window temperature [K]             ! ||   ||
REAL, DIMENSION(1) :: ZABS_SW_WIN        ! window absorbed shortwave radiation [W m-2]!\\     //
REAL, DIMENSION(1) :: ZABS_LW_WIN        ! absorbed infrared rad. [W m-2]             ! \\   //
REAL, DIMENSION(1) :: ZSHGC              ! window solar transmittance                 !  \\ //
REAL, DIMENSION(1) :: ZSHGC_SH           ! window + shading solar heat gain coef.     !   \\/
REAL, DIMENSION(1) :: ZUGG_WIN           ! window glass-to-glass U-factro [W m-2 K-1] !
REAL, DIMENSION(1) :: ZALB_WIN           ! window albedo                              ! ||   ||
REAL, DIMENSION(1) :: PU_WIN             ! window U-factor [K m W-2]                  ! ||   ||
REAL, DIMENSION(1) :: ZABS_WIN           ! window absortance                          ! ||   ||
REAL, DIMENSION(1) :: ZTRAN_WIN          ! window transmittance                       ! ||   ||
!                                                                                     ! ||   ||
! new argument for the UTCI calculation                                               ! ||   ||
REAL, DIMENSION(1)  :: ZEMIT_LW_GRND     ! LW flux emitted by the ground (W/m² ground)! ||   ||
REAL, DIMENSION(1)  :: ZEMIT_LW_FAC      ! LW flux emitted by the facade (W/m² ground)! ||   ||
REAL, DIMENSION(1)  :: ZT_RAD_IND        ! Indoor mean radiant temperature [K]        ! ||   ||
REAL, DIMENSION(1)  :: ZREF_SW_GRND      ! total solar rad reflected from ground      ! ||   ||
REAL, DIMENSION(1)  :: ZREF_SW_FAC       ! total solar rad reflected from facade      ! ||   ||
REAL, DIMENSION(1)  :: ZHU_BLD           ! Indoor relative humidity 0 < (-) < 1       ! ||   ||
!                                                                                     ! ||   ||
! solar panels                                                                        ! ||   ||
REAL, DIMENSION(1)  :: ZEMIS_PANEL       ! Emissivity of solar panel [-]              ! ||   ||
REAL, DIMENSION(1)  :: ZALB_PANEL        ! albedo of solar panel  [-]                 ! ||   ||
REAL, DIMENSION(1)  :: ZEFF_PANEL        ! Efficiency of solar panel [-]              ! ||   ||
REAL, DIMENSION(1)  :: ZTHER_PROD_PANEL  ! Thermal      Energy production             ! ||   ||
!                                        !     of solar panel on roofs (W/m2 panel)   ! ||   ||
REAL, DIMENSION(1)  :: ZPHOT_PROD_PANEL  ! Photovoltaic Energy production             ! ||   ||
!                                        !     of solar panel on roofs (W/m2 panel)   ! ||   ||
REAL, DIMENSION(1)  :: ZPROD_PANEL       ! Averaged     Energy production             ! ||   ||
!                                        !     of solar panel on roofs (W/m2 panel)   ! ||   ||
REAL, DIMENSION(1)  :: ZTHER_PROD_BLD    ! Thermal      Energy production             ! ||   ||
!                                        !     of solar panel on roofs (W/m2 bld  )   ! ||   ||
REAL, DIMENSION(1)  :: ZPHOT_PROD_BLD    ! Photovoltaic Energy production             ! ||   ||
!                                        !     of solar panel on roofs (W/m2 bld  )   !\\     //
REAL, DIMENSION(1)  :: ZPROD_BLD         ! Averaged     Energy production             ! \\   //
!                                        !     of solar panel on roofs (W/m2 bld  )   !  \\ //
REAL, DIMENSION(1)  :: ZTHER_PRODC_DAY=0.! Present day integrated thermal production  !   \\/
!                                        ! of energy (J/m2 panel). zero value at start!
REAL, DIMENSION(1)  :: ZH_PANEL          ! Sensible heat flux from solar panels       ! ||   ||
!                                        !     (W/m2 panel)                           ! ||   ||
REAL, DIMENSION(1)  :: ZRN_PANEL         ! Net radiation        of solar panel        ! ||   ||
!                                        !     (W/m2 panel)                           ! ||   ||
!                                                                                     ! ||   ||
REAL, DIMENSION(1)  :: ZRESIDENTIAL      ! Fraction of residential use in buildings(-)! ||   ||
REAL                :: ZDT_RES           ! target temperature change when unoccupied  ! ||   ||
!                                        ! (K) (residential buildings)                ! ||   ||
REAL                :: ZDT_OFF           ! target temperature change when unoccupied  ! ||   ||
!                                        ! (K) (office buildings)                     ! ||   ||
!                                                                                     ! ||   ||
! road watering                                                                       ! ||   ||
REAL, DIMENSION(1)  :: ZRD_START_MONTH   ! start month for watering of roads(included)! ||   ||
REAL, DIMENSION(1)  :: ZRD_END_MONTH     ! end   month for watering of roads(included)! ||   ||
REAL, DIMENSION(1)  :: ZRD_START_HOUR    ! start hour  for watering of roads(included)! ||   ||
REAL, DIMENSION(1)  :: ZRD_END_HOUR      ! end   hour  for watering of roads(excluded)! ||   ||
REAL, DIMENSION(1)  :: ZRD_24H_IRRIG     ! 24h quantity of water used                 ! ||   ||
!                                        ! for road watering (liter/m2)               ! ||   ||
REAL, DIMENSION(1)  :: ZIRRIG_ROAD       ! road irrigation during current time-step   ! ||   ||
!                                                                                     ! ||   ||
! new arguments for shading, schedule or natural ventilation                          ! ||   ||
LOGICAL,DIMENSION(1) :: LSHADE            ! Flag to use shading devices               ! ||   ||
REAL,   DIMENSION(1) :: ZSHADE            ! flag to activate shading devices          ! ||   ||
!                                         ! -> REAL for i/o 0. or 1                   ! ||   ||
LOGICAL,DIMENSION(1) :: GSHAD_DAY         ! has shading been necessary this day ?     ! ||   ||
LOGICAL,DIMENSION(1) :: GNATVENT_NIGHT    ! has natural ventilation been              ! ||   ||
!                                         ! necessary/possible this night ?           ! ||   ||
!                                                                                     ! ||   ||
REAL, DIMENSION(1)  :: ZN_FLOOR          ! Number of floors                           ! ||   ||
REAL, DIMENSION(1)  :: ZWALL_O_BLD       ! Wall area [m2_wall/m2_bld]                 ! ||   ||
REAL, DIMENSION(1)  :: ZGLAZ_O_BLD       ! Window area [m2_win/m2_bld]                ! ||   ||
REAL, DIMENSION(1)  :: ZMASS_O_BLD       ! Mass area [m2_mass/m2_bld]                 ! ||   ||
REAL, DIMENSION(1)  :: ZFLOOR_HW_RATIO   ! H/W ratio of 1 floor level                 ! ||   ||
REAL, DIMENSION(1)  :: ZF_FLOOR_MASS     ! View factor floor-mass                     ! ||   ||
REAL, DIMENSION(1)  :: ZF_FLOOR_WALL     ! View factor floor-wall                     !\\     //
REAL, DIMENSION(1)  :: ZF_FLOOR_WIN      ! View factor floor-window                   ! \\   //
REAL, DIMENSION(1)  :: ZF_FLOOR_ROOF     ! View factor floor-roof                     !  \\ //
REAL, DIMENSION(1)  :: ZF_WALL_FLOOR     ! View factor wall-floor                     !   \\/
REAL, DIMENSION(1)  :: ZF_WALL_MASS      ! View factor wall-mass                      !
REAL, DIMENSION(1)  :: ZF_WALL_WIN       ! View factor wall-win                       ! ||   ||
REAL, DIMENSION(1)  :: ZF_WIN_FLOOR      ! View factor win-floor                      ! ||   ||
REAL, DIMENSION(1)  :: ZF_WIN_MASS       ! View factor win-mass                       ! ||   ||
REAL, DIMENSION(1)  :: ZF_WIN_WALL       ! View factor win-wall                       ! ||   ||
REAL, DIMENSION(1)  :: ZF_MASS_FLOOR     ! View factor mass-floor                     ! ||   ||
REAL, DIMENSION(1)  :: ZF_MASS_WALL      ! View factor mass-wall                      ! ||   ||
REAL, DIMENSION(1)  :: ZF_MASS_WIN       ! View factor mass-window                    ! ||   ||
LOGICAL             :: LCANOPY           ! is canopy active ?                         ! ||   ||
CHARACTER(LEN=5)    :: CCH_BEM           ! TEB option for building outside conv. coef ! ||   ||
REAL, DIMENSION(1)  :: ZROUGH_ROOF       ! roof roughness coef.                       ! ||   ||
REAL, DIMENSION(1)  :: ZROUGH_WALL       ! wall roughness coef.                       ! ||   ||
REAL, DIMENSION(1)  :: ZF_WIN_WIN        ! indoor win to win view factor              ! ||   ||
                                                                                      ! ||   ||
! Average urban fluxes variables                                                      ! ||   ||
REAL,DIMENSION(1)  :: ZTS_TOWN          ! town surface temperature                    ! ||   ||
REAL,DIMENSION(1)  :: ZEMIS_TOWN        ! town equivalent emissivity                  ! ||   ||
REAL,DIMENSION(1)  :: ZALB_TOWN         ! town eqivalent albedo                       ! ||   ||
REAL,DIMENSION(1)  :: ZQ_TOWN           ! town eqivalent specific humidity            ! ||   ||
REAL,DIMENSION(1)  :: ZESNOW_ROOF  = 1.0! snow roof emissivity                        ! ||   ||
REAL,DIMENSION(1)  :: ZABS_LW_SNOW_ROOF ! abs. LW rad. by snow                        ! ||   ||
REAL,DIMENSION(1)  :: ZABS_LW_SNOW_ROAD ! abs. LW rad. by snow                        ! ||   ||
REAL,DIMENSION(1)  :: ZRN_STRLROOF      ! net radiation over structural roof          !
REAL,DIMENSION(1)  :: ZH_STRLROOF       ! sensible heat flux over structural roof     ! ||   ||
REAL,DIMENSION(1)  :: ZLE_STRLROOF      ! latent heat flux over structural roof       ! ||   ||
REAL,DIMENSION(1)  :: ZGFLUX_STRLROOF   ! flux through the structural roof            ! ||   ||
REAL,DIMENSION(1)  :: ZRUNOFF_STRLROOF  ! water runoff on the structural roof         ! ||   ||
REAL,DIMENSION(1)  :: ZLEW_ROOF         ! latent heat flux of snowfree roof           ! ||   ||
REAL,DIMENSION(1)  :: ZLESNOW_ROOF      ! latent heat flux over snow                  ! ||   ||
REAL,DIMENSION(1)  :: ZLEW_ROAD         ! latent heat flux of snowfree road           ! ||   ||
REAL,DIMENSION(1)  :: ZRN_GRND          ! net radiation over ground                   ! ||   ||
REAL,DIMENSION(1)  :: ZH_GRND           ! sensible heat flux over ground              ! ||   ||
REAL,DIMENSION(1)  :: ZLE_GRND          ! latent heat flux over ground                ! ||   ||
REAL,DIMENSION(1)  :: ZGFLUX_GRND       ! flux through the ground                     ! ||   ||
REAL,DIMENSION(1)  :: ZRN_TOWN          ! net radiation over town                     ! ||   ||
REAL,DIMENSION(1)  :: ZH_TOWN           ! sensible heat flux over town                ! ||   ||
REAL,DIMENSION(1)  :: ZLE_TOWN          ! latent heat flux over town                  ! ||   ||
REAL,DIMENSION(1)  :: ZGFLUX_TOWN       ! flux through the ground for town            ! ||   ||
REAL,DIMENSION(1)  :: ZEVAP_TOWN        ! evaporation (kg/m2/s)                       ! ||   ||
! Atmospheric Forcing variables                                                       ! ||   ||
REAL, DIMENSION(:,:), ALLOCATABLE :: ZTA    ! air temperature forcing (K)             ! ||   ||
REAL, DIMENSION(:,:), ALLOCATABLE :: ZQA    ! air humidity forcing (kg/m3)            ! ||   ||
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWIND  ! wind speed (m/s)                        ! ||   ||
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSCA_SW! diffuse solar radiation (on hor surf)   ! ||   ||
REAL, DIMENSION(:,:), ALLOCATABLE :: ZDIR_SW! direct  solar radiation (on hor surf)   ! ||   ||
REAL, DIMENSION(:,:), ALLOCATABLE :: ZLW    ! longwave radiation (on horizontal surf) ! ||   ||
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSNOW  ! snow precipitation  (kg/m2/s)           ! ||   ||
REAL, DIMENSION(:,:), ALLOCATABLE :: ZRAIN  ! liquid precipitation  (kg/m2/s)         !\\     //
REAL, DIMENSION(:,:), ALLOCATABLE :: ZPS    ! pressure at forcing level  (Pa)         ! \\   //
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCO2   ! CO2 concentration in the air  (kg/m3)   !  \\ //
REAL, DIMENSION(:,:), ALLOCATABLE :: ZDIR   ! wind direction                          !   \\/
REAL, DIMENSION(1)  :: ZCOEF                ! work array                              !
REAL, DIMENSION(1)  :: ZF1_o_B              ! Coefficient for sky model               ! ||   ||
REAL, DIMENSION(1,1):: ZTDIR_SW             ! total direct SW                         ! ||   ||
REAL, DIMENSION(1,1):: ZTSCA_SW             ! total diffuse SW                        ! ||   ||
REAL                :: ZTIME_BEG            ! Time at beginning of time step          ! ||   ||
REAL                :: ZTIME                ! Time at end       of time step          ! ||   ||
! -----------------------------------------------------------                         ! ||   ||
! Outputs                                                                             ! ||   ||
! -----------------------------------------------------------                         ! ||   ||
!                                                                                     ! ||   ||
CHARACTER(LEN=*), PARAMETER       :: T_ROOF1 = 'output/T_ROOF1.txt'                   ! ||   ||
CHARACTER(LEN=*), PARAMETER       :: T_CANYON = 'output/T_CANYON.txt'                 ! ||   ||
CHARACTER(LEN=*), PARAMETER       :: T_ROAD1 = 'output/T_ROAD1.txt'                   ! ||   ||
CHARACTER(LEN=*), PARAMETER       :: T_FLOOR1 = 'output/T_FLOOR1.txt'                   ! ||   ||
CHARACTER(LEN=*), PARAMETER       :: T_MASS1 = 'output/T_MASS1.txt'                   ! ||   ||
CHARACTER(LEN=*), PARAMETER       :: T_WALLA1= 'output/T_WALLA1.txt'                  ! ||   ||
CHARACTER(LEN=*), PARAMETER       :: T_WALLB1= 'output/T_WALLB1.txt'                  ! ||   ||
CHARACTER(LEN=*), PARAMETER       :: TI_BLD = 'output/TI_BLD.txt'                     ! ||   ||
CHARACTER(LEN=*), PARAMETER       :: Q_CANYON = 'output/Q_CANYON.txt'                 ! ||   ||
CHARACTER(LEN=*), PARAMETER       :: U_CANYON = 'output/U_CANYON.txt'                 ! ||   ||
CHARACTER(LEN=*), PARAMETER       :: H_TOWN = 'output/H_TOWN.txt'                     !\\     //
CHARACTER(LEN=*), PARAMETER       :: LE_TOWN = 'output/LE_TOWN.txt'                   ! \\   //
CHARACTER(LEN=*), PARAMETER       :: RN_TOWN = 'output/RN_TOWN.txt'                   !  \\ //
CHARACTER(LEN=*), PARAMETER       :: HVAC_COOL = 'output/HVAC_COOL.txt'               !   \\/
CHARACTER(LEN=*), PARAMETER       :: HVAC_HEAT = 'output/HVAC_HEAT.txt'               !
CHARACTER(LEN=*), PARAMETER       :: COP = 'output/COP.txt'
!
CHARACTER(LEN=*), PARAMETER       :: DIR_CANYON = 'output/DIR_CANYON.txt'
CHARACTER(LEN=*), PARAMETER       :: ALB_TOWN = 'output/ALB_TOWN.txt'
CHARACTER(LEN=*), PARAMETER       :: EMIS_TOWN = 'output/EMIS_TOWN.txt'
CHARACTER(LEN=*), PARAMETER       :: EVAP_TOWN = 'output/EVAP_TOWN.txt'
CHARACTER(LEN=*), PARAMETER       :: GFLUX_TOWN = 'output/GFLUX_TOWN.txt'
CHARACTER(LEN=*), PARAMETER       :: TS_TOWN = 'output/TS_TOWN.txt'
CHARACTER(LEN=*), PARAMETER       :: Q_TOWN = 'output/Q_TOWN.txt'
CHARACTER(LEN=*), PARAMETER       :: USTAR_TOWN = 'output/USTAR_TOWN.txt'
CHARACTER(LEN=*), PARAMETER       :: THER_PROD_PANEL = 'output/THER_PROD_PANEL.txt'
CHARACTER(LEN=*), PARAMETER       :: PHOT_PROD_PANEL = 'output/PHOT_PROD_PANEL.txt'
CHARACTER(LEN=*), PARAMETER       :: QF_TOWN = 'output/QF_TOWN.txt'
CHARACTER(LEN=*), PARAMETER       :: DQS_TOWN = 'output/DQS_TOWN.txt'
!
CHARACTER(LEN=*), PARAMETER       :: SWD = 'output/SWD.txt'
CHARACTER(LEN=*), PARAMETER       :: SWU = 'output/SWU.txt'
CHARACTER(LEN=*), PARAMETER       :: LWD = 'output/LWD.txt'
CHARACTER(LEN=*), PARAMETER       :: LWU = 'output/LWU.txt'
!
CHARACTER(LEN=*), PARAMETER       :: Forc_TA    = 'output/Forc_TA.txt'
CHARACTER(LEN=*), PARAMETER       :: Forc_PS    = 'output/Forc_PS.txt'
CHARACTER(LEN=*), PARAMETER       :: Forc_QA    = 'output/Forc_QA.txt'
CHARACTER(LEN=*), PARAMETER       :: Forc_WIND  = 'output/Forc_WIND.txt'
CHARACTER(LEN=*), PARAMETER       :: Forc_RAIN  = 'output/Forc_RAIN.txt'
CHARACTER(LEN=*), PARAMETER       :: Forc_SNOW  = 'output/Forc_SNOW.txt'

INTEGER                           :: ID_FORC
REAL, DIMENSION(1) :: ZSUMZEN
REAL, DIMENSION(1) :: ZSW
REAL                 :: ZBEGIN_TRAFFIC_TIME ! start traffic time (solar time, s)
REAL                 :: ZEND_TRAFFIC_TIME   ! end traffic time   (solar time, s)
INTEGER                           :: IYEAR2              ! current year at end of timestep(UTC)
INTEGER                           :: IMONTH2             ! current month at end of timestep(UTC)
INTEGER                           :: IDAY2               ! current day at end of timestep(UTC)
REAL                              :: ZTIME2              ! current time since start of the day at end of timestep (s)

logical :: exist ! TODO: testing remove

!===========================================================================
! Read from namelist file after default values are set.
!===========================================================================

OPEN(UNIT=99, FILE='input.nml', STATUS='old')
READ(99, NML = dimensions)

ALLOCATE(ZHC_ROOF    (1,NROOF_LAYER))
ALLOCATE(ZTC_ROOF    (1,NROOF_LAYER))
ALLOCATE(ZD_ROOF     (1,NROOF_LAYER))
ALLOCATE(ZT_ROOF     (1,NROOF_LAYER))

ALLOCATE(ZHC_ROAD    (1,NROAD_LAYER))
ALLOCATE(ZTC_ROAD    (1,NROAD_LAYER))
ALLOCATE(ZD_ROAD     (1,NROAD_LAYER))
ALLOCATE(ZT_ROAD     (1,NROAD_LAYER))

ALLOCATE(ZHC_WALL    (1,NWALL_LAYER))
ALLOCATE(ZTC_WALL    (1,NWALL_LAYER))
ALLOCATE(ZD_WALL     (1,NWALL_LAYER))
ALLOCATE(ZT_WALL   (1,NWALL_LAYER))
ALLOCATE(ZT_WALL_B   (1,NWALL_LAYER))

ALLOCATE(ZHC_FLOOR   (1,NFLOOR_LAYER))
ALLOCATE(ZTC_FLOOR   (1,NFLOOR_LAYER))
ALLOCATE(ZD_FLOOR    (1,NFLOOR_LAYER))
ALLOCATE(ZT_FLOOR    (1,NFLOOR_LAYER))
ALLOCATE(ZT_MASS     (1,NFLOOR_LAYER))

READ(99, NML = parameters)
CLOSE(99)

!===========================================================================
!===========================================================================
!===========================================================================
!===========================================================================
!===========================================================================
!===========================================================================
!                     *** DO NOT CHANGE VALUES BELOW ***
!===========================================================================
!===========================================================================
!===========================================================================
!===========================================================================
!===========================================================================
!===========================================================================
! -----------------------------------------------------------
! Inizializations (SYSTEM)
! -----------------------------------------------------------

ZT_WALL_B = ZT_WALL    ! wall layers temperatures
!

GSHAD_DAY = .FALSE. ! has shading been necessary this day ?
GNATVENT_NIGHT =.FALSE. ! has natural ventilation been necessary/possible this night ?

! irrigation not supported currently
LPAR_RD_IRRIG = .FALSE.
ZRD_START_MONTH = 6
ZRD_END_MONTH = 6
ZRD_START_HOUR = 6
ZRD_END_HOUR = 6
ZRD_24H_IRRIG = 0

IF (ZNATVENT(1) >= 0.0 .AND. ZNATVENT(1) < 0.5) THEN
  HNATVENT = 'NONE'
ELSEIF (ZNATVENT(1) >= 0.5 .AND. ZNATVENT(1) < 1.5) THEN
  HNATVENT = 'MANU'
ELSEIF (ZNATVENT(1) >= 1.5 .AND. ZNATVENT(1) <= 2.5) THEN
  HNATVENT = 'AUTO'        
ELSEIF (ZNATVENT(1) >= 2.5 .AND. ZNATVENT(1) <= 3.5) THEN
  HNATVENT = 'MECH'        
ELSE
  HNATVENT = 'NONE'        
ENDIF

IF (ZSHADE(1) >= 0.0 .AND. ZSHADE(1) < 0.5) THEN
  LSHADE = .FALSE.
ELSEIF (ZSHADE(1) >= 0.5 .AND. ZSHADE(1) <= 1.0) THEN
  LSHADE = .TRUE.
ELSE
  LSHADE = .FALSE.
ENDIF

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
!
! multi layer option
LCANOPY= .FALSE.  ! is canopy active ?  ** DO NOT CHANGE **
!
! initialization of physical constants
!
CALL INI_CSTS
!
!CALL INIT_SURFCONSPHY
!
CQSAT='NEW' ! saturation is computed relative to water above 0°C, and relative to ice below 0°C
!
!* various thresholds
!
XCISMIN = 0.5     ! minimum wind shear
XVMODMIN = 0.5    ! minimum wind speed
LALDTHRES = .FALSE.    ! activate aladin threshold for wind
XRIMAX = 0.2 ! Maximum richardson number for exchange coefficients computations
!
!* Open atmospheric forcing files
!
CALL OPEN_CLOSE_BIN_ASC_FORC('OPEN ','ASCII ','R')
!
! allocation of variables
!
CALL OL_ALLOC_ATM(1,1,1) ! INI, IBANDS, ISCAL
!
ZTIME = ZTIME_START
!
! sun position at run start
!
CALL SUNPOS(IYEAR, IMONTH, IDAY, ZTIME, ZLON, ZLAT, XTSUN, XZENITH, XAZIM)
!
! allocate local atmospheric variables
ALLOCATE(ZTA    (1,2))
ALLOCATE(ZQA    (1,2))
ALLOCATE(ZWIND  (1,2))
ALLOCATE(ZDIR_SW(1,2))
ALLOCATE(ZSCA_SW(1,2))
ALLOCATE(ZLW    (1,2))
ALLOCATE(ZSNOW  (1,2))
ALLOCATE(ZRAIN  (1,2))
ALLOCATE(ZPS    (1,2))
ALLOCATE(ZCO2   (1,2))
ALLOCATE(ZDIR   (1,2))
!
!* reads atmospheric forcing for first time-step
!
CALL OL_READ_ATM('ASCII ', 'ASCII ', 1,    &
                    ZTA,ZQA,ZWIND,ZDIR_SW,ZSCA_SW,ZLW,ZSNOW,ZRAIN,ZPS,&
                    ZCO2,ZDIR,.false. )
XCO2(:)  = ZCO2(:,1)
XRHOA(:) = ZPS(:,1) / ( ZTA(:,1)*XRD * ( 1.+((XRV/XRD)-1.)*ZQA(:,1) ) + XZREF(:)*XG )
!
! -----------------------------------------------------------
! Geometric parameters
! -----------------------------------------------------------
!
ZROAD         = (1. - ZBLD - ZGARDEN)
XZREF         = ZZREF
XUREF         = ZZREF
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
! -----------------------------------------------------------
! Window characteristics
! -----------------------------------------------------------
!
CALL WINDOW_DATA_STRUCT(1, ZSHGC, PU_WIN, ZALB_WIN, ZABS_WIN, ZUGG_WIN, ZTRAN_WIN)
!
! -----------------------------------------------------------
! Default at first time-step : no snow
! -----------------------------------------------------------
!
ZWSNOW_ROOF = 0.
ZTSNOW_ROOF = XUNDEF
ZRSNOW_ROOF = XUNDEF
ZASNOW_ROOF = XUNDEF
ZESNOW_ROOF = XUNDEF
ZTSSNOW_ROOF= XUNDEF
ZWSNOW_ROAD = 0.
ZTSNOW_ROAD = XUNDEF
ZRSNOW_ROAD = XUNDEF
ZASNOW_ROAD = XUNDEF
ZESNOW_ROAD = XUNDEF
ZTSSNOW_ROAD= XUNDEF
! -----------------------------------------------------------
! Outputs
! -----------------------------------------------------------
!
OPEN(UNIT=13, FILE = T_ROOF1,   ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=14, FILE = T_CANYON,  ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=15, FILE = T_ROAD1,   ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=16, FILE = T_WALLA1,  ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=17, FILE = T_WALLB1,  ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=18, FILE = TI_BLD,    ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=19, FILE = Q_CANYON,  ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=21, FILE = U_CANYON,  ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=22, FILE = H_TOWN,    ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=23, FILE = LE_TOWN,   ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=24, FILE = RN_TOWN,   ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=25, FILE = HVAC_COOL, ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=26, FILE = HVAC_HEAT, ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=37, FILE = COP,       ACCESS = 'APPEND',STATUS = 'REPLACE')
!
OPEN(UNIT=27, FILE = DIR_CANYON, ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=28, FILE = ALB_TOWN,   ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=29, FILE = EMIS_TOWN,  ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=30, FILE = EVAP_TOWN,  ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=31, FILE = GFLUX_TOWN, ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=32, FILE = TS_TOWN,    ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=33, FILE = Q_TOWN,     ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=34, FILE = USTAR_TOWN, ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=35, FILE = THER_PROD_PANEL, ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=36, FILE = PHOT_PROD_PANEL, ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=38, FILE = QF_TOWN, ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=39, FILE = DQS_TOWN, ACCESS = 'APPEND',STATUS = 'REPLACE')
!
OPEN(UNIT=40, FILE = SWD, ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=41, FILE = SWU, ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=42, FILE = LWD, ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=43, FILE = LWU, ACCESS = 'APPEND',STATUS = 'REPLACE')

OPEN(UNIT=50, FILE = Forc_TA, ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=51, FILE = Forc_PS, ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=52, FILE = Forc_QA, ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=53, FILE = Forc_WIND, ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=54, FILE = Forc_RAIN, ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=55, FILE = Forc_SNOW, ACCESS = 'APPEND',STATUS = 'REPLACE')

OPEN(UNIT=56, FILE = T_FLOOR1,   ACCESS = 'APPEND',STATUS = 'REPLACE')
OPEN(UNIT=57, FILE = T_MASS1,   ACCESS = 'APPEND',STATUS = 'REPLACE')
!
! -----------------------------------------------------------
! Temporal loops
! -----------------------------------------------------------
!
DO JFORC_STEP= 1,INB_STEP_ATM
   WRITE(*,FMT='(I5,A1,I5)') JFORC_STEP,'/',INB_STEP_ATM
   !
   ZTS_ROOF   = 0.
   ZTS_ROAD   = 0.
   ZTS_WALL_A = 0.
   ZTS_WALL_B = 0.
   ! read Forcing
   CALL OL_READ_ATM('ASCII ', 'ASCII ', JFORC_STEP,    &
                    ZTA,ZQA,ZWIND,ZDIR_SW,ZSCA_SW,ZLW,ZSNOW,ZRAIN,ZPS,&
                    ZCO2,ZDIR,.false. )
   !
  !COMPUTE SUM ZENITH angle between 2 timestepA
  ZSUMZEN(:)=0.0
  DO JSURF_STEP = 1,INB_ATM
    IDAY2  = IDAY
    ZTIME2 = ZTIME + (JSURF_STEP-1.)*XTSTEP_SURF
    IF (ZTIME2>86400.) THEN
      ZTIME2 = ZTIME2-86400
      IDAY2  = IDAY+1
    ENDIF
    CALL SUNPOS(IYEAR, IMONTH, IDAY2, ZTIME+(JSURF_STEP-1.)*XTSTEP_SURF, &
                ZLON, ZLAT, XTSUN, XZENITH, XAZIM)
    !
    ZSUMZEN(:)= ZSUMZEN(:) + MAX(COS(XZENITH(:)+0.1),0.)/(INB_ATM*1.0)
    !
  ENDDO
  WHERE ( ZSUMZEN<0.01 ) ZSUMZEN = 0.0

   !
   DO JSURF_STEP=1,INB_ATM
      !
      ! time interpolation of the forcing
      !
      !
      CALL SUNPOS(IYEAR, IMONTH, IDAY, ZTIME, ZLON, ZLAT, XTSUN, XZENITH, XAZIM)
      IYEAR2 = IYEAR
      IMONTH2= IMONTH
      IDAY2  = IDAY
      ZTIME2 = ZTIME+XTSTEP_SURF
      CALL ADD_FORECAST_TO_DATE_SURF(IYEAR2, IMONTH2, IDAY2, ZTIME2)
      CALL SUNPOS(IYEAR2, IMONTH2, IDAY2, ZTIME2, ZLON, ZLAT, XTSUN, XZENITH2, XAZIM)
       !
       TPTIME%TIME= ZTIME
       TPTIME%TDATE%YEAR =IYEAR
       TPTIME%TDATE%MONTH=IMONTH
       TPTIME%TDATE%DAY  =IDAY

      !
      !interpolation between beginning and end of current forcing time step
       ID_FORC = 1
    CALL OL_TIME_INTERP_ATM(JSURF_STEP,INB_ATM,                      &
                            ZTA(:,ID_FORC),ZTA(:,ID_FORC+1),         &
                            ZQA(:,ID_FORC),ZQA(:,ID_FORC+1),         &
                            ZWIND(:,ID_FORC),ZWIND(:,ID_FORC+1),     &
                            ZDIR_SW(:,ID_FORC),ZDIR_SW(:,ID_FORC+1), &
                            ZSCA_SW(:,ID_FORC),ZSCA_SW(:,ID_FORC+1), &
                            ZLW(:,ID_FORC),ZLW(:,ID_FORC+1),         &
                            ZSNOW(:,ID_FORC+1),ZRAIN(:,ID_FORC+1),   &
                            ZPS(:,ID_FORC),ZPS(:,ID_FORC+1),         &
                            ZCO2(:,ID_FORC),ZCO2(:,ID_FORC+1),       &
                            ZDIR(:,ID_FORC),ZDIR(:,ID_FORC+1),       &
                            XZENITH+0.1,ZSUMZEN                      )

      ZSW(:) = 0.
      DO JLOOP=1,SIZE(XDIR_SW,2)
        ZSW(:) = ZSW(:) + XDIR_SW(:,JLOOP) + XSCA_SW(:,JLOOP)
      END DO
      WHERE (ZSW(:)>0.)
        XZENITH  = MIN (XZENITH ,XPI/2.-0.01)
        XZENITH2 = MIN (XZENITH2,XPI/2.-0.01)
      ELSEWHERE
        XZENITH  = MAX (XZENITH ,XPI/2.)
        XZENITH2 = MAX (XZENITH2,XPI/2.)
      END WHERE
    
      ! updates time
      ZTIMEC= ZTIMEC+XTSTEP_SURF

       ! Exner functions
       !
       ZEXNS = (XPS/XP00)**(XRD/XCPD)
       ZEXNA = (XPA/XP00)**(XRD/XCPD)
       !
       ! specific humidity (conversion from kg/m3 to kg/kg)
       !
       ZQA_KGKG(:) = XQA(:) / XRHOA(:)

       ZTDIR_SW = XDIR_SW(1,1)
       ZTSCA_SW = XSCA_SW(1,1)
       KSW = 1  ! only one spectral band here
       !
       ! Sky model for diffuse radiation
       !
       !add directionnal contrib from scattered radiation
       CALL CIRCUMSOLAR_RAD(XDIR_SW(:,1), XSCA_SW(:,1), XZENITH, ZF1_o_B)
       ZTDIR_SW(:,1) = XDIR_SW(:,1) + XSCA_SW(:,1) * ZF1_o_B
       ZTSCA_SW(:,1) = XSCA_SW(:,1) * (1. - ZF1_o_B)


       ZBEGIN_TRAFFIC_TIME = 21600.
       ZEND_TRAFFIC_TIME   = 64800.
       !
       WHERE( XTSUN>ZBEGIN_TRAFFIC_TIME  .AND.  XTSUN<ZEND_TRAFFIC_TIME  )
         XH_TRAFFIC  (:) = ZH_TRAFFIC(:)
         XLE_TRAFFIC (:) = ZLE_TRAFFIC(:)
       ELSEWHERE
         XH_TRAFFIC  (:) = 0.
         XLE_TRAFFIC (:) = 0.   
       END WHERE

       !
       !
       ! Wind speed
       ZVMOD = SQRT(XU**2+XV**2)
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

!*****************************************************************************
!*****************************************************************************
!                  Call of physical routines of TEB is here                  !
!*****************************************************************************
!*****************************************************************************

    CALL TEB_GARDEN_STRUCT (LGARDEN, LGREENROOF, LSOLAR_PANEL,                &
                     HZ0H, HIMPLICIT_WIND, HROAD_DIR, HWALL_OPT, TPTIME,      &
                     XTSUN, ZT_CANYON, ZQ_CANYON, ZU_CANYON,                  &
                     ZT_LOWCAN, ZQ_LOWCAN, ZU_LOWCAN, ZZ_LOWCAN,              &
                     ZTI_BLD,                                                 &
                     ZT_ROOF, ZT_ROAD, ZT_WALL, ZT_WALL_B,                  &
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
                     XH_TRAFFIC, XLE_TRAFFIC, ZH_INDUSTRY, ZLE_INDUSTRY,      &
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
                     ZHU_BLD, ZTIME, LSHADE, ZSHADE, GSHAD_DAY,           &
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
!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
!

    ZTIME = ZTIME + XTSTEP_SURF
    CALL ADD_FORECAST_TO_DATE_SURF(IYEAR, IMONTH, IDAY, ZTIME)


   END DO

!  Instantaneous diagnostics
ZTS_ROOF = ZT_ROOF(1,1)
ZTS_ROAD = ZT_ROAD(1,1)
ZTS_WALL_A = ZT_WALL(1,1)
ZTS_WALL_B = ZT_WALL_B(1,1)

! Calucalte the total (direct + diffuse radiaitons)town albedo
IF (ZTDIR_SW(1,1) == 0.) THEN
  ZALB_TOWN  = 0.
ELSE
  ZALB_TOWN = (ZDIR_ALB_TOWN * ZTDIR_SW(1,1) + ZSCA_ALB_TOWN * ZTSCA_SW(1,1)) / (ZTDIR_SW(1,1) + ZTSCA_SW(1,1))
END IF
! TEB does no modify the wind direction
ZDIR_CANYON = ZDIR(1,1)
! Town specific humidity assumed to be same as canyon humidity
ZQ_TOWN = ZQ_CANYON

   !
   WRITE(13,*) ZTS_ROOF
   WRITE(14,*) ZT_CANYON
   WRITE(15,*) ZTS_ROAD
   WRITE(16,*) ZTS_WALL_A
   WRITE(17,*) ZTS_WALL_B
   WRITE(18,*) ZTI_BLD
   WRITE(19,*) ZQ_CANYON
   WRITE(21,*) ZU_CANYON
   WRITE(22,*) ZH_TOWN
   WRITE(23,*) ZLE_TOWN
   WRITE(24,*) ZRN_TOWN
   IF (CBEM=='BEM') THEN
     ! The heating and cooling energy demand are converted
     ! from W/m²(bld) to W/m²(urb).
     ! FIXME: this does not seem correct anymore, why?
     ! it seems that ZHVAC_COOL is already W/m²(urb)...  
     ZHVAC_COOL_TOT = ZBLD * ZHVAC_COOL
     ZHVAC_HEAT_TOT = ZBLD * ZHVAC_HEAT
     WRITE(25,*) ZHVAC_COOL
     WRITE(26,*) ZHVAC_HEAT
     WRITE(37,*) ZCOP
   END IF
!
    WRITE(27,*) ZDIR_CANYON
    WRITE(28,*) ZALB_TOWN
    WRITE(29,*) ZEMIS_TOWN
    WRITE(30,*) ZEVAP_TOWN
    WRITE(31,*) ZGFLUX_TOWN
    WRITE(32,*) ZTS_TOWN
    WRITE(33,*) ZQ_TOWN
    WRITE(34,*) ZUSTAR_TOWN
    WRITE(38,*) ZQF_TOWN
    WRITE(39,*) ZDQS_TOWN
    IF (LSOLAR_PANEL) THEN
      WRITE(35,*) ZTHER_PROD_PANEL
      WRITE(36,*) ZPHOT_PROD_PANEL
    END IF
    WRITE(40,*) ZTDIR_SW(1,1) + ZTSCA_SW(1,1)       
                      ! PDIR_SW(:,JSWB)                    + PSCA_SW(:,JSWB)
    WRITE(41,*) ZTDIR_SW(1,1) * ZDIR_ALB_TOWN(1) + ZTSCA_SW(1,1) * ZSCA_ALB_TOWN(1)
                      ! PDIR_SW(:,JSWB) * PDIR_ALB(:,JSWB) + PSCA_SW(:,JSWB) * PSCA_ALB(:,JSWB) 
    WRITE(42,*) XLW
                      ! PLW(:)
    WRITE(43,*) ZEMIS_TOWN * XSTEFAN *ZTS_TOWN**4 + (1.-ZEMIS_TOWN)*XLW
      ! PEMIS(:)*XSTEFAN*PTRAD(:)**4 + (1.-PEMIS(:))*PLW(:)

    ! Interpolated forcing variables
    WRITE(50,*) XTA
    WRITE(51,*) XPS
    WRITE(52,*) XQA
    WRITE(53,*) ZVMOD
    WRITE(54,*) XRAIN
    WRITE(55,*) XSNOW

    WRITE(56,*) ZT_FLOOR(1,1)
    WRITE(57,*) ZT_MASS(1,1)

!
END DO
!
! --------------------------------------------------------------------------------------
!
!
!  DEALLOCATE variables
DEALLOCATE(ZTA)
DEALLOCATE(ZQA)
DEALLOCATE(ZWIND)
DEALLOCATE(ZDIR_SW)
DEALLOCATE(ZSCA_SW)
DEALLOCATE(ZLW)
DEALLOCATE(ZSNOW)
DEALLOCATE(ZRAIN)
DEALLOCATE(ZPS)
DEALLOCATE(ZCO2)
DEALLOCATE(ZDIR)
!
DEALLOCATE(ZHC_WALL)
DEALLOCATE(ZTC_WALL)
DEALLOCATE(ZD_WALL)
DEALLOCATE(ZT_WALL)
DEALLOCATE(ZT_WALL_B)
DEALLOCATE(ZHC_ROOF)
DEALLOCATE(ZTC_ROOF)
DEALLOCATE(ZD_ROOF)
DEALLOCATE(ZT_ROOF)
DEALLOCATE(ZHC_ROAD)
DEALLOCATE(ZTC_ROAD)
DEALLOCATE(ZD_ROAD)
DEALLOCATE(ZT_ROAD)
!
CALL OPEN_CLOSE_BIN_ASC_FORC('CLOSE ','ASCII ','R')
CLOSE(13)
CLOSE(14)
CLOSE(15)
CLOSE(16)
CLOSE(17)
CLOSE(18)
CLOSE(19)
CLOSE(21)
CLOSE(22)
CLOSE(23)
CLOSE(24)
CLOSE(25)
CLOSE(26)
CLOSE(37)
CLOSE(27)
CLOSE(28)
CLOSE(29)
CLOSE(30)
CLOSE(31)
CLOSE(32)
CLOSE(33)
CLOSE(34)
CLOSE(35)
CLOSE(36)
CLOSE(38)
CLOSE(39)
CLOSE(40)
CLOSE(41)
CLOSE(42)
CLOSE(43)

CLOSE(50)
CLOSE(51)
CLOSE(52)
CLOSE(53)
CLOSE(54)
CLOSE(55)

CLOSE(56)
CLOSE(57)


!
    WRITE(*,*) ' '
    WRITE(*,*) '    --------------------------'
    WRITE(*,*) '    |  DRIVER ENDS CORRECTLY |'
    WRITE(*,*) '    --------------------------'
    WRITE(*,*) ' '
!
! --------------------------------------------------------------------------------------
!
END PROGRAM DRIVER
