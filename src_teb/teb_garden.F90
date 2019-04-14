!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE TEB_GARDEN (OGARDEN, OGREENROOF, OSOLAR_PANEL,                 &
                     HZ0H, HIMPLICIT_WIND, HROAD_DIR, HWALL_OPT, TPTIME,      &
                     PTSUN, PT_CANYON, PQ_CANYON, PU_CANYON,                  &
                     PT_LOWCAN, PQ_LOWCAN, PU_LOWCAN, PZ_LOWCAN, PTI_BLD,     &
                     PT_ROOF, PT_ROAD, PT_WALL_A, PT_WALL_B,                  &
                     PWS_ROOF,PWS_ROAD,                                       &
                     HSNOW_ROOF,                                              &
                     PWSNOW_ROOF, PTSNOW_ROOF, PRSNOW_ROOF, PASNOW_ROOF,      &
                     PTSSNOW_ROOF, PESNOW_ROOF,                               &
                     HSNOW_ROAD,                                              &
                     PWSNOW_ROAD, PTSNOW_ROAD, PRSNOW_ROAD, PASNOW_ROAD,      &
                     PTSSNOW_ROAD, PESNOW_ROAD,                               &
                     PPEW_A_COEF, PPEW_B_COEF,                                &
                     PPEW_A_COEF_LOWCAN, PPEW_B_COEF_LOWCAN,                  &
                     PPS, PPA, PEXNS, PEXNA,                                  &
                     PTA, PQA, PRHOA, PCO2,                                   &
                     PLW_RAD, PDIR_SW, PSCA_SW, PSW_BANDS, KSW,               &
                     PZENITH, PAZIM,                                          &
                     PRR, PSR,                                                &
                     PZREF, PUREF, PVMOD,                                     &
                     PH_TRAFFIC, PLE_TRAFFIC, PH_INDUSTRY, PLE_INDUSTRY,      &
                     PTSTEP, PZ0_TOWN, PBLD, PGARDEN, PROAD_DIR, PROAD,       &
                     PFRAC_GR,                                                &
                     PBLD_HEIGHT, PWALL_O_HOR, PCAN_HW_RATIO,                 &
                     PROAD_O_GRND, PGARDEN_O_GRND, PWALL_O_GRND,              &
                     PALB_ROOF, PEMIS_ROOF,                                   &
                     PHC_ROOF,PTC_ROOF,PD_ROOF,                               &
                     PALB_ROAD, PEMIS_ROAD, PSVF_ROAD,                        &
                     PHC_ROAD,PTC_ROAD,PD_ROAD,                               &
                     PALB_WALL, PEMIS_WALL, PSVF_WALL,                        &
                     PSVF_GARDEN,                                             &
                     PHC_WALL,PTC_WALL,PD_WALL,                               &
                     PRN_ROOF, PH_ROOF, PLE_ROOF, PLEW_ROOF, PGFLUX_ROOF,     &
                     PRUNOFF_ROOF,                                            &
                     PRN_ROAD, PH_ROAD, PLE_ROAD, PLEW_ROAD, PGFLUX_ROAD,     &
                     PRUNOFF_ROAD,                                            &
                     PRN_WALL_A, PH_WALL_A, PLE_WALL_A, PGFLUX_WALL_A,        &
                     PRN_WALL_B, PH_WALL_B, PLE_WALL_B, PGFLUX_WALL_B,        &
                     PRN_GARDEN,PH_GARDEN,PLE_GARDEN, PGFLUX_GARDEN,          &
                     PRUNOFF_GARDEN, PDRAIN_GARDEN, PIRRIG_GARDEN,            &
                     PRN_GREENROOF,PH_GREENROOF,PLE_GREENROOF, PGFLUX_GREENROOF, &
                     PRN_STRLROOF,PH_STRLROOF,PLE_STRLROOF, PGFLUX_STRLROOF,  &
                     PRUNOFF_STRLROOF,                                        &
                     PRN_BLT,PH_BLT,PLE_BLT, PGFLUX_BLT,                      &
                     PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,    &
                     PMELT_ROOF,                                              &
                     PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,    &
                     PMELT_ROAD,                                              &
                     PRN_GRND, PH_GRND, PLE_GRND, PGFLUX_GRND,                &
                     PRN_TOWN, PH_TOWN, PLE_TOWN, PGFLUX_TOWN, PEVAP_TOWN,    &
                     PRUNOFF_TOWN, PSFCO2,                                    &
                     PUW_GRND, PUW_ROOF, PDUWDU_GRND, PDUWDU_ROOF,            &
                     PUSTAR_TOWN, PCD, PCDN, PCH_TOWN, PRI_TOWN,              &
                     PTS_TOWN, PEMIS_TOWN, PDIR_ALB_TOWN, PSCA_ALB_TOWN,      &
                     PRESA_TOWN, PDQS_TOWN, PQF_TOWN, PQF_BLD,                &
                     PFLX_BLD, PAC_ROAD, PAC_GARDEN, PAC_GREENROOF,           &
                     PAC_ROAD_WAT, PAC_GARDEN_WAT, PAC_GREENROOF_WAT,         &
                     PABS_SW_ROOF,PABS_LW_ROOF,                               &
                     PABS_SW_SNOW_ROOF,PABS_LW_SNOW_ROOF,                     &
                     PABS_SW_ROAD,PABS_LW_ROAD,                               &
                     PABS_SW_SNOW_ROAD,PABS_LW_SNOW_ROAD,                     &
                     PABS_SW_WALL_A,PABS_LW_WALL_A,                           &
                     PABS_SW_WALL_B,PABS_LW_WALL_B,                           &
                     PABS_SW_PANEL,PABS_LW_PANEL,                             &
                     PABS_SW_GARDEN,PABS_LW_GARDEN,                           &
                     PABS_SW_GREENROOF,PABS_LW_GREENROOF,                     &
                     PG_GREENROOF_ROOF,PRUNOFF_GREENROOF,PDRAIN_GREENROOF,    &
                     PIRRIG_GREENROOF, HCOOL_COIL, PF_WATER_COND, HHEAT_COIL, &
                     HNATVENT, KDAY, PAUX_MAX, PT_FLOOR,                      &
                     PT_MASS, PH_BLD_COOL, PT_BLD_COOL, PH_BLD_HEAT,          &
                     PLE_BLD_COOL, PLE_BLD_HEAT, PH_WASTE, PLE_WASTE,         &
                     PF_WASTE_CAN, PHVAC_COOL, PHVAC_HEAT, PQIN, PQIN_FRAD,   &
                     PQIN_FLAT, PGR, PEFF_HEAT, PINF,                         &
                     PTCOOL_TARGET, PTHEAT_TARGET, PHR_TARGET, PT_WIN2,       &
                     PQI_BLD, PV_VENT, PCAP_SYS_HEAT, PCAP_SYS_RAT, PT_ADP,   &
                     PM_SYS_RAT, PCOP_RAT, PCAP_SYS, PM_SYS, PCOP, PQ_SYS,    &
                     PT_SYS, PTR_SW_WIN, PFAN_POWER, PHC_FLOOR, PTC_FLOOR,    &
                     PD_FLOOR, PT_WIN1, PABS_SW_WIN, PABS_LW_WIN, PSHGC,      &
                     PSHGC_SH, PUGG_WIN, PALB_WIN, PABS_WIN, PEMIT_LW_FAC,    &
                     PEMIT_LW_GRND, PT_RAD_IND, PREF_SW_GRND, PREF_SW_FAC,    &
                     PHU_BLD, PTIME, OSHADE, OSHAD_DAY, ONATVENT_NIGHT, HBEM, &
                     PN_FLOOR, PWALL_O_BLD, PGLAZ_O_BLD, PMASS_O_BLD,         &
                     PFLOOR_HW_RATIO, PF_FLOOR_MASS, PF_FLOOR_WALL,           &
                     PF_FLOOR_WIN, PF_FLOOR_ROOF, PF_WALL_FLOOR, PF_WALL_MASS,&
                     PF_WALL_WIN, PF_WIN_FLOOR, PF_WIN_MASS, PF_WIN_WALL,     &
                     PF_MASS_FLOOR, PF_MASS_WALL, PF_MASS_WIN, OCANOPY,       &
                     PTRAN_WIN, HCH_BEM, PROUGH_ROOF, PROUGH_WALL, PF_WIN_WIN,&
                     OPAR_RD_IRRIG, PRD_START_MONTH, PRD_END_MONTH,           &
                     PRD_START_HOUR, PRD_END_HOUR, PRD_24H_IRRIG, PIRRIG_ROAD,&
                     PEMIS_PANEL, PALB_PANEL, PEFF_PANEL, PFRAC_PANEL,        &
                     PRESIDENTIAL,                                            &
                     PTHER_PROD_PANEL, PPHOT_PROD_PANEL, PPROD_PANEL,         &
                     PTHER_PROD_BLD  , PPHOT_PROD_BLD  , PPROD_BLD  ,         &
                     PTHER_PRODC_DAY, PH_PANEL, PRN_PANEL,                    &
                     PDT_RES, PDT_OFF,                                        &
                     PCUR_TCOOL_TARGET, PCUR_THEAT_TARGET, PCUR_QIN           )
!   ##########################################################################
!
!!****  *TEB_GARDEN*  
!!
!!    PURPOSE
!!    -------
!
!     
!!**  METHOD
!     ------
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!      A. Lemonsu          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!    Original    05/2009
!!                04/2012 add PTRAN_WIN
!!    modified    08/2012 HCH_BEM/ ROUGH_WALL and ROUGH_ROOF for buildind conv coef.
!!    modified    10/2012 add PF_WIN_WIN as arg
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
!
!
!
!
USE MODD_TYPE_DATE_SURF,    ONLY: DATE_TIME
USE MODD_CSTS,              ONLY: XTT, XSTEFAN
USE MODD_SURF_PAR,          ONLY: XUNDEF
USE MODD_SNOW_PAR,          ONLY: XEMISSN, XANSMAX
!
USE MODE_THERMOS
USE MODE_SURF_SNOW_FRAC
!
USE MODI_SOLAR_PANEL
USE MODI_GARDEN_PROPERTIES
USE MODI_GREENROOF_PROPERTIES
USE MODI_WINDOW_SHADING_AVAILABILITY
USE MODI_URBAN_SOLAR_ABS
USE MODI_URBAN_LW_COEF
USE MODI_GARDEN
USE MODI_GREENROOF
USE MODI_TEB
USE MODI_AVG_URBAN_FLUXES
USE MODI_BLD_OCC_CALENDAR
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    Declarations of arguments
!
!
!
!
 LOGICAL,              INTENT(IN)    :: OGARDEN           ! Flag to use a garden    model inside the canyon
 LOGICAL,              INTENT(IN)    :: OGREENROOF        ! Flag to use a greenroof model on roofs
 LOGICAL,              INTENT(IN)    :: OSOLAR_PANEL      ! Flag to use a Solar Panel model on roofs
 CHARACTER(LEN=6)    , INTENT(IN)    :: HZ0H              ! TEB option for z0h roof & road
!                                                         ! 'MASC95' : Mascart et al 1995
!                                                         ! 'BRUT82' : Brustaert     1982
!                                                         ! 'KAND07' : Kanda         2007
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND       ! wind implicitation option
!                                                         ! 'OLD' = direct
!                                                         ! 'NEW' = Taylor serie, order 1
 CHARACTER(LEN=4)    , INTENT(IN)    :: HROAD_DIR          ! TEB option for road direction
                                                          ! 'UNIF' : uniform roads
                                                          ! 'ORIE' : specified
                                                          !   road orientation 
 CHARACTER(LEN=4)    , INTENT(IN)    :: HWALL_OPT          ! TEB option for walls representation
                                                          ! 'UNIF' : uniform walls
                                                          ! 'TWO ' : 2 opposite  walls
TYPE(DATE_TIME)     , INTENT(IN)    :: TPTIME             ! current date and time from teb
REAL, DIMENSION(:),   INTENT(IN)    :: PTSUN              ! solar time   (s from midnight)

!                                                         
REAL, DIMENSION(:)  , INTENT(INOUT) :: PT_CANYON          ! canyon air temperature
REAL, DIMENSION(:)  , INTENT(INOUT) :: PQ_CANYON          ! canyon air specific humidity
REAL, DIMENSION(:)  , INTENT(IN)    :: PU_CANYON          ! canyon hor. wind
REAL, DIMENSION(:)  , INTENT(IN)    :: PU_LOWCAN          ! wind near the road
REAL, DIMENSION(:)  , INTENT(IN)    :: PT_LOWCAN          ! temp. near the road
REAL, DIMENSION(:)  , INTENT(IN)    :: PQ_LOWCAN          ! hum. near the road
REAL, DIMENSION(:)  , INTENT(IN)    :: PZ_LOWCAN          ! height of atm. var. near the road
REAL, DIMENSION(:)  , INTENT(INOUT) :: PTI_BLD            ! inside building temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_ROOF            ! roof layers temperatures
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_ROAD            ! road layers temperatures
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_WALL_A          ! wall layers temperatures
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_WALL_B          ! wall layers temperatures
REAL, DIMENSION(:)  , INTENT(INOUT) :: PWS_ROOF           ! roof water reservoir
REAL, DIMENSION(:)  , INTENT(INOUT) :: PWS_ROAD           ! road water reservoir
 CHARACTER(LEN=*)    , INTENT(IN)    :: HSNOW_ROOF         ! snow roof scheme 'NONE', 'D95 ', '1-L '
 CHARACTER(LEN=*)    , INTENT(IN)    :: HSNOW_ROAD         ! snow road scheme 'NONE', 'D95 ', '1-L '
REAL, DIMENSION(:,:), INTENT(INOUT) :: PWSNOW_ROOF        ! snow layers reservoir
REAL, DIMENSION(:,:), INTENT(INOUT) :: PTSNOW_ROOF        ! snow layers temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRSNOW_ROOF        ! snow layers density
REAL, DIMENSION(:)  , INTENT(INOUT) :: PASNOW_ROOF        ! snow albedo
REAL, DIMENSION(:)  , INTENT(INOUT) :: PESNOW_ROOF        ! snow emissivity
REAL, DIMENSION(:)  , INTENT(INOUT) :: PTSSNOW_ROOF       ! snow surface temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PWSNOW_ROAD        ! snow layers reservoir
REAL, DIMENSION(:,:), INTENT(INOUT) :: PTSNOW_ROAD        ! snow layers temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRSNOW_ROAD        ! snow layers density
REAL, DIMENSION(:)  , INTENT(INOUT) :: PASNOW_ROAD        ! snow albedo
REAL, DIMENSION(:)  , INTENT(INOUT) :: PESNOW_ROAD        ! snow emissivity
REAL, DIMENSION(:)  , INTENT(INOUT) :: PTSSNOW_ROAD       ! snow surface temperature
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEW_A_COEF        ! implicit coefficients
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEW_B_COEF        ! for wind coupling
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEW_A_COEF_LOWCAN ! implicit coefficients for wind coupling
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEW_B_COEF_LOWCAN ! between low canyon wind and road
REAL, DIMENSION(:)  , INTENT(IN)    :: PPS                ! pressure at the surface
REAL, DIMENSION(:)  , INTENT(IN)    :: PPA                ! pressure at the first atmospheric level
REAL, DIMENSION(:)  , INTENT(IN)    :: PEXNS              ! surface exner function
REAL, DIMENSION(:)  , INTENT(IN)    :: PTA                ! temperature at the lowest level
REAL, DIMENSION(:)  , INTENT(IN)    :: PQA                ! specific humidity at the lowest level
REAL, DIMENSION(:)  , INTENT(IN)    :: PVMOD              ! module of the horizontal wind
REAL, DIMENSION(:)  , INTENT(IN)    :: PEXNA              ! exner function at the lowest level
REAL, DIMENSION(:)  , INTENT(IN)    :: PRHOA              ! air density at the lowest level
REAL, DIMENSION(:)  , INTENT(IN)    :: PCO2               ! CO2 concentration in the air    (kg/m3)
REAL, DIMENSION(:)  , INTENT(IN)    :: PLW_RAD            ! atmospheric infrared radiation
REAL, DIMENSION(:,:), INTENT(IN)    :: PDIR_SW            ! incoming direct solar rad on an horizontal surface
REAL, DIMENSION(:,:), INTENT(IN)    :: PSCA_SW            ! scattered incoming solar rad.
REAL, DIMENSION(:)  , INTENT(IN)    :: PSW_BANDS          ! mean wavelength of each shortwave band (m)
INTEGER,              INTENT(IN)    :: KSW                ! number of short-wave spectral bands
REAL, DIMENSION(:)  , INTENT(IN)    :: PZENITH            ! solar zenithal angle
REAL, DIMENSION(:)  , INTENT(IN)    :: PAZIM              ! solar azimuthal angle
                                                          ! (radian form N, clockwise)
REAL, DIMENSION(:)  , INTENT(IN)    :: PRR                ! rain rate
REAL, DIMENSION(:)  , INTENT(IN)    :: PSR                ! snow rate
REAL, DIMENSION(:)  , INTENT(IN)    :: PH_TRAFFIC         ! anthropogenic sensible heat fluxes due to traffic
REAL, DIMENSION(:)  , INTENT(IN)    :: PLE_TRAFFIC        ! anthropogenic latent heat fluxes due to traffic
REAL, DIMENSION(:)  , INTENT(IN)    :: PH_INDUSTRY        ! anthropogenic sensible heat fluxes due to factories
REAL, DIMENSION(:)  , INTENT(IN)    :: PLE_INDUSTRY       ! anthropogenic latent heat fluxes due to factories
REAL, DIMENSION(:)  , INTENT(IN)    :: PZREF              ! reference height of the first atm level (temperature)
REAL, DIMENSION(:)  , INTENT(IN)    :: PUREF              ! reference height of the first atm level (wind)
REAL                , INTENT(IN)    :: PTSTEP             ! time step
REAL, DIMENSION(:)  , INTENT(IN)    :: PZ0_TOWN           ! town roughness length for momentum
REAL, DIMENSION(:)  , INTENT(IN)    :: PBLD               ! fraction of buildings
REAL, DIMENSION(:)  , INTENT(IN)    :: PGARDEN            ! fraction of green areas
REAL, DIMENSION(:)  , INTENT(IN)    :: PROAD_DIR          ! road direction (deg from North, clockwise)
REAL, DIMENSION(:)  , INTENT(IN)    :: PROAD              ! fraction of roads
REAL, DIMENSION(:)  , INTENT(IN)    :: PFRAC_GR           ! fraction of green roofs
REAL, DIMENSION(:)  , INTENT(IN)    :: PBLD_HEIGHT        ! buildings h
REAL, DIMENSION(:)  , INTENT(IN)    :: PWALL_O_HOR        ! wall surf. / hor. surf.
REAL, DIMENSION(:)  , INTENT(IN)    :: PCAN_HW_RATIO      ! canyon    h/W
REAL, DIMENSION(:)  , INTENT(IN)    :: PROAD_O_GRND       ! road surf.   / (road+garden surf.) 
REAL, DIMENSION(:)  , INTENT(IN)    :: PGARDEN_O_GRND     ! garden surf. / (road+garden surf.) 
REAL, DIMENSION(:)  , INTENT(IN)    :: PWALL_O_GRND       ! wall surf.   / (road+garden surf.) 
REAL, DIMENSION(:)  , INTENT(IN)    :: PALB_ROOF          ! roof albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PEMIS_ROOF         ! roof emissivity
REAL, DIMENSION(:,:), INTENT(IN)    :: PHC_ROOF           ! heat capacity for roof layers
REAL, DIMENSION(:,:), INTENT(IN)    :: PTC_ROOF           ! thermal conductivity for roof layers
REAL, DIMENSION(:,:), INTENT(IN)    :: PD_ROOF            ! depth of roof layers
REAL, DIMENSION(:)  , INTENT(IN)    :: PALB_ROAD          ! road albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PEMIS_ROAD         ! road emissivity
REAL, DIMENSION(:,:), INTENT(IN)    :: PHC_ROAD           ! heat capacity for road layers
REAL, DIMENSION(:,:), INTENT(IN)    :: PTC_ROAD           ! thermal conductivity for road layers
REAL, DIMENSION(:,:), INTENT(IN)    :: PD_ROAD            ! depth of road layers
REAL, DIMENSION(:)  , INTENT(IN)    :: PSVF_ROAD          ! road sky view factor
REAL, DIMENSION(:)  , INTENT(IN)    :: PALB_WALL          ! wall albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PEMIS_WALL         ! wall emissivity
REAL, DIMENSION(:,:), INTENT(IN)    :: PHC_WALL           ! heat capacity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)    :: PTC_WALL           ! thermal conductivity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)    :: PD_WALL            ! depth of wall layers
REAL, DIMENSION(:)  , INTENT(IN)    :: PSVF_WALL          ! wall sky view factor
REAL, DIMENSION(:)  , INTENT(IN)    :: PSVF_GARDEN        ! green area sky view factor
     !
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_ROOF           ! net radiation over roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_ROOF            ! sensible heat flux over roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_ROOF           ! latent heat flux over roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLEW_ROOF          ! latent heat flux over roof (snow)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX_ROOF        ! flux through the roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRUNOFF_ROOF       ! runoff over the ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_ROAD           ! net radiation over road
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_ROAD            ! sensible heat flux over road
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_ROAD           ! latent heat flux over road
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLEW_ROAD          ! latent heat flux over road (snow)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX_ROAD        ! flux through the road
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRUNOFF_ROAD       ! runoff over the ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_WALL_A         ! net radiation over wall
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_WALL_A          ! sensible heat flux over wall
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_WALL_A         ! latent heat flux over wall
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX_WALL_A      ! flux through the wall
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_WALL_B         ! net radiation over wall
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_WALL_B          ! sensible heat flux over wall
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_WALL_B         ! latent heat flux over wall
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX_WALL_B      ! flux through the wall
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_GARDEN         ! net radiation over green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_GARDEN          ! sensible heat flux over green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_GARDEN         ! latent heat flux over green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX_GARDEN      ! flux through the green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRUNOFF_GARDEN     ! runoff over garden
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDRAIN_GARDEN      ! drainage over garden
REAL, DIMENSION(:)  , INTENT(OUT)   :: PIRRIG_GARDEN      ! water supply for irrigation over garden (kg/m2/s)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_GREENROOF      ! net radiation over greenroofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_GREENROOF       ! sensible heat flux over greenroofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_GREENROOF      ! latent heat flux over greenroofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX_GREENROOF   ! flux through the greenroofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_STRLROOF       ! net radiation over structural roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_STRLROOF        ! sensible heat flux over structural roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_STRLROOF       ! latent heat flux over structural roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX_STRLROOF    ! flux through the structural roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRUNOFF_STRLROOF   ! runoff aver structural roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_BLT            ! net radiation over built surf 
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_BLT             ! sensible heat flux over built surf 
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_BLT            ! latent heat flux over built surf 
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX_BLT         ! flux through the built surf 
!
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRNSNOW_ROOF       ! net radiation over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PHSNOW_ROOF        ! sensible heat flux over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLESNOW_ROOF       ! latent heat flux over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGSNOW_ROOF        ! flux under the snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PMELT_ROOF         ! snow melt
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRNSNOW_ROAD       ! net radiation over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PHSNOW_ROAD        ! sensible heat flux over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLESNOW_ROAD       ! latent heat flux over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGSNOW_ROAD        ! flux under the snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PMELT_ROAD         ! snow melt
!
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_GRND           ! net radiation over ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_GRND            ! sensible heat flux over ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_GRND           ! latent heat flux over ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX_GRND        ! flux through the ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_TOWN           ! net radiation over town
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_TOWN            ! sensible heat flux over town
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_TOWN           ! latent heat flux over town
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX_TOWN        ! flux through the ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PEVAP_TOWN         ! evaporation flux (kg/m2/s)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRUNOFF_TOWN       ! runoff over the ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PSFCO2             ! flux of CO2       (m/s*kg_CO2/kg_air)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PUW_GRND           ! momentum flux for ground built surf
REAL, DIMENSION(:)  , INTENT(OUT)   :: PUW_ROOF           ! momentum flux for roofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDUWDU_GRND        !
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDUWDU_ROOF        !
REAL, DIMENSION(:)  , INTENT(OUT)   :: PUSTAR_TOWN        ! friciton velocity over town
REAL, DIMENSION(:)  , INTENT(OUT)   :: PCD                ! town averaged drag coefficient
REAL, DIMENSION(:)  , INTENT(OUT)   :: PCDN               ! town averaged neutral drag coefficient
REAL, DIMENSION(:)  , INTENT(OUT)   :: PCH_TOWN           ! town averaged heat transfer coefficient
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRI_TOWN           ! town averaged Richardson number
REAL, DIMENSION(:)  , INTENT(OUT)   :: PTS_TOWN           ! town surface temperature
REAL, DIMENSION(:)  , INTENT(OUT)   :: PEMIS_TOWN         ! town equivalent emissivity
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDIR_ALB_TOWN      ! town equivalent direct albedo
REAL, DIMENSION(:)  , INTENT(OUT)   :: PSCA_ALB_TOWN      ! town equivalent diffuse albedo
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRESA_TOWN         ! town aerodynamical resistance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDQS_TOWN          ! heat storage inside town
REAL, DIMENSION(:)  , INTENT(OUT)   :: PQF_TOWN           ! total anthropogenic heat
REAL, DIMENSION(:)  , INTENT(OUT)   :: PQF_BLD            ! anthropogenic heat flux of domestic heating
REAL, DIMENSION(:)  , INTENT(OUT)   :: PFLX_BLD           ! heat flx from inside bld through its structure
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_ROAD           ! road conductance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_GARDEN         ! green area conductance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_GREENROOF      ! green roof conductance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_ROAD_WAT       ! road conductance for latent heat
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_GARDEN_WAT     ! green area conductance for latent heat
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_GREENROOF_WAT  ! green roof conductance for latent heat
!
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_SW_ROOF       ! absorbed solar rad by roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_SW_SNOW_ROOF  ! absorbed solar rad by snow on roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_LW_ROOF       ! absorbed IR rad by roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_LW_SNOW_ROOF  ! absorbed IR rad by snow on roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_SW_ROAD       ! absorbed solar rad by road
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_SW_SNOW_ROAD  ! absorbed solar rad by snow on road
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_LW_ROAD       ! absorbed IR rad by road
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_LW_SNOW_ROAD  ! absorbed IR rad by snow on road
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_SW_WALL_A     ! absorbed solar rad by wall
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_LW_WALL_A     ! absorbed IR rad by wall
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_SW_WALL_B     ! absorbed solar rad by wall
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_LW_WALL_B     ! absorbed IR rad by wall
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_SW_PANEL      ! absorbed solar rad by solar panels
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_LW_PANEL      ! absorbed IR rad by solar panels
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_SW_GARDEN     ! absorbed solar rad by green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_LW_GARDEN     ! absorbed IR rad by green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_SW_GREENROOF  ! absorbed solar rad by green roofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_LW_GREENROOF  ! absorbed IR rad by green roofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PG_GREENROOF_ROOF  ! heat flux between base of greenroof
!                                                         !    and structural roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRUNOFF_GREENROOF  ! greenroof surface runoff
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDRAIN_GREENROOF   ! greenroof total vertical drainage
REAL, DIMENSION(:)  , INTENT(OUT)   :: PIRRIG_GREENROOF   ! water supply for greenroof ground irrigation
!
! new arguments created after BEM
!
 CHARACTER(LEN=6)    , INTENT(IN)    :: HCOOL_COIL
 CHARACTER(LEN=6)    , INTENT(IN)    :: HHEAT_COIL
REAL, DIMENSION(:)  , INTENT(IN)    :: PF_WATER_COND      ! fraction of evaporation for the condensers
 CHARACTER(LEN=4), DIMENSION(:), INTENT(IN) :: HNATVENT
INTEGER             , INTENT(IN)    :: KDAY               ! Simulation day
REAL, DIMENSION(:)  , INTENT(INOUT) :: PAUX_MAX           ! Auxiliar variable for autosize calcs
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_FLOOR           ! Floor layers temperatures [K]
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_MASS            ! Internal mass layers temperatures [K]
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_BLD_COOL        ! Sensible cooling energy demand  
                                                          ! of the building [W m-2(bld)]
REAL, DIMENSION(:)  , INTENT(OUT)   :: PT_BLD_COOL        ! Total cooling energy demand  
                                                          ! of the building [W m-2(bld)]
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_BLD_HEAT        ! Heating energy demand       
                                                          ! of the building [W m-2(bld)]
REAL, DIMENSION(:  ), INTENT(OUT)   :: PLE_BLD_COOL       ! Latent cooling energy demand 
                                                          ! of the building [W m-2(bld)]
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_BLD_HEAT       ! Latent heating energy demand 
                                                          ! of the building [W m-2(bld)]
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_WASTE           ! Sensible waste heat from HVAC system
                                                          ! [W m-2(tot)]
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_WASTE          ! Latent waste heat from HVAC system
                                                          ! [W m-2(tot)]
REAL, DIMENSION(:)  , INTENT(IN)    :: PF_WASTE_CAN       ! fraction of waste heat released into the canyon
REAL, DIMENSION(:)  , INTENT(OUT)   :: PHVAC_COOL         ! Energy consumption of the cooling system
                                                          ! [W m-2(bld)]
REAL, DIMENSION(:)  , INTENT(OUT)   :: PHVAC_HEAT         ! Energy consumption of the heating system
                                                          ! [W m-2(bld)]
REAL, DIMENSION(:)  , INTENT(IN)    :: PQIN               ! Internal heat gains [W m-2(floor)]
REAL, DIMENSION(:)  , INTENT(IN)    :: PQIN_FRAD          ! Radiant fraction of internal heat gains
REAL, DIMENSION(:)  , INTENT(IN)    :: PQIN_FLAT          ! Latent franction of internal heat gains
REAL, DIMENSION(:)  , INTENT(IN)    :: PGR                ! Glazing ratio
REAL, DIMENSION(:)  , INTENT(IN)    :: PEFF_HEAT          ! Efficiency of the heating system
REAL, DIMENSION(:)  , INTENT(IN)    :: PINF               ! Infiltration flow rate [AC/H]
REAL, DIMENSION(:)  , INTENT(IN)    :: PTCOOL_TARGET      ! Cooling setpoint of HVAC system [K]
REAL, DIMENSION(:)  , INTENT(IN)    :: PTHEAT_TARGET      ! Heating setpoint of HVAC system [K]
REAL, DIMENSION(:)  , INTENT(IN)    :: PHR_TARGET         ! Relative humidity setpoint
REAL, DIMENSION(:)  , INTENT(INOUT) :: PT_WIN2            ! Indoor window temperature [K]
REAL, DIMENSION(:)  , INTENT(INOUT) :: PQI_BLD            ! Indoor air specific humidity [kg kg-1]
REAL, DIMENSION(:)  , INTENT(IN)    :: PV_VENT            ! Ventilation flow rate [AC/H]
REAL, DIMENSION(:)  , INTENT(IN)    :: PCAP_SYS_HEAT      ! Capacity of the heating system 
                                                          ! [W m-2(bld)]
REAL, DIMENSION(:)  , INTENT(INOUT) :: PCAP_SYS_RAT       ! Rated capacity of the cooling system
                                                          ! [W m-2(bld)]
REAL, DIMENSION(:)  , INTENT(IN)    :: PT_ADP             ! Apparatus dewpoint temperature of the
                                                          ! cooling coil [K]
REAL, DIMENSION(:)  , INTENT(INOUT) :: PM_SYS_RAT         ! Rated HVAC mass flow rate 
                                                          ! [kg s-1 m-2(bld)]
REAL, DIMENSION(:)  , INTENT(IN)    :: PCOP_RAT           ! Rated COP of the cooling system
REAL, DIMENSION(:)  , INTENT(OUT)   :: PCAP_SYS           ! Actual capacity of the cooling system
                                                          ! [W m-2(bld)] 
REAL, DIMENSION(:)  , INTENT(OUT)   :: PM_SYS             ! Actual HVAC mass flow rate 
                                                          ! [kg s-1 m-2(bld)]
REAL, DIMENSION(:)  , INTENT(OUT)   :: PCOP               ! COP of the cooling system
REAL, DIMENSION(:)  , INTENT(OUT)   :: PQ_SYS             ! Supply air specific humidity [kg kg-1]
REAL, DIMENSION(:)  , INTENT(OUT)   :: PT_SYS             ! Supply air temperature [K]
REAL, DIMENSION(:)  , INTENT(OUT)   :: PTR_SW_WIN         ! Solar radiation transmitted throught
                                                          ! windows [W m-2(bld)]
REAL, DIMENSION(:)  , INTENT(OUT)   :: PFAN_POWER         ! HVAC fan power
REAL, DIMENSION(:,:), INTENT(IN)    :: PHC_FLOOR          ! heat capacity for road layers
REAL, DIMENSION(:,:), INTENT(IN)    :: PTC_FLOOR          ! thermal conductivity for 
                                                          ! road layers
REAL, DIMENSION(:,:), INTENT(IN)    :: PD_FLOOR           ! depth of road layers

REAL, DIMENSION(:)  , INTENT(INOUT) :: PT_WIN1            ! outdoor window temperature [K]
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_SW_WIN        ! window absorbed shortwave radiation [W m-2] 
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_LW_WIN        ! absorbed infrared rad. [W m-2]
REAL, DIMENSION(:)  , INTENT(IN)    :: PSHGC              ! window solar transmittance
REAL, DIMENSION(:)  , INTENT(IN)    :: PSHGC_SH           ! window + shading solar heat gain coef.
REAL, DIMENSION(:)  , INTENT(IN)    :: PUGG_WIN           ! window glass-to-glass U-factro [W m-2 K-1]
REAL, DIMENSION(:)  , INTENT(OUT)   :: PALB_WIN           ! window albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PABS_WIN           ! window absortance
REAL, DIMENSION(:)  , INTENT(IN)    :: PTRAN_WIN          ! window transmittance
!
! new argument for the UTCI calculation
REAL, DIMENSION(:)  , INTENT(OUT)    :: PEMIT_LW_GRND     ! LW flux emitted by the ground (W/m2 ground)
REAL, DIMENSION(:)  , INTENT(OUT)    :: PEMIT_LW_FAC      ! LW flux emitted by the facade (W/m2 ground)
REAL, DIMENSION(:)  , INTENT(OUT)    :: PT_RAD_IND        ! Indoor mean radiant temperature [K]
REAL, DIMENSION(:)  , INTENT(OUT)    :: PREF_SW_GRND      ! total solar rad reflected from ground
REAL, DIMENSION(:)  , INTENT(OUT)    :: PREF_SW_FAC       ! total solar rad reflected from facade
REAL, DIMENSION(:)  , INTENT(OUT)    :: PHU_BLD           ! Indoor relative humidity 0 < (-) < 1
!
! new arguments for shading, schedule or natural ventilation
REAL                , INTENT(IN)     :: PTIME             ! current time since midnight (UTC, s)
LOGICAL,DIMENSION(:), INTENT(IN)     :: OSHADE
LOGICAL,DIMENSION(:), INTENT(INOUT)  :: OSHAD_DAY         ! has shading been necessary this day ?
LOGICAL,DIMENSION(:), INTENT(INOUT)  :: ONATVENT_NIGHT    ! has natural ventilation been necessary/possible this night ?
 CHARACTER(LEN=3)    , INTENT(IN)     :: HBEM              ! Building Energy model 'DEF' or 'BEM'
!
REAL, DIMENSION(:)  , INTENT(IN)     :: PN_FLOOR          ! Number of floors     
REAL, DIMENSION(:)  , INTENT(IN)     :: PWALL_O_BLD       ! Wall area [m2_wall/m2_bld]
REAL, DIMENSION(:)  , INTENT(IN)     :: PGLAZ_O_BLD       ! Window area [m2_win/m2_bld]
REAL, DIMENSION(:)  , INTENT(IN)     :: PMASS_O_BLD       ! Mass area [m2_mass/m2_bld]
REAL, DIMENSION(:)  , INTENT(IN)     :: PFLOOR_HW_RATIO   ! H/W ratio of 1 floor level
REAL, DIMENSION(:)  , INTENT(IN)     :: PF_FLOOR_MASS     ! View factor floor-mass
REAL, DIMENSION(:)  , INTENT(IN)     :: PF_FLOOR_WALL     ! View factor floor-wall
REAL, DIMENSION(:)  , INTENT(IN)     :: PF_FLOOR_WIN      ! View factor floor-window
REAL, DIMENSION(:)  , INTENT(IN)     :: PF_FLOOR_ROOF     ! View factor floor-roof
REAL, DIMENSION(:)  , INTENT(IN)     :: PF_WALL_FLOOR     ! View factor wall-floor
REAL, DIMENSION(:)  , INTENT(IN)     :: PF_WALL_MASS      ! View factor wall-mass
REAL, DIMENSION(:)  , INTENT(IN)     :: PF_WALL_WIN       ! View factor wall-win
REAL, DIMENSION(:)  , INTENT(IN)     :: PF_WIN_FLOOR      ! View factor win-floor
REAL, DIMENSION(:)  , INTENT(IN)     :: PF_WIN_MASS       ! View factor win-mass
REAL, DIMENSION(:)  , INTENT(IN)     :: PF_WIN_WALL       ! View factor win-wall
REAL, DIMENSION(:)  , INTENT(IN)     :: PF_MASS_FLOOR     ! View factor mass-floor
REAL, DIMENSION(:)  , INTENT(IN)     :: PF_MASS_WALL      ! View factor mass-wall
REAL, DIMENSION(:)  , INTENT(IN)     :: PF_MASS_WIN       ! View factor mass-window
LOGICAL             , INTENT(IN)     :: OCANOPY           ! is canopy active ?
 CHARACTER(LEN=5)    , INTENT(IN)     :: HCH_BEM         ! TEB option for building outside conv. coef
REAL, DIMENSION(:)  , INTENT(IN)     :: PROUGH_ROOF     ! roof roughness coef.
REAL, DIMENSION(:)  , INTENT(IN)     :: PROUGH_WALL     ! wall roughness coef.
REAL, DIMENSION(:)  , INTENT(IN)     :: PF_WIN_WIN      ! indoor win to win view factor
!
!* arguments for road watering
!
LOGICAL,              INTENT(IN)      :: OPAR_RD_IRRIG  ! True if there is potentially irrigation
REAL, DIMENSION(:),   INTENT(IN)      :: PRD_START_MONTH! First month of irrigation
REAL, DIMENSION(:),   INTENT(IN)      :: PRD_END_MONTH  ! Last month of irrigation
REAL, DIMENSION(:),   INTENT(IN)      :: PRD_START_HOUR ! First hour of irrigation (included)
REAL, DIMENSION(:),   INTENT(IN)      :: PRD_END_HOUR   ! Last hour of irrigation  (excluded)
REAL, DIMENSION(:),   INTENT(IN)      :: PRD_24H_IRRIG  ! diurnal averaged irrigation flux (kg/m2/s)
REAL, DIMENSION(:),   INTENT(OUT)     :: PIRRIG_ROAD    ! man-made watering of roads (kg/m2/s)
!
!* arguments for solar panels
!
REAL, DIMENSION(:),   INTENT(IN)      :: PEMIS_PANEL    ! Emissivity of solar panels (-)
REAL, DIMENSION(:),   INTENT(IN)      :: PALB_PANEL     ! Albedo     of solar panels (-)
REAL, DIMENSION(:),   INTENT(IN)      :: PEFF_PANEL     ! Efficiency of solar panel (-)
REAL, DIMENSION(:),   INTENT(IN)      :: PFRAC_PANEL    ! Fraction of solar panel on roofs (-)
REAL, DIMENSION(:),   INTENT(OUT)     :: PTHER_PROD_PANEL ! Thermal      Energy production of solar panel on roofs (W/m2 panel)
REAL, DIMENSION(:),   INTENT(OUT)     :: PPHOT_PROD_PANEL ! Photovoltaic Energy production of solar panel on roofs (W/m2 panel)
REAL, DIMENSION(:),   INTENT(OUT)     :: PPROD_PANEL      ! Averaged     Energy production of solar panel on roofs (W/m2 panel)
REAL, DIMENSION(:),   INTENT(OUT)     :: PTHER_PROD_BLD   ! Thermal      Energy production of solar panel on roofs (W/m2 bld  )
REAL, DIMENSION(:),   INTENT(OUT)     :: PPHOT_PROD_BLD   ! Photovoltaic Energy production of solar panel on roofs (W/m2 bld  )
REAL, DIMENSION(:),   INTENT(OUT)     :: PPROD_BLD        ! Averaged     Energy production of solar panel on roofs (W/m2 bld  )
REAL, DIMENSION(:),   INTENT(OUT)     :: PTHER_PRODC_DAY  ! Present day integrated thermal production of energy    (J/m2 panel)
REAL, DIMENSION(:),   INTENT(OUT)     :: PH_PANEL       ! Sensible heat flux from solar panels      (W/m2 panel)
REAL, DIMENSION(:),   INTENT(OUT)     :: PRN_PANEL      ! Net radiation        of solar panel       (W/m2 panel)
!
!* arguments for building's use description
!
REAL, DIMENSION(:),   INTENT(IN)      :: PRESIDENTIAL      ! Fraction of residential use (-)
REAL,                 INTENT(IN)      :: PDT_RES           ! target temperature change when unoccupied (K) (residential buildings)
REAL,                 INTENT(IN)      :: PDT_OFF           ! target temperature change when unoccupied (K) (office buildings)
REAL, DIMENSION(:),   INTENT(OUT)     :: PCUR_TCOOL_TARGET ! Cooling target temperature at current time (K)
REAL, DIMENSION(:),   INTENT(OUT)     :: PCUR_THEAT_TARGET ! Heating target temperature at current time (K)
REAL, DIMENSION(:),   INTENT(OUT)     :: PCUR_QIN          ! Internal heat gains        at current time (W/m2 floor)
!
!
!*      0.2    Declarations of local variables
!
REAL, DIMENSION(SIZE(PTA)) :: ZTA            ! air temperature extrapolated at roof level
REAL, DIMENSION(SIZE(PTA)) :: ZQA            ! air humidity extrapolated at roof level
!
REAL, DIMENSION(SIZE(PTA)) :: ZDN_ROOF       ! snow fraction on roofs
REAL, DIMENSION(SIZE(PTA)) :: ZDN_ROAD       ! snow fraction on roads
REAL, DIMENSION(SIZE(PTA)) :: ZDF_ROOF       ! free-snow fraction on roofs
REAL, DIMENSION(SIZE(PTA)) :: ZDF_ROAD       ! free-snow fraction on roads
REAL, DIMENSION(SIZE(PTA)) :: ZDELT_ROAD     ! fraction of water on roads
REAL, DIMENSION(SIZE(PTA)) :: ZDELT_ROOF     ! fraction of water on roofs
REAL, DIMENSION(SIZE(PTA)) :: ZAC_ROOF       ! roof conductance
REAL, DIMENSION(SIZE(PTA)) :: ZAC_ROOF_WAT   ! roof water conductance
REAL, DIMENSION(SIZE(PTA)) :: ZAC_WALL       ! wall conductance
REAL, DIMENSION(SIZE(PTA)) :: ZAC_TOP        ! top conductance
REAL, DIMENSION(SIZE(PTA)) :: ZQSAT_ROAD     ! hum of saturation for roads
REAL, DIMENSION(SIZE(PTA)) :: ZQSAT_GARDEN   ! hum of saturation for green areas
REAL, DIMENSION(SIZE(PTA)) :: ZQSAT_ROOF     ! hum of saturation for roofs
REAL, DIMENSION(SIZE(PTA)) :: ZQSAT_GREENROOF! hum of saturation for green roofs
!
! coefficients for LW computations over snow (from previous time-step)
!
REAL, DIMENSION(SIZE(PTA)) :: ZTSSNOW_ROOF   ! roof snow temp at previous time-step
REAL, DIMENSION(SIZE(PTA)) :: ZTSSNOW_ROAD   ! road snow temp at previous time-step
REAL, DIMENSION(SIZE(PTA)) :: ZESNOW_ROOF    ! snow emissivity at previous time-step
REAL, DIMENSION(SIZE(PTA)) :: ZESNOW_ROAD    ! snow emissivity at previous time-step
!
! incoming shortwave radiation
!
REAL, DIMENSION(SIZE(PTA)) :: ZDIR_SW             ! direct  solar rad
REAL, DIMENSION(SIZE(PTA)) :: ZSCA_SW             ! diffuse solar rad
INTEGER                    :: JSWB
!
! albedo & emissivity
!
REAL, DIMENSION(SIZE(PTA)) :: ZALB_GARDEN    ! albedo     for green areas
REAL, DIMENSION(SIZE(PTA)) :: ZEMIS_GARDEN   ! emissivity for green areas
REAL, DIMENSION(SIZE(PTA)) :: ZALB_GREENROOF ! albedo     for green roofs
REAL, DIMENSION(SIZE(PTA)) :: ZEMIS_GREENROOF! emissivity for green roofs
!
REAL, DIMENSION(SIZE(PTA)) :: ZALBNIR_TVEG_GARDEN      ! nearIR  veg tot albedo
REAL, DIMENSION(SIZE(PTA)) :: ZALBVIS_TVEG_GARDEN      ! visible veg tot albedo
REAL, DIMENSION(SIZE(PTA)) :: ZALBNIR_TSOIL_GARDEN     ! nearIR  soil tot albedo
REAL, DIMENSION(SIZE(PTA)) :: ZALBVIS_TSOIL_GARDEN     ! visible soil tot albedo
!
REAL, DIMENSION(SIZE(PTA)) :: ZALBNIR_TVEG_GREENROOF      ! nearIR  veg tot albedo
REAL, DIMENSION(SIZE(PTA)) :: ZALBVIS_TVEG_GREENROOF      ! visible veg tot albedo
REAL, DIMENSION(SIZE(PTA)) :: ZALBNIR_TSOIL_GREENROOF     ! nearIR  soil tot albedo
REAL, DIMENSION(SIZE(PTA)) :: ZALBVIS_TSOIL_GREENROOF     ! visible soil tot albedo
!
! radiation received by surfaces
!
REAL, DIMENSION(SIZE(PTA)) :: ZREC_SW_ROAD        ! solar rad received by roads
REAL, DIMENSION(SIZE(PTA)) :: ZREC_SW_WALL_A      ! solar rad received by walls
REAL, DIMENSION(SIZE(PTA)) :: ZREC_SW_WALL_B      ! solar rad received by walls
REAL, DIMENSION(SIZE(PTA)) :: ZREC_SW_GARDEN      ! solar rad received by gardens
REAL, DIMENSION(SIZE(PTA)) :: ZREC_SW_SNOW_ROAD   ! solar rad received by snow on roads
!
REAL, DIMENSION(SIZE(PTA)) :: ZREC_LW_GARDEN      ! IR rad received by gardens
!
REAL, DIMENSION(SIZE(PTA)) :: ZSW_RAD_GARDEN      ! solar radiation reaching urban green areas
REAL, DIMENSION(SIZE(PTA)) :: ZREC_SW_ROOF        ! solar rad received by roofs in presence of solar panels
!
! coefficients for LW contributions
!
REAL, DIMENSION(SIZE(PTA)) :: ZLW_WA_TO_WB        ! LW contrib. wall       -> opposite wall
REAL, DIMENSION(SIZE(PTA)) :: ZLW_WA_TO_R         ! LW contrib. wall       -> road
REAL, DIMENSION(SIZE(PTA)) :: ZLW_WB_TO_R         ! LW contrib. wall       -> road
REAL, DIMENSION(SIZE(PTA)) :: ZLW_WA_TO_G         ! LW contrib. wall       -> green
REAL, DIMENSION(SIZE(PTA)) :: ZLW_WB_TO_G         ! LW contrib. wall       -> green
REAL, DIMENSION(SIZE(PTA)) :: ZLW_WA_TO_WIN       ! Radiative heat trasfer coeff wall-window 
REAL, DIMENSION(SIZE(PTA)) :: ZLW_WB_TO_WIN       ! Radiative heat trasfer coeff wall-window 
                                                  ! [W K-1 m-2] 
REAL, DIMENSION(SIZE(PTA)) :: ZLW_WA_TO_NR        ! LW contrib. wall       -> road(snow)
REAL, DIMENSION(SIZE(PTA)) :: ZLW_WB_TO_NR        ! LW contrib. wall       -> road(snow)
REAL, DIMENSION(SIZE(PTA)) :: ZLW_R_TO_WA         ! LW contrib. road       -> wall
REAL, DIMENSION(SIZE(PTA)) :: ZLW_R_TO_WB         ! LW contrib. road       -> wall
REAL, DIMENSION(SIZE(PTA)) :: ZLW_R_TO_WIN        ! Radiative heat trasfer coeff road-window 
                                                  ! [W K-1 m-2] 
REAL, DIMENSION(SIZE(PTA)) :: ZLW_G_TO_WA         ! LW contrib. green      -> wall
REAL, DIMENSION(SIZE(PTA)) :: ZLW_G_TO_WB         ! LW contrib. green      -> wall
REAL, DIMENSION(SIZE(PTA)) :: ZLW_G_TO_WIN        ! L.W. interactions GARDEN areas->road
                                                  !
REAL, DIMENSION(SIZE(PTA)) :: ZLW_WIN_TO_WA       ! Radiative heat trasfer coeff window-wall
                                                  ! [W K-1 m-2] 
REAL, DIMENSION(SIZE(PTA)) :: ZLW_WIN_TO_WB       ! Radiative heat trasfer coeff window-wall
                                                  ! [W K-1 m-2] 
REAL, DIMENSION(SIZE(PTA)) :: ZLW_WIN_TO_R        ! Radiative heat trasfer coeff window-road 
                                                  ! [W K-1 m-2]
REAL, DIMENSION(SIZE(PTA)) :: ZLW_WIN_TO_G        ! Radiative heat trasfer coeff window-garden 
                                                  ! [W K-1 m-2]
REAL, DIMENSION(SIZE(PTA)) :: ZLW_WIN_TO_NR       ! Radiative heat trasfer coeff window-road(snow) 
                                                  ! [W K-1 m-2] 
REAL, DIMENSION(SIZE(PTA)) :: ZLW_S_TO_WA         ! LW contrib. sky        -> wall
REAL, DIMENSION(SIZE(PTA)) :: ZLW_S_TO_WB         ! LW contrib. sky        -> wall
REAL, DIMENSION(SIZE(PTA)) :: ZLW_S_TO_R          ! LW contrib. sky        -> road
REAL, DIMENSION(SIZE(PTA)) :: ZLW_S_TO_G          ! LW contrib. sky        -> green
REAL, DIMENSION(SIZE(PTA)) :: ZLW_S_TO_WIN        ! Radiative heat trasfer coeff window-sky 
                                                  ! [W K-1 m-2]
REAL, DIMENSION(SIZE(PTA)) :: ZLW_S_TO_NR         ! LW contrib. sky        -> road(snow)
REAL, DIMENSION(SIZE(PTA)) :: ZLW_NR_TO_WA        ! LW contrib. road(snow) -> wall
REAL, DIMENSION(SIZE(PTA)) :: ZLW_NR_TO_WB        ! LW contrib. road(snow) -> wall
REAL, DIMENSION(SIZE(PTA)) :: ZLW_NR_TO_WIN       ! L.W. interactions snow(road)->GARDEN areas
                                                  !
REAL, DIMENSION(SIZE(PTA)) :: ZREC_SW_WIN         ! solar received by windows [W m-2(win)]
REAL, DIMENSION(SIZE(PTA)) :: ZT_SKY              ! sky temperature [K]
REAL, DIMENSION(SIZE(PTA)) :: ZREC_LW_ROOF        ! Incoming LW on roofs in presence of solar panels
!
! local variable at previous time-step
!
REAL, DIMENSION(SIZE(PTA)) :: ZPET_A_COEF          
REAL, DIMENSION(SIZE(PTA)) :: ZPET_B_COEF          
REAL, DIMENSION(SIZE(PTA)) :: ZPEQ_A_COEF          
REAL, DIMENSION(SIZE(PTA)) :: ZPEQ_B_COEF          
!
REAL, DIMENSION(SIZE(PTA)) :: ZUW_ROAD            ! momentum flux for roads
REAL, DIMENSION(SIZE(PTA)) :: ZUW_GARDEN          ! momentum flux for green areas
REAL, DIMENSION(SIZE(PTA)) :: ZUW_GREENROOF       ! momentum flux for green roofs
REAL, DIMENSION(SIZE(PTA)) :: ZDUWDU_ROAD         !
!
REAL, DIMENSION(SIZE(PTA)) :: ZAC_AGG_GARDEN      ! aggreg. aeodynamic resistance for green areas
REAL, DIMENSION(SIZE(PTA)) :: ZHU_AGG_GARDEN      ! aggreg. relative humidity for green areas
REAL, DIMENSION(SIZE(PTA)) :: ZAC_AGG_GREENROOF   ! aggreg. aeodynamic resistance for green roofs
REAL, DIMENSION(SIZE(PTA)) :: ZHU_AGG_GREENROOF   ! aggreg. relative humidity for green roofs
!
!  surfaces relative fractions
!
REAL, DIMENSION(SIZE(PTA)) :: ZROOF_FRAC          ! roof, wall and
REAL, DIMENSION(SIZE(PTA)) :: ZWALL_FRAC          ! road fractions
REAL, DIMENSION(SIZE(PTA)) :: ZROAD_FRAC          ! of exchange surf.
REAL, DIMENSION(SIZE(PTA)) :: ZGARDEN_FRAC        !                  
REAL, DIMENSION(SIZE(PTA)) :: ZTOTS_O_HORS        ! total canyon+roof surface
!                                                 ! over horizontal surface
REAL, DIMENSION(SIZE(PTA)) :: ZWALL_O_ROAD        ! wall surface over road surface
REAL, DIMENSION(SIZE(PTA)) :: ZWALL_O_GRND        ! wall surface over (road+green area) surface
!
! surface temperatures
!
REAL, DIMENSION(SIZE(PTA)) :: ZTS_GARDEN          ! surface temperature of urban green areas at t
REAL, DIMENSION(SIZE(PTA)) :: ZTS_GREENROOF       ! surface temperature of urban greenroofs at t
REAL, DIMENSION(SIZE(PTA)) :: ZMTC_O_GR_R1        ! mean thermal conductivity over distance 
!                                                 ! between two layers (bottom GR & roof)
!
! fluxes from green surfaces
!
REAL, DIMENSION(SIZE(PTA)) :: ZEVAP_GARDEN        ! evaporation (kg/m2/s)
REAL, DIMENSION(SIZE(PTA)) :: ZSFCO2_GARDEN       ! CO2 fluxes (m/s*kg_CO2/kg_air)
REAL, DIMENSION(SIZE(PTA)) :: ZEMIT_LW_GARDEN     ! LW flux emitted by the garden (W/m2 garden)
REAL, DIMENSION(SIZE(PTA)) :: ZEVAP_GREENROOF     ! evaporation over greenroofs (kg/m2/s)
REAL, DIMENSION(SIZE(PTA)) :: ZSFCO2_GREENROOF    ! CO2 fluxes over greenroofs (m/s*kg_CO2/kg_air)
!
! fluxes from built surfaces
REAL, DIMENSION(SIZE(PTA)) :: ZEMIT_LW_ROAD       ! LW flux emitted by the road (W/m2 road)
!
! fluxes from/to solar panel
REAL, DIMENSION(SIZE(PTA)) :: ZEMIT_LWDN_PANEL    ! LW flux emitted DOWNWARDS by the solar panel (W/m2 panel)
REAL, DIMENSION(SIZE(PTA)) :: ZEMIT_LWUP_PANEL    ! LW flux emitted UPWARDS   by the solar panel (W/m2 panel)
REAL, DIMENSION(SIZE(PTA)) :: ZEMIT_LW_ROOF       ! LW flux emitted UPWARDS   by the roof        (W/m2 roof )
!
!new local variables for shading
REAL, DIMENSION(SIZE(PTA)) :: ZE_SHADING          ! energy not ref., nor absorbed, nor
                                                  ! trans. by glazing [Wm-2(win)]
LOGICAL, DIMENSION(SIZE(PTA)) :: GSHADE           ! describes if one encounters the
!                                                 ! conditions to close windows
!                                              
!
INTEGER :: JJ

REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*      1.     Initializations
!              ---------------
!
!*      1.0    broadband radiative fluxes
!              --------------------------
!
IF (LHOOK) CALL DR_HOOK('TEB_GARDEN',0,ZHOOK_HANDLE)
ZDIR_SW(:) = 0.
ZSCA_SW(:) = 0.
!
DO JSWB=1,KSW
  DO JJ=1,SIZE(PDIR_SW,1)
    ZDIR_SW(JJ) = ZDIR_SW(JJ) + PDIR_SW(JJ,JSWB)
    ZSCA_SW(JJ) = ZSCA_SW(JJ) + PSCA_SW(JJ,JSWB)
  ENDDO
END DO
!
!
!*      1.1    surfaces relative fractions
!              ---------------------------
!
DO JJ=1,SIZE(PROAD)
  IF (PROAD(JJ) .GT. 0.) THEN
    ZTOTS_O_HORS(JJ) = 1. + PWALL_O_HOR(JJ)
    ZROOF_FRAC(JJ)   = PBLD(JJ)        / ZTOTS_O_HORS(JJ)
    ZWALL_FRAC(JJ)   = PWALL_O_HOR(JJ) / ZTOTS_O_HORS(JJ)
    ZROAD_FRAC(JJ)   = PROAD(JJ)       / ZTOTS_O_HORS(JJ)
    ZGARDEN_FRAC(JJ) = PGARDEN(JJ)     / ZTOTS_O_HORS(JJ)
    ZWALL_O_ROAD(JJ) = ZWALL_FRAC(JJ) /  ZROAD_FRAC(JJ)
    ZWALL_O_GRND(JJ) = ZWALL_FRAC(JJ) / (ZROAD_FRAC(JJ)+ZGARDEN_FRAC(JJ))
  ELSEIF (PBLD(JJ) .EQ. 0. .AND. PROAD(JJ) .EQ. 0.) THEN
    ZTOTS_O_HORS(JJ) = 1. + PWALL_O_HOR(JJ)
    ZROOF_FRAC(JJ)   = 0.
    ZWALL_FRAC(JJ)   = PWALL_O_HOR(JJ) / ZTOTS_O_HORS(JJ)
    ZROAD_FRAC(JJ)   = 0.
    ZGARDEN_FRAC(JJ) = PGARDEN(JJ)     / ZTOTS_O_HORS(JJ)
    ZWALL_O_ROAD(JJ) = 0.
    ZWALL_O_GRND(JJ) = ZWALL_FRAC(JJ) / ZGARDEN_FRAC(JJ)
  ENDIF
ENDDO
!
!-------------------------------------------------------------------------------
!
!*      2.     Snow-covered surfaces relative effects
!              --------------------------------------
!
!*      2.1    Snow-covered surfaces relative fractions (at previous time-step)
!              ----------------------------------------
 CALL SNOW_FRAC_ROAD(PWSNOW_ROAD(:,1),PSR(:)>0.,ZDN_ROAD,ZDF_ROAD)
 CALL SNOW_FRAC_ROOF(PWSNOW_ROOF(:,1),PSR(:)>0.,ZDN_ROOF,ZDF_ROOF)
!
!* new snow albedo
!
WHERE (PWSNOW_ROAD(:,1)==0. .AND. PSR(:)>0.) PASNOW_ROAD(:) = XANSMAX
WHERE (PWSNOW_ROOF(:,1)==0. .AND. PSR(:)>0.) PASNOW_ROOF(:) = XANSMAX
!
!*      2.2    If snow was not present at previous time-step but is falling
!              ------------------------------------------------------------
!
WHERE (PWSNOW_ROAD(:,1)==0. .AND. PSR(:)>0.)
  PASNOW_ROAD(:) = XANSMAX
  PESNOW_ROAD(:) = XEMISSN
  PTSSNOW_ROAD(:)= MIN(PT_ROAD(:,1), XTT)
END WHERE
WHERE (PWSNOW_ROOF(:,1)==0. .AND. PSR(:)>0.)
  PASNOW_ROOF(:) = XANSMAX
  PESNOW_ROOF(:) = XEMISSN
  PTSSNOW_ROOF(:)= MIN(PT_ROOF(:,1), XTT)
END WHERE
!
!*      2.3    Radiative snow variables at previous time-step
!              ----------------------------------------------
!
ZESNOW_ROOF  (:) = PESNOW_ROOF  (:)
ZESNOW_ROAD  (:) = PESNOW_ROAD  (:)
ZTSSNOW_ROOF (:) = PTSSNOW_ROOF (:)
ZTSSNOW_ROAD (:) = PTSSNOW_ROAD (:)
!
!-------------------------------------------------------------------------------
!
!*      3.     Extrapolation of atmospheric T and q at roof level (for fluxes computation)
!              --------------------------------------------------
!
ZTA(:) = PTA(:) * PEXNS(:) / PEXNA(:)
ZQA(:) = PQA(:) * QSAT(PTA(:),PPS(:)) / QSAT(ZTA(:),PPA(:))
!
!-------------------------------------------------------------------------------
!
!*      4.     Grid-averaged albedo and emissivity of green areas
!              --------------------------------------------------
!
ZALB_GARDEN   = XUNDEF
ZEMIS_GARDEN  = XUNDEF
ZTS_GARDEN    = XUNDEF
!
IF (OGARDEN) THEN
 CALL GARDEN_PROPERTIES(PDIR_SW, PSCA_SW, PSW_BANDS, KSW,                    &
                        ZTS_GARDEN, ZEMIS_GARDEN, ZALB_GARDEN, PTA=PT_LOWCAN )
ENDIF
!
! for greenroofs :
!
ZALB_GREENROOF   = XUNDEF
ZEMIS_GREENROOF  = XUNDEF
ZTS_GREENROOF    = XUNDEF
!
IF (OGREENROOF) THEN
 CALL GREENROOF_PROPERTIES(PDIR_SW, PSCA_SW, PSW_BANDS, KSW,                       &
                           ZTS_GREENROOF, ZEMIS_GREENROOF, ZALB_GREENROOF, PTA=PTA )
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      5.     Occupation of buildings
!              -----------------------
!
!* when building in unoccupied, target temperature is modified
!
IF (HBEM=="BEM") THEN
  CALL BLD_OCC_CALENDAR(TPTIME%TDATE%YEAR,TPTIME%TDATE%MONTH,TPTIME%TDATE%DAY,PTSUN,      &
                        PRESIDENTIAL,PTCOOL_TARGET, PTHEAT_TARGET, PQIN,                  &
                        PDT_RES,PDT_OFF,1.,PCUR_TCOOL_TARGET, PCUR_THEAT_TARGET, PCUR_QIN )
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      5.     Solar radiation
!              ---------------
!
!
!*       5.1   Windows management
!              ------------------
!
!* checks if one encouters the conditions for closing the windows at day
!  (needs to be equiped of solar protections and not to be confortable inside)
!  
GSHADE(:) = .FALSE.
IF (HBEM=="BEM") THEN
 CALL WINDOW_SHADING_AVAILABILITY(OSHADE, PTI_BLD, PCUR_TCOOL_TARGET,GSHADE)
!
!* Reinitialize shading of windows when changing day
!
  WHERE (PTSUN .LT. PTSTEP + 1E-3) OSHAD_DAY(:) = .FALSE.
END IF
!
!
!*       5.2   computes solar radiation exchanges
!              ----------------------------------
!
 CALL URBAN_SOLAR_ABS(HBEM, HROAD_DIR, HWALL_OPT,                  &
                     ZDIR_SW, ZSCA_SW, PZENITH, PAZIM,             &
                     PBLD, PGARDEN, PROAD_DIR, PROAD, PFRAC_GR,    &
                     PWALL_O_HOR, PCAN_HW_RATIO,                   &
                     PALB_ROOF,                                    &
                     PALB_ROAD, PSVF_ROAD, PALB_WALL, PSVF_WALL,   &
                     PFRAC_PANEL, PALB_PANEL,                      &
                     ZALB_GARDEN, PSVF_GARDEN,                     &
                     ZALB_GREENROOF,                               &
                     PASNOW_ROOF, PASNOW_ROAD,                     &
                     ZDN_ROOF, ZDF_ROOF, ZDN_ROAD, ZDF_ROAD,       &
                     PGR, PABS_WIN, PSHGC, PSHGC_SH, PALB_WIN,     &
                     PABS_SW_ROOF, PABS_SW_ROAD,                   &
                     PABS_SW_WALL_A, PABS_SW_WALL_B,               &
                     PABS_SW_GARDEN, PABS_SW_GREENROOF,            &
                     PABS_SW_SNOW_ROOF, PABS_SW_SNOW_ROAD,         &
                     PABS_SW_PANEL,                                &
                     ZREC_SW_ROAD,  ZREC_SW_SNOW_ROAD,             &
                     ZREC_SW_WALL_A, ZREC_SW_WALL_B,               &
                     ZREC_SW_GARDEN, ZREC_SW_ROOF,                 &
                     PDIR_ALB_TOWN,PSCA_ALB_TOWN,                  &
                     ZSW_RAD_GARDEN, PABS_SW_WIN, ZREC_SW_WIN,     &
                     PTRAN_WIN,                                    &
                     PREF_SW_GRND, PREF_SW_FAC,                    &
                     PTR_SW_WIN, ZE_SHADING, OSHAD_DAY,            &
                     GSHADE                                        )
!
!-------------------------------------------------------------------------------
!
!*      6.     LW properties
!              -------------
!
  CALL URBAN_LW_COEF(PGR, PBLD, PLW_RAD,                                &
                     PEMIS_ROAD, PSVF_ROAD, PEMIS_WALL, PSVF_WALL,      &
                     ZEMIS_GARDEN, PROAD, PGARDEN,                      &
                     PESNOW_ROAD,                                       &
                     PTSSNOW_ROAD, PT_WALL_A(:,1), PT_WALL_B(:,1),      &
                     PT_ROAD(:,1), ZTS_GARDEN, PT_WIN1,                 &
                     ZLW_WA_TO_WB, ZLW_WA_TO_R, ZLW_WB_TO_R,            &
                     ZLW_WA_TO_NR, ZLW_WB_TO_NR,ZLW_WA_TO_G,ZLW_WB_TO_G,&
                     ZLW_WA_TO_WIN, ZLW_WB_TO_WIN,                      &
                     ZLW_R_TO_WA, ZLW_R_TO_WB, ZLW_R_TO_WIN,            &
                     ZLW_G_TO_WA, ZLW_G_TO_WB, ZLW_G_TO_WIN,            &
                     ZLW_S_TO_WA, ZLW_S_TO_WB, ZLW_S_TO_R,              &
                     ZLW_S_TO_NR, ZLW_S_TO_G,ZLW_S_TO_WIN,              &
                     ZLW_WIN_TO_WA, ZLW_WIN_TO_WB,                      &
                     ZLW_WIN_TO_R, ZLW_WIN_TO_NR, ZLW_WIN_TO_G,         &
                     ZLW_NR_TO_WA, ZLW_NR_TO_WB, ZLW_NR_TO_WIN          )
!
!-------------------------------------------------------------------------------
!
! The subroutine is splitted in 2 because of compilation optimization issues
 CALL TEB_GARDEN2
 CALL TEB_GARDEN3
!
IF (LHOOK) CALL DR_HOOK('TEB_GARDEN',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
 CONTAINS
!-------------------------------------------------------------------------------
SUBROUTINE TEB_GARDEN2
!
!*      7.     Terms of radiation absorption
!              -----------------------------
!
!*      7.1    IR rad received by gardens (snow free and snow covered separately)
!              --------------------------
!
IF (OGARDEN) THEN
  ZT_SKY  (:) = (PLW_RAD(:)/XSTEFAN)**0.25
  ZREC_LW_GARDEN(:) = (ZLW_S_TO_G  (:) * (ZT_SKY(:)     - ZTS_GARDEN(:))   &
                     + ZLW_WA_TO_G (:) * (PT_WALL_A(:,1)- ZTS_GARDEN(:))   &
                     + ZLW_WB_TO_G (:) * (PT_WALL_B(:,1)- ZTS_GARDEN(:))   &
                     + ZLW_WIN_TO_G(:) * (PT_WIN1(:)    - ZTS_GARDEN(:)) ) &
                     / ZEMIS_GARDEN(:) + XSTEFAN * ZTS_GARDEN(:)**4
ELSE
  ZREC_LW_GARDEN      (:) = XUNDEF
END IF
!
!*     7.2     Effect of solar panels on incoming LW on roofs
!              ----------------------------------------------
!
IF (OSOLAR_PANEL) THEN
  !
  ! solar panels downwards surface is supposed to be at air temperature
  ! and to be have an emissivity of 1.
  !
  ZEMIT_LWDN_PANEL = XSTEFAN * PTA**4
  !
  ! note that, for the time being, one considers that the solar panel 
  ! intercept radiation both above roof and greenroofs (if any)
  ZREC_LW_ROOF = (1.-PFRAC_PANEL(:)) * PLW_RAD            &
                +    PFRAC_PANEL(:)  * ZEMIT_LWDN_PANEL
ELSE
  ZEMIT_LWDN_PANEL = XUNDEF
  ZREC_LW_ROOF     = PLW_RAD
END IF
!
!-------------------------------------------------------------------------------
!
!*      8.     Treatment of green areas
!              ------------------------
!
!*      8.1    Implicit coeefs for T and Q
!              ---------------------------
!
!* explicit coupling for the time being.
!  canopy may need implicitation if there is a lot a garden in the grid mesh
!
ZPET_A_COEF(:) = 0.
ZPET_B_COEF(:) = PT_LOWCAN(:) / PEXNS(:)
ZPEQ_A_COEF(:) = 0.
ZPEQ_B_COEF(:) = PQ_LOWCAN(:)
!
!*      8.2    Call ISBA for green areas
!              -------------------------
!
IF (OGARDEN) THEN
!
  CALL GARDEN(HIMPLICIT_WIND, TPTIME, PTSUN, PPEW_A_COEF_LOWCAN, PPEW_B_COEF_LOWCAN, &
              ZPET_A_COEF, ZPEQ_A_COEF, ZPET_B_COEF, ZPEQ_B_COEF,                  &
              PTSTEP, PZ_LOWCAN,                                                   &
              PT_LOWCAN, PQ_LOWCAN, PEXNS, PRHOA, PCO2, PPS, PRR, PSR, PZENITH,    &
              ZREC_SW_GARDEN, ZREC_LW_GARDEN, PU_LOWCAN,                           &
              PRN_GARDEN,PH_GARDEN,PLE_GARDEN,PGFLUX_GARDEN, ZSFCO2_GARDEN,        &
              ZEVAP_GARDEN, ZUW_GARDEN,PRUNOFF_GARDEN,                             &
              PAC_GARDEN,ZQSAT_GARDEN,ZTS_GARDEN,                                  &
              ZAC_AGG_GARDEN, ZHU_AGG_GARDEN,                                      &
              PDRAIN_GARDEN, PIRRIG_GARDEN                                         )  

  PAC_GARDEN_WAT(:) = PAC_GARDEN(:)
  PABS_SW_GARDEN(:) = (1.-ZALB_GARDEN(:)) * ZREC_SW_GARDEN
  PABS_LW_GARDEN(:) = ZEMIS_GARDEN(:) * ZREC_LW_GARDEN(:) - &
                      XSTEFAN * ZEMIS_GARDEN(:) * ZTS_GARDEN(:)**4 
  ZEMIT_LW_GARDEN(:) = XSTEFAN * ZTS_GARDEN(:)**4 + &
                      (1 - ZEMIS_GARDEN(:)) / ZEMIS_GARDEN(:) * PABS_LW_GARDEN(:)   

ELSE
!
 PRN_GARDEN     (:) = 0.
 PH_GARDEN      (:) = 0.
 PLE_GARDEN     (:) = 0.
 PGFLUX_GARDEN  (:) = 0.
 ZUW_GARDEN     (:) = 0.
 PAC_GARDEN     (:) = 0.
 PGFLUX_GARDEN  (:) = 0.
 ZEVAP_GARDEN   (:) = 0.
 ZSFCO2_GARDEN  (:) = 0.
 ZQSAT_GARDEN   (:) = XUNDEF
 ZTS_GARDEN     (:) = XUNDEF
 ZAC_AGG_GARDEN (:) = XUNDEF
 ZHU_AGG_GARDEN (:) = XUNDEF
 PAC_GARDEN_WAT (:) = XUNDEF
 PABS_SW_GARDEN (:) = XUNDEF
 PABS_LW_GARDEN (:) = XUNDEF
 ZEMIT_LW_GARDEN(:) = 0.
 PRUNOFF_GARDEN (:) = 0.
 PDRAIN_GARDEN  (:) = 0.
 PIRRIG_GARDEN  (:) = 0.
!
ENDIF
!
!*      8.3    Call ISBA for greenroofs
!              -------------------------
!
IF (OGREENROOF) THEN
!
  !
  CALL GREENROOF(HIMPLICIT_WIND, TPTIME, PTSUN, PPEW_A_COEF, PPEW_B_COEF,            &
                ZPET_A_COEF, ZPEQ_A_COEF, ZPET_B_COEF, ZPEQ_B_COEF,                  &
                PTSTEP, PZREF, PUREF,                                                &
                PTA, PQA, PEXNS, PEXNA,PRHOA, PCO2, PPS, PRR, PSR, PZENITH,          &
                ZREC_SW_ROOF, ZREC_LW_ROOF, PVMOD,                                   &
                PRN_GREENROOF,PH_GREENROOF,PLE_GREENROOF,PGFLUX_GREENROOF,           &
                ZSFCO2_GREENROOF,ZEVAP_GREENROOF, ZUW_GREENROOF,                     &
                PAC_GREENROOF,ZQSAT_GREENROOF,ZTS_GREENROOF,                         &
                ZAC_AGG_GREENROOF, ZHU_AGG_GREENROOF,PG_GREENROOF_ROOF,              &
                PRUNOFF_GREENROOF, PDRAIN_GREENROOF, PIRRIG_GREENROOF                )  
!
  PAC_GREENROOF_WAT(:) = PAC_GREENROOF(:)
  PABS_SW_GREENROOF(:) = (1.-ZALB_GREENROOF(:)) * ZREC_SW_ROOF
  PABS_LW_GREENROOF(:) = ZEMIS_GREENROOF * ZREC_LW_ROOF - XSTEFAN * ZEMIS_GREENROOF * ZTS_GREENROOF**4

ELSE
!
 PRN_GREENROOF    (:) = 0.
 PH_GREENROOF     (:) = 0.
 PLE_GREENROOF    (:) = 0.
 PGFLUX_GREENROOF (:) = 0.
 ZUW_GREENROOF    (:) = 0.
 PAC_GREENROOF    (:) = 0.
 ZEVAP_GREENROOF  (:) = 0.
 ZSFCO2_GREENROOF (:) = 0.
 PRUNOFF_GREENROOF(:) = 0.
 PDRAIN_GREENROOF (:) = 0.
 PIRRIG_GREENROOF (:) = 0.
 ZQSAT_GREENROOF  (:) = XUNDEF
 ZTS_GREENROOF    (:) = XUNDEF
 ZAC_AGG_GREENROOF(:) = XUNDEF
 ZHU_AGG_GREENROOF(:) = XUNDEF
 PABS_SW_GREENROOF(:) = XUNDEF
 PABS_LW_GREENROOF(:) = XUNDEF
 ZMTC_O_GR_R1     (:) = XUNDEF
 PG_GREENROOF_ROOF(:) = XUNDEF
!
ENDIF

END SUBROUTINE TEB_GARDEN2
!-------------------------------------------------------------------------------
SUBROUTINE TEB_GARDEN3
!
!*     9.      Treatment of built covers
!              -------------------------
!
  CALL TEB  (HZ0H, HIMPLICIT_WIND, HWALL_OPT, HBEM, TPTIME, PTSUN,    &
             PT_CANYON, PQ_CANYON, PU_CANYON,                         &
             PT_LOWCAN, PQ_LOWCAN, PU_LOWCAN, PZ_LOWCAN,              &
             PTI_BLD,                                                 &
             PT_ROOF, PT_ROAD, PT_WALL_A, PT_WALL_B,                  &
             PWS_ROOF, PWS_ROAD,                                      &
             HSNOW_ROOF,                                              &
             PWSNOW_ROOF, PTSNOW_ROOF, PRSNOW_ROOF, PASNOW_ROOF,      &
             PTSSNOW_ROOF, PESNOW_ROOF,                               &
             HSNOW_ROAD,                                              &
             PWSNOW_ROAD, PTSNOW_ROAD, PRSNOW_ROAD, PASNOW_ROAD,      &
             PTSSNOW_ROAD, PESNOW_ROAD,                               &
             PPEW_A_COEF, PPEW_B_COEF,                                &
             PPEW_A_COEF_LOWCAN, PPEW_B_COEF_LOWCAN,                  &
             PPS, PPA, PEXNS, PEXNA,                                  &
             PTA, PQA, PRHOA,                                         &
             PLW_RAD,                                                 &
             PRR, PSR,                                                &
             PZREF, PUREF, PVMOD,                                     &
             PH_TRAFFIC, PLE_TRAFFIC, PH_INDUSTRY, PLE_INDUSTRY,      &
             PTSTEP,                                                  &
             PZ0_TOWN,                                                &
             PBLD, PGARDEN, PROAD, PFRAC_GR,                          &
             PBLD_HEIGHT, PWALL_O_HOR, PCAN_HW_RATIO, PWALL_O_GRND,   &
             ZDF_ROOF, ZDN_ROOF, ZDF_ROAD, ZDN_ROAD,                  &
             ZQSAT_ROOF, ZQSAT_ROAD, ZDELT_ROOF, ZDELT_ROAD,          &
             PEMIS_ROOF,                                              &
             PHC_ROOF,PTC_ROOF,PD_ROOF,                               &
             PEMIS_ROAD,                                              &
             PHC_ROAD,PTC_ROAD,PD_ROAD,                               &
             PEMIS_WALL,                                              &
             ZTS_GARDEN,                                              &
             PHC_WALL,PTC_WALL,PD_WALL,                               &
             PRN_ROOF, PH_ROOF, PLE_ROOF, PLEW_ROOF, PGFLUX_ROOF,     &
             PRUNOFF_ROOF,                                            &
             PRN_GREENROOF, PH_GREENROOF, PLE_GREENROOF,              &
             PGFLUX_GREENROOF, ZUW_GREENROOF,                         &
             PRUNOFF_GREENROOF, PDRAIN_GREENROOF,                     &
             PRN_STRLROOF, PH_STRLROOF, PLE_STRLROOF, PGFLUX_STRLROOF,&
             PRUNOFF_STRLROOF,                                        &
             PRN_ROAD, PH_ROAD, PLE_ROAD, PLEW_ROAD, PGFLUX_ROAD,     &
             PRUNOFF_ROAD,                                            &
             PRN_WALL_A, PH_WALL_A, PLE_WALL_A, PGFLUX_WALL_A,        &
             PRN_WALL_B, PH_WALL_B, PLE_WALL_B, PGFLUX_WALL_B,        &
             PRN_BLT,PH_BLT,PLE_BLT,PGFLUX_BLT,                       &
             PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,    &
             PMELT_ROOF,                                              &
             PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,    &
             PMELT_ROAD,                                              &
             PG_GREENROOF_ROOF,                                       &
             ZUW_ROAD, PUW_ROOF, ZDUWDU_ROAD, PDUWDU_ROOF,            &
             PUSTAR_TOWN, PCD, PCDN, PCH_TOWN, PRI_TOWN,              &
             PRESA_TOWN, PDQS_TOWN, PQF_TOWN, PQF_BLD, PFLX_BLD,      &
             ZAC_ROOF, PAC_ROAD, ZAC_WALL, ZAC_TOP, PAC_GARDEN,       &
             ZAC_ROOF_WAT, PAC_ROAD_WAT,                              &
             PABS_SW_ROOF, PABS_LW_ROOF,                              &
             PABS_SW_SNOW_ROOF, PABS_LW_SNOW_ROOF,                    &
             PABS_SW_ROAD, PABS_LW_ROAD,                              &
             PABS_SW_SNOW_ROAD, PABS_LW_SNOW_ROAD,                    &
             PABS_SW_WALL_A, PABS_LW_WALL_A,                          &
             PABS_SW_WALL_B, PABS_LW_WALL_B,                          &
             ZLW_WA_TO_WB,                                            &
             ZLW_WA_TO_R, ZLW_WB_TO_R,                                &
             ZLW_WA_TO_NR, ZLW_WB_TO_NR,                              &
             ZLW_R_TO_WA, ZLW_R_TO_WB,                                &
             ZLW_G_TO_WA, ZLW_G_TO_WB,                                &
             ZLW_S_TO_WA, ZLW_S_TO_WB, ZLW_S_TO_R,                    &
             ZLW_S_TO_NR, ZLW_NR_TO_WA, ZLW_NR_TO_WB,                 &
             ZLW_NR_TO_WIN, ZLW_WA_TO_WIN, ZLW_WB_TO_WIN,             &
             ZLW_G_TO_WIN,                                            &
             ZLW_R_TO_WIN, ZLW_S_TO_WIN, ZLW_WIN_TO_WA, ZLW_WIN_TO_WB,&
             ZLW_WIN_TO_R, ZLW_WIN_TO_NR,                             &
             HNATVENT,                                                &
             HCOOL_COIL, PF_WATER_COND, HHEAT_COIL,                   &
             .FALSE., KDAY, PAUX_MAX, PT_FLOOR,                       & 
             PT_MASS, PH_BLD_COOL, PT_BLD_COOL,                       &    
             PH_BLD_HEAT, PLE_BLD_COOL, PLE_BLD_HEAT,                 &
             PH_WASTE, PLE_WASTE, PF_WASTE_CAN, PHVAC_COOL,           &
             PHVAC_HEAT, PCUR_QIN, PQIN_FRAD, PQIN_FLAT,              &
             PGR, PEFF_HEAT, PINF, PCUR_TCOOL_TARGET,                 &
             PCUR_THEAT_TARGET, PHR_TARGET, PT_WIN2, PQI_BLD,         &
             PV_VENT, PCAP_SYS_HEAT, PCAP_SYS_RAT, PT_ADP,            &
             PM_SYS_RAT, PCOP_RAT, PCAP_SYS, PM_SYS, PCOP,            &
             PQ_SYS, PT_SYS, PTR_SW_WIN, PFAN_POWER, PHC_FLOOR,       &
             PTC_FLOOR, PD_FLOOR, PT_WIN1, PABS_SW_WIN, PABS_LW_WIN,  &
             PUGG_WIN, PEMIT_LW_FAC, ZEMIT_LW_ROAD, PT_RAD_IND,       &
             PHU_BLD, PTIME, ZE_SHADING, ONATVENT_NIGHT,              &
             PN_FLOOR, PWALL_O_BLD, PGLAZ_O_BLD,                      &
             PMASS_O_BLD, PFLOOR_HW_RATIO,                            &
             PF_FLOOR_MASS, PF_FLOOR_WALL, PF_FLOOR_WIN,              &
             PF_FLOOR_ROOF, PF_WALL_FLOOR, PF_WALL_MASS, PF_WALL_WIN, &
             PF_WIN_FLOOR, PF_WIN_MASS, PF_WIN_WALL, PF_MASS_FLOOR,   &
             PF_MASS_WALL, PF_MASS_WIN, OCANOPY,                      &
             HCH_BEM, PROUGH_ROOF, PROUGH_WALL, PF_WIN_WIN,           &
             OPAR_RD_IRRIG, PRD_START_MONTH, PRD_END_MONTH,           &
             PRD_START_HOUR, PRD_END_HOUR, PRD_24H_IRRIG, PIRRIG_ROAD )
!
!-------------------------------------------------------------------------------
!
!*    10.      Treatment of solar panels
!              -------------------------
!
IF (OSOLAR_PANEL) THEN
  ! 
  !* LW radiation coming upwards from roofs
  !
  ZEMIT_LW_ROOF =  ZREC_LW_ROOF                                  &
     - (        PFRAC_GR(:)  *               PABS_LW_GREENROOF(:)&
          + (1.-PFRAC_GR(:)) * ZDF_ROOF(:) * PABS_LW_ROOF(:)     &
          + (1.-PFRAC_GR(:)) * ZDN_ROOF(:) * PABS_LW_SNOW_ROOF(:)&
       )
  !
  ! note that, for the time being, one considers that the solar panel 
  ! intercept radiation both above roof and greenroofs (if any)
  CALL SOLAR_PANEL(PTSTEP, PTSUN, PRESIDENTIAL,                   &
                   ZEMIT_LW_ROOF,ZEMIT_LWDN_PANEL, PLW_RAD,       &
                   PABS_SW_PANEL, PTA, PN_FLOOR, PFRAC_PANEL,     &
                   PEMIS_PANEL, PALB_PANEL, PEFF_PANEL,           &
                   PABS_LW_PANEL, PH_PANEL, PRN_PANEL,            &
                   PTHER_PRODC_DAY,                               &
                   PTHER_PROD_PANEL, PPHOT_PROD_PANEL,PPROD_PANEL,&
                   PTHER_PROD_BLD,   PPHOT_PROD_BLD,  PPROD_BLD   )
ELSE
  PABS_LW_PANEL    = XUNDEF
  PTHER_PROD_PANEL = XUNDEF
  PPHOT_PROD_PANEL = XUNDEF
  PPROD_PANEL      = XUNDEF
  PTHER_PROD_BLD   = XUNDEF
  PPHOT_PROD_BLD   = XUNDEF
  PPROD_BLD        = XUNDEF
  PH_PANEL         = XUNDEF
  PRN_PANEL        = XUNDEF
END IF
!-------------------------------------------------------------------------------
!
!*     11.     Aggregation
!              -----------
!
 CALL AVG_URBAN_FLUXES(PTS_TOWN, PEMIS_TOWN,                                    &
                     PT_CANYON, PQ_CANYON,                                     &
                     PT_LOWCAN, PQ_LOWCAN,                                     &
                     PT_ROOF(:,1),PT_ROAD(:,1),PT_WALL_A(:,1), PT_WALL_B(:,1), &
                     ZTS_GARDEN,                                               &
                     ZTA, ZQA, PRHOA, PPS,                                     &
                     PH_TRAFFIC, PLE_TRAFFIC, PH_INDUSTRY, PLE_INDUSTRY,       &
                     PBLD, PROAD, PGARDEN, PWALL_O_HOR, ZWALL_O_GRND,          &
                     PFRAC_GR,                                                 &
                     PEMIS_ROOF, ZESNOW_ROOF, ZEMIS_GREENROOF,                 &
                     PLW_RAD,                                                  &
                     PABS_LW_ROOF, PABS_LW_WALL_A, PABS_LW_WALL_B,             &
                     PABS_LW_ROAD, PABS_LW_GARDEN, PABS_LW_GREENROOF,          &
                     PABS_LW_SNOW_ROOF, PABS_LW_SNOW_ROAD,                     &
                     ZAC_ROOF, ZAC_ROOF_WAT,                                   &
                     ZAC_WALL, PAC_ROAD, PAC_ROAD_WAT, ZAC_TOP,                &
                     PAC_GARDEN,                                               &
                     ZQSAT_GARDEN, ZAC_AGG_GARDEN, ZHU_AGG_GARDEN,             &
                     ZQSAT_ROOF, ZQSAT_ROAD,                                   &
                     ZDELT_ROOF, ZDELT_ROAD,                                   &
                     ZROOF_FRAC, ZWALL_FRAC, ZROAD_FRAC, ZGARDEN_FRAC,         &
                     ZTOTS_O_HORS,                                             &
                     ZDF_ROOF, ZDN_ROOF, ZDF_ROAD, ZDN_ROAD,                   &
                     PRN_ROOF, PH_ROOF, PLE_ROOF, PGFLUX_ROOF,                 &
                     PRN_ROAD, PH_ROAD, PLE_ROAD, PGFLUX_ROAD,                 &
                     PRN_GARDEN, PH_GARDEN, PLE_GARDEN, PGFLUX_GARDEN,         &
                     PRN_WALL_A, PH_WALL_A, PLE_WALL_A, PGFLUX_WALL_A,         &
                     PRN_WALL_B, PH_WALL_B, PLE_WALL_B, PGFLUX_WALL_B,         &
                     PLEW_ROOF, PLESNOW_ROOF,                                  &
                     PLEW_ROAD, PLESNOW_ROAD, PHSNOW_ROAD,                     &
                     ZEVAP_GARDEN, ZEVAP_GREENROOF,                            &
                     PRN_GRND, PH_GRND, PLE_GRND, PGFLUX_GRND,                 &
                     PRN_TOWN, PH_TOWN, PLE_TOWN, PGFLUX_TOWN, PEVAP_TOWN,     &
                     PRUNOFF_GARDEN,PRUNOFF_ROAD,PRUNOFF_STRLROOF,             &
                     PRUNOFF_GREENROOF, PDRAIN_GREENROOF, PRUNOFF_TOWN,        &
                     PABS_LW_PANEL, PEMIS_PANEL, PFRAC_PANEL, PRN_PANEL,       &
                     PH_PANEL,                                                 &
                     PH_WASTE, PLE_WASTE, PF_WASTE_CAN,                        &
                     PABS_LW_WIN, PT_WIN1, PGR, ZEMIT_LW_ROAD, ZEMIT_LW_GARDEN,&
                     PEMIT_LW_GRND, HBEM, PSVF_ROAD, PSVF_GARDEN, PSVF_WALL,   &
                     PGARDEN_O_GRND, PROAD_O_GRND,                             &
                     PEMIS_ROAD, PESNOW_ROAD, PEMIS_WALL, ZEMIS_GARDEN, OCANOPY)
!
PSFCO2(:) = PGARDEN(:) * ZSFCO2_GARDEN(:) + PBLD(:) * PFRAC_GR(:) * ZSFCO2_GREENROOF(:) ! no CO2 flux from built and road yet.
!
!-------------------------------------------------------------------------------
!
!*     12.     Momentum flux for ground built surfaces
!              ---------------------------------------
!
PUW_GRND (:)     = (PROAD(:)*ZUW_ROAD(:) + PGARDEN(:)*ZUW_GARDEN(:)) / (PROAD(:)+PGARDEN(:))
!
PDUWDU_GRND (:)  = 0.
!
END SUBROUTINE TEB_GARDEN3
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE TEB_GARDEN
