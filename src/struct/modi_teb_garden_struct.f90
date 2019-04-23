!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
MODULE MODI_TEB_GARDEN_STRUCT
!
INTERFACE
!
    SUBROUTINE TEB_GARDEN_STRUCT (OGARDEN, OGREENROOF, OSOLAR_PANEL,          &
                     HZ0H, HIMPLICIT_WIND, HROAD_DIR, HWALL_OPT, TPTIME,      &
                     PTSUN, PT_CAN, PQ_CAN, PU_CAN,                           &
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
                     PG_GREENROOF_ROOF,                                       &
                     HCOOL_COIL, PF_WATER_COND, HHEAT_COIL,                   &
                     HNATVENT, PNATVENT, KDAY, PAUX_MAX, PT_FLOOR,            &
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
                     PHU_BLD, PTIME, OSHADE, PSHADE, OSHAD_DAY, ONATVENT_NIGHT, HBEM, &
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
!!****  *TEB_GARDEN_STRUCT*  
!!
!!    PURPOSE
!!    -------
!
! To put all individual variables in the corresponding fortran structure
! for the call to TEB_GARDEN routine.
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
!!      V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!    Original    03/2018
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TYPE_DATE_SURF
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
REAL, DIMENSION(:)  , INTENT(INOUT) :: PT_CAN             ! canyon air temperature
REAL, DIMENSION(:)  , INTENT(INOUT) :: PQ_CAN             ! canyon air specific humidity
REAL, DIMENSION(:)  , INTENT(IN)    :: PU_CAN             ! canyon hor. wind
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
!
! new arguments created after BEM
!
 CHARACTER(LEN=12)   , INTENT(IN)    :: HCOOL_COIL
 CHARACTER(LEN=6)    , INTENT(IN)    :: HHEAT_COIL
REAL, DIMENSION(:)  , INTENT(IN)    :: PF_WATER_COND      ! fraction of evaporation for the condensers
 CHARACTER(LEN=4), DIMENSION(:), INTENT(IN) :: HNATVENT
REAL, DIMENSION(:)  , INTENT(IN)    :: PNATVENT           ! flag to describe surventilation system for i/o 
                                                          ! 0 for NONE, 1 for MANU and 2 for AUTO
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
REAL, DIMENSION(:)  , INTENT(IN)     :: PSHADE            ! flag to activate shading devices -> REAL for i/o 0. or 1     
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
REAL, DIMENSION(:),   INTENT(INOUT)   :: PTHER_PRODC_DAY  ! Present day integrated thermal production of energy    (J/m2 panel)
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
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE TEB_GARDEN_STRUCT
!
END INTERFACE
!
END MODULE MODI_TEB_GARDEN_STRUCT
