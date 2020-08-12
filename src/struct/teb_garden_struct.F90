!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL licence
!SFX_LIC version 2.1. See Licence_CeCILL_V2.1-en.txt and Licence_CeCILL_V2.1-fr.txt  
!SFX_LIC for details.
!     #########
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
                     HNATVENT, PNATVENT, KDAY, PAUX_MAX, PT_FLOOR,                      &
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
                     PCUR_TCOOL_TARGET, PCUR_THEAT_TARGET, PCUR_QIN,          &
                     PDIAGVENT, PFOPEN, PN50, PVENT_BEHAV_ANYWAY,             &
                     PVENT_BEHAV_ADAPTI, PTDESV, PFVSHO,                      &
                     PISMECH, PMECHRATE                                       )
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
USE MODD_SURF_PAR, ONLY : XUNDEF
!
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_TEB_PANEL_n, ONLY : TEB_PANEL_t
USE MODD_TEB_IRRIG_n, ONLY : TEB_IRRIG_t
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_t
!
USE MODI_ALLOC_TEB_STRUCT
USE MODI_DEALLOC_TEB_STRUCT
USE MODI_TEB_GARDEN
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
! New arguments from 8.2  (Ventilation)
REAL, DIMENSION(:),   INTENT(IN) :: PDIAGVENT           ! Is ventilation active?
REAL, DIMENSION(:),   INTENT(IN) :: PFOPEN              ! Fraction of windows to be opened
REAL, DIMENSION(:),   INTENT(IN) :: PN50                ! Airtightness [AC/H at 50
REAL, DIMENSION(:),   INTENT(IN) :: PVENT_BEHAV_ANYWAY  ! Fraction of windows opened in any case
REAL, DIMENSION(:),   INTENT(IN) :: PVENT_BEHAV_ADAPTI  ! Fraction of windows that could be opened
REAL, DIMENSION(:),   INTENT(IN) :: PTDESV              ! Desing temperature for ventilation [K]
REAL, DIMENSION(:),   INTENT(IN) :: PFVSHO              ! Ventilation independent of indoor/outdoor temperature
REAL, DIMENSION(:),   INTENT(IN) :: PISMECH             ! Presence of mechanical ventilation. OFF: 0.
REAL, DIMENSION(:),   INTENT(IN) :: PMECHRATE           ! Air exchange rate due to mechanical ventilation [1/h]
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
TYPE(TEB_OPTIONS_t) :: TOP
TYPE(TEB_t) :: T
TYPE(BEM_OPTIONS_t) :: BOP
TYPE(BEM_t) :: B
TYPE(TEB_PANEL_t) :: TPN
TYPE(TEB_IRRIG_t) :: TIR
TYPE(DIAG_MISC_TEB_t) :: DMT
!
!-------------------------------------------------------------------------------
!
!* All OUT (not IN or INOUT) dummy arguments are set to value XUNDEF at the beginning of the routine.
!
PRN_ROOF         = XUNDEF  ! net radiation over roof
PH_ROOF          = XUNDEF  ! sensible heat flux over roof
PLE_ROOF         = XUNDEF  ! latent heat flux over roof
PLEW_ROOF        = XUNDEF  ! latent heat flux over roof (snow)
PGFLUX_ROOF      = XUNDEF  ! flux through the roof
PRUNOFF_ROOF     = XUNDEF  ! runoff over the ground
PRN_ROAD         = XUNDEF  ! net radiation over road
PH_ROAD          = XUNDEF  ! sensible heat flux over road
PLE_ROAD         = XUNDEF  ! latent heat flux over road
PLEW_ROAD        = XUNDEF  ! latent heat flux over road (snow)
PGFLUX_ROAD      = XUNDEF  ! flux through the road
PRUNOFF_ROAD     = XUNDEF  ! runoff over the ground
PRN_WALL_A       = XUNDEF  ! net radiation over wall
PH_WALL_A        = XUNDEF  ! sensible heat flux over wall
PLE_WALL_A       = XUNDEF  ! latent heat flux over wall
PGFLUX_WALL_A    = XUNDEF  ! flux through the wall
PRN_WALL_B       = XUNDEF  ! net radiation over wall
PH_WALL_B        = XUNDEF  ! sensible heat flux over wall
PLE_WALL_B       = XUNDEF  ! latent heat flux over wall
PGFLUX_WALL_B    = XUNDEF  ! flux through the wall
PRN_STRLROOF     = XUNDEF  ! net radiation over structural roof
PH_STRLROOF      = XUNDEF  ! sensible heat flux over structural roof
PLE_STRLROOF     = XUNDEF  ! latent heat flux over structural roof
PGFLUX_STRLROOF  = XUNDEF  ! flux through the structural roof
PRUNOFF_STRLROOF = XUNDEF  ! runoff aver structural roof
PRN_BLT          = XUNDEF  ! net radiation over built surf 
PH_BLT           = XUNDEF  ! sensible heat flux over built surf 
PLE_BLT          = XUNDEF  ! latent heat flux over built surf 
PGFLUX_BLT       = XUNDEF  ! flux through the built surf 
!
PRNSNOW_ROOF     = XUNDEF  ! net radiation over snow
PHSNOW_ROOF      = XUNDEF  ! sensible heat flux over snow
PLESNOW_ROOF     = XUNDEF  ! latent heat flux over snow
PGSNOW_ROOF      = XUNDEF  ! flux under the snow
PMELT_ROOF       = XUNDEF  ! snow melt
PRNSNOW_ROAD     = XUNDEF  ! net radiation over snow
PHSNOW_ROAD      = XUNDEF  ! sensible heat flux over snow
PLESNOW_ROAD     = XUNDEF  ! latent heat flux over snow
PGSNOW_ROAD      = XUNDEF  ! flux under the snow
PMELT_ROAD       = XUNDEF  ! snow melt
!
PRN_GRND         = XUNDEF  ! net radiation over ground
PH_GRND          = XUNDEF  ! sensible heat flux over ground
PLE_GRND         = XUNDEF  ! latent heat flux over ground
PGFLUX_GRND      = XUNDEF  ! flux through the ground
PRN_TOWN         = XUNDEF  ! net radiation over town
PH_TOWN          = XUNDEF  ! sensible heat flux over town
PLE_TOWN         = XUNDEF  ! latent heat flux over town
PGFLUX_TOWN      = XUNDEF  ! flux through the ground
PEVAP_TOWN       = XUNDEF  ! evaporation flux (kg/m2/s)
PRUNOFF_TOWN     = XUNDEF  ! runoff over the ground
PSFCO2           = XUNDEF  ! flux of CO2       (m/s*kg_CO2/kg_air)
PUW_GRND         = XUNDEF  ! momentum flux for ground built surf
PUW_ROOF         = XUNDEF  ! momentum flux for roofs
PDUWDU_GRND      = XUNDEF  !
PDUWDU_ROOF      = XUNDEF  !
PUSTAR_TOWN      = XUNDEF  ! friciton velocity over town
PCD              = XUNDEF  ! town averaged drag coefficient
PCDN             = XUNDEF  ! town averaged neutral drag coefficient
PCH_TOWN         = XUNDEF  ! town averaged heat transfer coefficient
PRI_TOWN         = XUNDEF  ! town averaged Richardson number
PTS_TOWN         = XUNDEF  ! town surface temperature
PEMIS_TOWN       = XUNDEF  ! town equivalent emissivity
PDIR_ALB_TOWN    = XUNDEF  ! town equivalent direct albedo
PSCA_ALB_TOWN    = XUNDEF  ! town equivalent diffuse albedo
PRESA_TOWN       = XUNDEF  ! town aerodynamical resistance
PDQS_TOWN        = XUNDEF  ! heat storage inside town
PQF_TOWN         = XUNDEF  ! total anthropogenic heat
PQF_BLD          = XUNDEF  ! anthropogenic heat flux of domestic heating
PFLX_BLD         = XUNDEF  ! heat flx from inside bld through its structure
PAC_ROAD         = XUNDEF  ! road conductance
PAC_GARDEN       = XUNDEF  ! green area conductance
PAC_GREENROOF    = XUNDEF  ! green roof conductance
PAC_ROAD_WAT     = XUNDEF  ! road conductance for latent heat
PAC_GARDEN_WAT   = XUNDEF  ! green area conductance for latent heat
PAC_GREENROOF_WAT= XUNDEF  ! green roof conductance for latent heat
!
PABS_SW_ROOF     = XUNDEF  ! absorbed solar rad by roof
PABS_SW_SNOW_ROOF= XUNDEF  ! absorbed solar rad by snow on roof
PABS_LW_ROOF     = XUNDEF  ! absorbed IR rad by roof
PABS_LW_SNOW_ROOF= XUNDEF  ! absorbed IR rad by snow on roof
PABS_SW_ROAD     = XUNDEF  ! absorbed solar rad by road
PABS_SW_SNOW_ROAD= XUNDEF  ! absorbed solar rad by snow on road
PABS_LW_ROAD     = XUNDEF  ! absorbed IR rad by road
PABS_LW_SNOW_ROAD= XUNDEF  ! absorbed IR rad by snow on road
PABS_SW_WALL_A   = XUNDEF  ! absorbed solar rad by wall
PABS_LW_WALL_A   = XUNDEF  ! absorbed IR rad by wall
PABS_SW_WALL_B   = XUNDEF  ! absorbed solar rad by wall
PABS_LW_WALL_B   = XUNDEF  ! absorbed IR rad by wall
PABS_SW_PANEL    = XUNDEF  ! absorbed solar rad by solar panels
PABS_LW_PANEL    = XUNDEF  ! absorbed IR rad by solar panels
PABS_SW_GARDEN   = XUNDEF  ! absorbed solar rad by green areas
PABS_LW_GARDEN   = XUNDEF  ! absorbed IR rad by green areas
PABS_SW_GREENROOF= XUNDEF  ! absorbed solar rad by green roofs
PABS_LW_GREENROOF= XUNDEF  ! absorbed IR rad by green roofs
PG_GREENROOF_ROOF= XUNDEF  ! heat flux between base of greenroof
!                                                         !    and structural roof

PH_BLD_COOL      = XUNDEF  ! Sensible cooling energy demand  
                                                          ! of the building [W m-2(bld)]
PT_BLD_COOL      = XUNDEF  ! Total cooling energy demand  
                                                          ! of the building [W m-2(bld)]
PH_BLD_HEAT      = XUNDEF  ! Heating energy demand       
                                                          ! of the building [W m-2(bld)]
PLE_BLD_COOL     = XUNDEF  ! Latent cooling energy demand 
                                                          ! of the building [W m-2(bld)]
PLE_BLD_HEAT     = XUNDEF  ! Latent heating energy demand 
                                                          ! of the building [W m-2(bld)]
PH_WASTE         = XUNDEF  ! Sensible waste heat from HVAC system
                                                          ! [W m-2(tot)]
PLE_WASTE        = XUNDEF  ! Latent waste heat from HVAC system
                                                          ! [W m-2(tot)]
PHVAC_COOL       = XUNDEF  ! Energy consumption of the cooling system
                                                          ! [W m-2(bld)]
PHVAC_HEAT       = XUNDEF  ! Energy consumption of the heating system
                                                          ! [W m-2(bld)]
PCAP_SYS         = XUNDEF  ! Actual capacity of the cooling system
                                                          ! [W m-2(bld)] 
PM_SYS           = XUNDEF  ! Actual HVAC mass flow rate 
                                                          ! [kg s-1 m-2(bld)]
PCOP             = XUNDEF  ! COP of the cooling system
PQ_SYS           = XUNDEF  ! Supply air specific humidity [kg kg-1]
PT_SYS           = XUNDEF  ! Supply air temperature [K]
PTR_SW_WIN       = XUNDEF  ! Solar radiation transmitted throught
                                                          ! windows [W m-2(bld)]
PFAN_POWER       = XUNDEF  ! HVAC fan power
PABS_SW_WIN      = XUNDEF  ! window absorbed shortwave radiation [W m-2] 
PABS_LW_WIN      = XUNDEF  ! absorbed infrared rad. [W m-2]
PALB_WIN         = XUNDEF  ! window albedo
!
! new argument for the UTCI calculation
PEMIT_LW_GRND    = XUNDEF ! LW flux emitted by the ground (W/m2 ground)
PEMIT_LW_FAC     = XUNDEF ! LW flux emitted by the facade (W/m2 ground)
PT_RAD_IND       = XUNDEF ! Indoor mean radiant temperature [K]
PREF_SW_GRND     = XUNDEF ! total solar rad reflected from ground
PREF_SW_FAC      = XUNDEF ! total solar rad reflected from facade
PHU_BLD          = XUNDEF ! Indoor relative humidity 0 < (-) < 1
PIRRIG_ROAD      = XUNDEF ! man-made watering of roads (kg/m2/s)
!
!* arguments for solar panels
!
PTHER_PROD_PANEL = XUNDEF! Thermal      Energy production of solar panel on roofs (W/m2 panel)
PPHOT_PROD_PANEL = XUNDEF! Photovoltaic Energy production of solar panel on roofs (W/m2 panel)
PPROD_PANEL      = XUNDEF! Averaged     Energy production of solar panel on roofs (W/m2 panel)
PTHER_PROD_BLD   = XUNDEF! Thermal      Energy production of solar panel on roofs (W/m2 bld  )
PPHOT_PROD_BLD   = XUNDEF! Photovoltaic Energy production of solar panel on roofs (W/m2 bld  )
PPROD_BLD        = XUNDEF! Averaged     Energy production of solar panel on roofs (W/m2 bld  )
PH_PANEL         = XUNDEF! Sensible heat flux from solar panels      (W/m2 panel)
PRN_PANEL        = XUNDEF! Net radiation        of solar panel       (W/m2 panel)
!
!* arguments for building's use description
!
PCUR_TCOOL_TARGET = XUNDEF! Cooling target temperature at current time (K)
PCUR_THEAT_TARGET = XUNDEF! Heating target temperature at current time (K)
PCUR_QIN          = XUNDEF! Internal heat gains        at current time (W/m2 floor)
!

!-------------------------------------------------------------------------------
!
!* Allocate all the local structure variables
!
CALL ALLOC_TEB_STRUCT(SIZE(PT_ROOF,2),SIZE(PT_ROAD,2),SIZE(PT_WALL_A,2),SIZE(PT_FLOOR,2), &
                      TOP,T,BOP,B,TPN,TIR,DMT) 
!
! This is of course not optimized here, because at each time step, the structures are allocated, set and deallocated at the end.
!  This is done so only to show that the structures can be independant from the rest of the code above (here the driver.F90 program).
! The structures could be saved from one time-step to another by putting them in a dedicated module for example.
!
!-------------------------------------------------------------------------------
!
! BEM OPTIONS structure
!
BOP%NFLOOR_LAYER = SIZE(PT_FLOOR,2)
BOP%CCOOL_COIL   = HCOOL_COIL
BOP%CHEAT_COIL   = HHEAT_COIL
BOP%LAUTOSIZE    = .FALSE.

!-------------------------------------------------------------------------------
!
! BEM structure
!
!
! Floor parameters
!
B%XHC_FLOOR(1,:) = PHC_FLOOR(1,:)    ! floor layers heat capacity        (J/K/m3)
B%XTC_FLOOR(1,:) = PTC_FLOOR(1,:)    ! floor layers thermal conductivity (W/K/m)
B%XD_FLOOR (1,:) = PD_FLOOR (1,:)    ! depth of floor layers             (m)
!
! HVAC parameters
!
B%XTCOOL_TARGET = PTCOOL_TARGET ! cooling setpoint of indoor air
B%XTHEAT_TARGET = PTHEAT_TARGET ! heating setpoint of indoor air
B%XF_WASTE_CAN = PF_WASTE_CAN   ! fraction of waste heat released into the canyon
B%XEFF_HEAT = PEFF_HEAT         ! efficiency of the heating system
!
! Indoor parameters
!
B%XTI_BLD = PTI_BLD                 ! building interior temperature    (K)
B%XT_FLOOR(1,:) = PT_FLOOR(1,:)     ! floor layer temperatures         (K)
B%XT_MASS (1,:) = PT_MASS (1,:)     ! Air cooled building internal th. mass temperature (K)
!
B%XQIN = PQIN                   ! internal heat gains [W m-2(floor)]
B%XQIN_FRAD = PQIN_FRAD         ! radiant fraction of internal heat gains
B%XSHGC = PSHGC                 ! solar heat gain coef. of windows
B%XSHGC_SH = PSHGC_SH           ! solar heat gain coef. of windows + shading
B%XTRAN_WIN = PTRAN_WIN         ! window transmittance (-)
B%XGR = PGR                     ! glazing ratio
B%XINF = PINF                   ! infiltration/ventilation flow rate [AC/H]
!
! New parameters
!
B%XF_WATER_COND = PF_WATER_COND   ! fraction of evaporation for condensers (cooling system)
B%XAUX_MAX = PAUX_MAX             ! Auxiliar variable for autosize calcs
B%XQIN_FLAT = PQIN_FLAT           ! Latent franction of internal heat gains
B%XHR_TARGET = PHR_TARGET         ! Relative humidity setpoint
B%XT_WIN2 = PT_WIN2               ! Indoor window temperature [K]
B%XQI_BLD = PQI_BLD               ! Indoor air specific humidity [kg kg-1]
B%XV_VENT = PV_VENT               ! Ventilation flow rate [AC/H]
B%XCAP_SYS_HEAT = PCAP_SYS_HEAT   ! Capacity of the heating system 
                                  ! [W m-2(bld)]
B%XCAP_SYS_RAT = PCAP_SYS_RAT     ! Rated capacity of the cooling system
                                  ! [W m-2(bld)]
B%XT_ADP = PT_ADP                 ! Apparatus dewpoint temperature of the
                                  ! cooling coil [K]
B%XM_SYS_RAT = PM_SYS_RAT         ! Rated HVAC mass flow rate 
                                  ! [kg s-1 m-2(bld)]
B%XCOP_RAT = PCOP_RAT             ! Rated COP of the cooling system
B%XT_WIN1 = PT_WIN1               ! outdoor window temperature [K]
B%XABS_WIN = PABS_WIN             ! window absortance
B%XUGG_WIN = PUGG_WIN             ! Window glass-to-glass U-factor [K m W-2]
B%LSHADE = OSHADE                 ! flag to activate shading devices -> LOGICAL in the code
B%XSHADE = PSHADE                 ! flag to activate shading devices -> REAL for i/o 0. or 1.
B%CNATVENT = HNATVENT             ! flag to activate natural ventilation 'NONE', 'MANU', 'AUTO'
B%XNATVENT = PNATVENT             ! flag to describe surventilation system for i/o 
                                  ! 0 for NONE, 1 for MANU and 2 for AUTO
B%LSHAD_DAY = OSHAD_DAY           !Has shading been necessary this day ?
B%LNATVENT_NIGHT = ONATVENT_NIGHT !Has nocturnal surventilation been necessary and possible this night ?
  !
  !indoor relative surfaces and view factors
B%XN_FLOOR        = PN_FLOOR               ! Number of floors     
B%XGLAZ_O_BLD     = PGLAZ_O_BLD         ! Window area [m2_win/m2_bld]
B%XMASS_O_BLD     = PMASS_O_BLD         ! Mass area [m2_mass/m2_bld]
B%XFLOOR_HW_RATIO = PFLOOR_HW_RATIO ! H/W ratio of 1 floor level
B%XF_FLOOR_MASS   = PF_FLOOR_MASS     ! View factor floor-mass
B%XF_FLOOR_WALL   = PF_FLOOR_WALL     ! View factor floor-wall
B%XF_FLOOR_WIN    = PF_FLOOR_WIN       ! View factor floor-window
B%XF_FLOOR_ROOF   = PF_FLOOR_ROOF     ! View factor floor-roof
B%XF_WALL_FLOOR   = PF_WALL_FLOOR     ! View factor wall-floor
B%XF_WALL_MASS    = PF_WALL_MASS      ! View factor wall-mass
B%XF_WALL_WIN     = PF_WALL_WIN         ! View factor wall-win
B%XF_WIN_FLOOR    = PF_WIN_FLOOR       ! View factor win-floor
B%XF_WIN_MASS     = PF_WIN_MASS         ! View factor win-mass
B%XF_WIN_WALL     = PF_WIN_WALL         ! View factor win-wall
B%XF_WIN_WIN      = PF_WIN_WIN           ! indoor View factor win-win
B%XF_MASS_FLOOR   = PF_MASS_FLOOR     ! View factor mass-floor
B%XF_MASS_WALL    = PF_MASS_WALL       ! View factor mass-wall
B%XF_MASS_WIN     = PF_MASS_WIN         ! View factor mass-window
!
! New arguments from 8.2  (Ventilation)
B%XDIAGVENT = PDIAGVENT                     ! Is ventilation active?
B%XFOPEN = PFOPEN                           ! Fraction of windows to be opened
B%XN50 = PN50                               ! Airtightness [AC/H at 50
B%XVENT_BEHAV_ANYWAY = PVENT_BEHAV_ANYWAY   ! Fraction of windows opened in any case
B%XVENT_BEHAV_ADAPTI = PVENT_BEHAV_ADAPTI   ! Fraction of windows that could be opened
B%XTDESV = PTDESV                           ! Desing temperature for ventilation [K]
B%XFVSHO = PFVSHO                           ! Ventilation independent of indoor/outdoor temperature
B%XISMECH = PISMECH                         ! Presence of mechanical ventilation. OFF: 0.
B%XMECHRATE = PMECHRATE                     ! Air exchange rate due to mechanical ventilation [1/h]
!
!-------------------------------------------------------------------------------
!
! TEB OPTIONS structure
!
TOP%LGARDEN   = OGARDEN         ! T: Urban green areas (call ISBA from TEB)
                                ! F: No urban green areas
TOP%CROAD_DIR = HROAD_DIR       ! TEB option for road directions
                                ! 'UNIF' : no specific direction
                                ! 'ORIE' : many road ORIEntations
                                ! ( one per TEB patch)
TOP%CWALL_OPT = HWALL_OPT       ! TEB option for walls
                                ! 'UNIF' : uniform walls
                                ! 'TWO ' : two separated walls
TOP%CZ0H      = HZ0H            ! TEB option for z0h roof & road
                                ! 'MASC95' : Mascart et al 1995
                                ! 'BRUT82' : Brustaert     1982
                                ! 'KAND07' : Kanda         2007
TOP%CCH_BEM   = HCH_BEM         ! BEM option for roof/wall outside convective coefficient
                                ! 'DOE-2' : DOE-2 model from
                                ! EnergyPlus Engineering reference, p65
TOP%CBEM      = HBEM            ! TEB option for the building energy model
                                ! 'DEF':  DEFault version force-restore model from Masson et al. 2002
                                ! 'BEM':  Building Energy Model Bueno et al. 2011

TOP%LGREENROOF   = OGREENROOF   ! T: green roofs (call ISBA from TEB)
TOP%LSOLAR_PANEL = OSOLAR_PANEL ! T: solar panels on roofs
! 
! type of initialization of vegetation: from cover types (ecoclimap) or parameters prescribed
!
!
! Number of layers
!
TOP%NROOF_LAYER = SIZE(PT_ROOF,2)    ! number of layers in roofs
TOP%NROAD_LAYER = SIZE(PT_ROAD,2)    ! number of layers in roads
TOP%NWALL_LAYER = SIZE(PT_WALL_A,2)  ! number of layers in walls
!
! Date:
!
TOP%TTIME = TPTIME        ! current date and time
!
! Time-step:
!
TOP%XTSTEP = PTSTEP       ! time step for TEB
!
!-------------------------------------------------------------------------------
!
! TEB structure
!
! Geometric Parameters:
!
T%XROAD_DIR      = PROAD_DIR      ! Road direction (deg from North, clockwise)
T%XGARDEN        = PGARDEN        ! fraction of veg in the streets   (-)
T%XGREENROOF     = PFRAC_GR       ! fraction of greenroofs on roofs  (-)
T%XBLD           = PBLD           ! fraction of buildings            (-)
T%XROAD          = PROAD          ! fraction of roads                (-)
T%XCAN_HW_RATIO  = PCAN_HW_RATIO  ! canyon    h/W                    (-)
T%XBLD_HEIGHT    = PBLD_HEIGHT    ! buildings height 'h'             (m)
T%XWALL_O_HOR    = PWALL_O_HOR    ! wall surf. / hor. surf.          (-)
T%XROAD_O_GRND   = PROAD_O_GRND   ! road surf. / (road + garden surf.) (-)
T%XGARDEN_O_GRND = PGARDEN_O_GRND ! gard. surf. / (road + garden surf.)(-)
T%XWALL_O_GRND   = PWALL_O_GRND   ! wall surf. / (road + garden surf.) (-)
T%XWALL_O_BLD    = PWALL_O_BLD    ! wall surf. / bld surf. (-)
T%XZ0_TOWN       = PZ0_TOWN       ! roughness length for momentum    (m)
T%XSVF_ROAD      = PSVF_ROAD      ! road sky view factor             (-)
T%XSVF_GARDEN    = PSVF_GARDEN    ! green area sky view factor       (-)
T%XSVF_WALL      = PSVF_WALL      ! wall sky view factor             (-)
!
! Roof parameters
!
T%XALB_ROOF  = PALB_ROOF     ! roof albedo                      (-)
T%XEMIS_ROOF = PEMIS_ROOF    ! roof emissivity                  (-)
T%XHC_ROOF   = PHC_ROOF      ! roof layers heat capacity        (J/K/m3)
T%XTC_ROOF   = PTC_ROOF      ! roof layers thermal conductivity (W/K/m)
T%XD_ROOF    = PD_ROOF       ! depth of roof layers             (m)
T%XROUGH_ROOF= PROUGH_ROOF   ! roof roughness coef
!
!
! Road parameters
!
T%XALB_ROAD  = PALB_ROAD     ! road albedo                      (-)
T%XEMIS_ROAD = PEMIS_ROAD    ! road emissivity                  (-)
T%XHC_ROAD   = PHC_ROAD      ! road layers heat capacity        (J/K/m3)
T%XTC_ROAD   = PTC_ROAD      ! road layers thermal conductivity (W/K/m)
T%XD_ROAD    = PD_ROAD       ! depth of road layers             (m)
!
! Wall parameters
!
T%XALB_WALL  = PALB_WALL     ! wall albedo                      (-)
T%XEMIS_WALL = PEMIS_WALL    ! wall emissivity                  (-)
T%XHC_WALL   = PHC_WALL      ! wall layers heat capacity        (J/K/m3)
T%XTC_WALL   = PTC_WALL      ! wall layers thermal conductivity (W/K/m)
T%XD_WALL    = PD_WALL       ! depth of wall layers             (m)
T%XROUGH_WALL= PROUGH_WALL   ! wall roughness coef
!
! Building's use type
!
T%XRESIDENTIAL = PRESIDENTIAL ! fraction of Residential use      (-)
T%XDT_RES = PDT_RES           ! target temperature change when unoccupied (K) (residential buildings)
T%XDT_OFF = PDT_OFF           ! target temperature change when unoccupied (K) (offices buildings)
  
!
! anthropogenic fluxes
!
T%XH_TRAFFIC   = PH_TRAFFIC   ! anthropogenic sensible
!                             ! heat fluxes due to traffic       (W/m2)
T%XLE_TRAFFIC  = PLE_TRAFFIC  ! anthropogenic latent
!                             ! heat fluxes due to traffic       (W/m2)
T%XH_INDUSTRY  = PH_INDUSTRY  ! anthropogenic sensible                   
!                             ! heat fluxes due to factories     (W/m2)
T%XLE_INDUSTRY = PLE_INDUSTRY ! anthropogenic latent
!                             ! heat fluxes due to factories     (W/m2)
!
! Prognostic variables:
!
T%XWS_ROOF = PWS_ROOF      ! roof water reservoir             (kg/m2)
T%XWS_ROAD = PWS_ROAD      ! road water reservoir             (kg/m2)
T%XT_ROOF  = PT_ROOF       ! roof layer temperatures          (K)
T%XT_ROAD  = PT_ROAD       ! road layer temperatures          (K)
T%XT_WALL_A= PT_WALL_A     ! wall layer temperatures          (K)
T%XT_WALL_B= PT_WALL_B     ! wall layer temperatures          (K)
!
!
!
! Semi-prognostic variables:
!
T%XT_CANYON   = PT_CAN  ! canyon air temperature           (K)
T%XQ_CANYON   = PQ_CAN  ! canyon air specific humidity     (kg/kg)
!
!
! Prognostic snow on roofs:
!
T%TSNOW_ROOF%SCHEME = HSNOW_ROOF    ! snow scheme
T%TSNOW_ROOF%WSNOW  = PWSNOW_ROOF   ! snow (& liq. water) content      (kg/m2) 
T%TSNOW_ROOF%T      = PTSNOW_ROOF   ! temperature                      (K)
T%TSNOW_ROOF%RHO    = PRSNOW_ROOF   ! density                          (kg m-3)
T%TSNOW_ROOF%ALB    = PASNOW_ROOF   ! albedo                           (-)
T%TSNOW_ROOF%EMIS   = PESNOW_ROOF   ! emissivity                       (-)
T%TSNOW_ROOF%TS     = PTSSNOW_ROOF  ! Surface temperature              (K)
!
! Prognostic snow on roads:
!
T%TSNOW_ROAD%SCHEME = HSNOW_ROAD    ! snow scheme
T%TSNOW_ROAD%WSNOW  = PWSNOW_ROAD   ! snow (& liq. water) content      (kg/m2) 
T%TSNOW_ROAD%T      = PTSNOW_ROAD   ! temperature                      (K)
T%TSNOW_ROAD%RHO    = PRSNOW_ROAD   ! density                          (kg m-3)
T%TSNOW_ROAD%ALB    = PASNOW_ROAD   ! albedo                           (-)
T%TSNOW_ROAD%EMIS   = PESNOW_ROAD   ! emissivity                       (-)
T%TSNOW_ROAD%TS     = PTSSNOW_ROAD  ! Surface temperature              (K)

!-------------------------------------------------------------------------------
!
! solar panels
!
TPN%XEMIS_PANEL      = PEMIS_PANEL      ! emissivity of solar panels               (-)
TPN%XALB_PANEL       = PALB_PANEL       ! albedo of soolar panels                  (-)
TPN%XEFF_PANEL       = PEFF_PANEL       ! Efficiency coefficient of solar panels   (-)
TPN%XFRAC_PANEL      = PFRAC_PANEL      ! Fraction of SOLAR panels on roofs        (-)
TPN%XTHER_PRODC_DAY  = PTHER_PRODC_DAY  ! Present day integrated thermal production of energy    (J/m2 panel)


!-------------------------------------------------------------------------------
!
! Irrigation
!
TIR%LPAR_RD_IRRIG   = OPAR_RD_IRRIG   ! flag to use prescribed irrigation for roads
TIR%XRD_START_MONTH = PRD_START_MONTH ! roads : start month for irrigation (included)
TIR%XRD_END_MONTH   = PRD_END_MONTH   ! roads : end   month for irrigation (included)
TIR%XRD_START_HOUR  = PRD_START_HOUR  ! roads : start solar hour for irrigation (included, hour)
TIR%XRD_END_HOUR    = PRD_END_HOUR    ! roads : end   solar hour for irrigation (excluded, hour)
TIR%XRD_24H_IRRIG   = PRD_24H_IRRIG   ! roads : total irrigation over 24 hours (kg/m2)
!
!-------------------------------------------------------------------------------
!
! TEB diagnostics structure
!
DMT%XZ0_TOWN = PZ0_TOWN   ! town roughness length
!-------------------------------------------------------------------------------
!
CALL TEB_GARDEN           (TOP, T, BOP, B, TPN, TIR, DMT,                                                          &
                           HIMPLICIT_WIND, PTSUN, PT_CAN, PQ_CAN, PU_CAN, PT_LOWCAN, PQ_LOWCAN,                    &
                           PU_LOWCAN, PZ_LOWCAN, PPEW_A_COEF, PPEW_B_COEF, PPEW_A_COEF_LOWCAN,                     &
                           PPEW_B_COEF_LOWCAN, PPS, PPA, PEXNS, PEXNA, PTA, PQA, PRHOA, PCO2,                      &
                           PLW_RAD, PDIR_SW, PSCA_SW, PSW_BANDS, KSW, PZENITH, PAZIM, PRR, PSR,                    &
                           PZREF, PUREF, PVMOD, PH_TRAFFIC, PLE_TRAFFIC, PTSTEP, PLEW_ROOF, PLEW_ROAD,             &
                           PLE_WALL_A, PLE_WALL_B, PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,PMELT_ROOF,&
                           PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD, PMELT_ROAD, PRN_GRND, PH_GRND,    &
                           PLE_GRND, PGFLUX_GRND, PRN_TOWN, PH_TOWN, PLE_TOWN, PGFLUX_TOWN, PEVAP_TOWN,            &
                           PSFCO2, PUW_GRND, PUW_ROOF, PDUWDU_GRND, PDUWDU_ROOF,                                   &
                           PUSTAR_TOWN, PCD, PCDN, PCH_TOWN, PRI_TOWN, PTS_TOWN, PEMIS_TOWN,                       &
                           PDIR_ALB_TOWN, PSCA_ALB_TOWN, PRESA_TOWN, PAC_ROAD, PAC_GARDEN, PAC_GREENROOF,          &
                           PAC_ROAD_WAT, PAC_GARDEN_WAT, PAC_GREENROOF_WAT, KDAY, PEMIT_LW_FAC, PEMIT_LW_GRND,     &
                           PT_RAD_IND, PREF_SW_GRND, PREF_SW_FAC, PHU_BLD, PTIME, PPROD_BLD                        )
!
!-------------------------------------------------------------------------------
!
! update of BEM prognostic variables
!
 PTI_BLD = B%XTI_BLD 
 PQI_BLD = B%XQI_BLD 
 PT_FLOOR = B%XT_FLOOR
 PT_MASS = B%XT_MASS
 PT_WIN1 = B%XT_WIN1
 PT_WIN2 = B%XT_WIN2
!
 PAUX_MAX = B%XAUX_MAX
 PCAP_SYS_RAT = B%XCAP_SYS_RAT
 PM_SYS_RAT = B%XM_SYS_RAT
!
 PALB_WIN = B%XALB_WIN             ! window albedo
!
 OSHAD_DAY = B%LSHAD_DAY
 ONATVENT_NIGHT = B%LNATVENT_NIGHT
!
! update of TEB prognostic variables
!
 PT_ROOF = T%XT_ROOF
 PT_ROAD = T%XT_ROAD
 PT_WALL_A = T%XT_WALL_A
 PT_WALL_B = T%XT_WALL_B
 PWS_ROOF  = T%XWS_ROOF
 PWS_ROAD  = T%XWS_ROAD
!
 PWSNOW_ROOF = T%TSNOW_ROOF%WSNOW
 PTSNOW_ROOF = T%TSNOW_ROOF%T
 PRSNOW_ROOF = T%TSNOW_ROOF%RHO
 PASNOW_ROOF = T%TSNOW_ROOF%ALB
 PESNOW_ROOF = T%TSNOW_ROOF%EMIS
 PTSSNOW_ROOF = T%TSNOW_ROOF%TS
!
 PWSNOW_ROAD = T%TSNOW_ROAD%WSNOW
 PTSNOW_ROAD = T%TSNOW_ROAD%T
 PRSNOW_ROAD = T%TSNOW_ROAD%RHO
 PASNOW_ROAD = T%TSNOW_ROAD%ALB
 PESNOW_ROAD = T%TSNOW_ROAD%EMIS
 PTSSNOW_ROAD = T%TSNOW_ROAD%TS
!
! Warning: Can temperature and humidity, XT_CANYON and XQ_CANYON are present in the structure, but 
!          are not updated in teb_garden and routines called by it. This is the PT_CAN and PQ_CAN
!          arguments of teb_garden routines that are physically updated. 
!          So we update the variables in the structure, but this has no consequences, because
!          these variables are not used in teb_garden.
!
 T%XT_CANYON = PT_CAN
 T%XQ_CANYON = PQ_CAN
!
!-------------------------------------------------------------------------------
!
! output diagnostics in TEB structure
!
!-------------------------------------------------------------------------------
!
! solar panels
!
 PTHER_PRODC_DAY = TPN%XTHER_PRODC_DAY   ! Present day integrated thermal production of energy    (J/m2 panel)
!
!-------------------------------------------------------------------------------
!
! TEB diagnostics structure
!
 PQF_BLD  = DMT%XQF_BLD   ! domestic heating
 PFLX_BLD = DMT%XFLX_BLD  ! heat flux from bld
 PQF_TOWN = DMT%XQF_TOWN  ! total anthropogenic heat
 PDQS_TOWN= DMT%XDQS_TOWN ! storage inside building
!
 PH_WALL_A   = DMT%XH_WALL_A ! wall sensible heat flux          (W/m2)
 PH_WALL_B   = DMT%XH_WALL_B ! wall sensible heat flux          (W/m2)
 PH_ROOF     = DMT%XH_ROOF   ! roof sensible heat flux          (W/m2)
 PH_ROAD     = DMT%XH_ROAD   ! road sensible heat flux          (W/m2)
 PH_STRLROOF = DMT%XH_STRLROOF ! structural roof sens. heat flux  (W/m2)
 PH_BLT      = DMT%XH_BLT      ! built surf sensible heat flux    (W/m2)
 PRN_WALL_A  = DMT%XRN_WALL_A  ! net radiation at wall            (W/m2)
 PRN_WALL_B  = DMT%XRN_WALL_B  ! net radiation at wall            (W/m2)
 PRN_ROOF    = DMT%XRN_ROOF    ! net radiation at roof            (W/m2)
 PRN_ROAD    = DMT%XRN_ROAD    ! net radiation at road            (W/m2)
 PRN_STRLROOF= DMT%XRN_STRLROOF!net radiation at structural roofs(W/m2)
 PRN_BLT     = DMT%XRN_BLT     ! net radiation at built surf      (W/m2)
 PGFLUX_WALL_A = DMT%XGFLUX_WALL_A !net wall conduction flux        (W/m2)
 PGFLUX_WALL_B = DMT%XGFLUX_WALL_B !net wall conduction flux        (W/m2)
 PGFLUX_ROOF   = DMT%XGFLUX_ROOF   ! net roof conduction flux         (W/m2)                                         
 PGFLUX_ROAD   = DMT%XGFLUX_ROAD   ! net road conduction flux         (W/m2)
 PGFLUX_STRLROOF = DMT%XGFLUX_STRLROOF !net structural roof cond flux (W/m2)
 PGFLUX_BLT      = DMT%XGFLUX_BLT! net built surf conduction flux   (W/m2)
 PLE_ROOF        = DMT%XLE_ROOF  ! roof latent heat flux            (W/m2)
 PLE_ROAD        = DMT%XLE_ROAD ! road latent heat flux            (W/m2)
 PLE_STRLROOF    = DMT%XLE_STRLROOF !structural roof latent heat flux (W/m2)
 PLE_BLT         = DMT%XLE_BLT      ! built surf latent heat flux      (W/m2)
!
 PRUNOFF_TOWN    = DMT%XRUNOFF_TOWN  ! aggregated water runoff for town      (kg/m2/s)
 PRUNOFF_ROAD    = DMT%XRUNOFF_ROAD  ! water runoff for roads                (kg/m2/s)
 PIRRIG_ROAD     = DMT%XIRRIG_ROAD  ! road man-made watering rate           (kg/m2/s)
 PRUNOFF_ROOF    = DMT%XRUNOFF_ROOF  ! aggregated water runoff for roofs     (kg/m2/s)
 PRUNOFF_STRLROOF= DMT%XRUNOFF_STRLROOF  ! water runoff for structural roofs     (kg/m2/s)
!
 PABS_SW_ROOF     = DMT%XABS_SW_ROOF ! absorbed shortwave radiation over roofs
 PABS_SW_SNOW_ROOF= DMT%XABS_SW_SNOW_ROOF ! absorbed longwave radiation over roofs
 PABS_LW_ROOF     = DMT%XABS_LW_ROOF ! absorbed shortwave radiation over roofs
 PABS_LW_SNOW_ROOF= DMT%XABS_LW_SNOW_ROOF! absorbed longwave radiation over roofs
 PABS_SW_ROAD     = DMT%XABS_SW_ROAD ! absorbed shortwave radiation over roads
 PABS_SW_SNOW_ROAD= DMT%XABS_SW_SNOW_ROAD ! absorbed longwave radiation over roads
 PABS_LW_ROAD     = DMT%XABS_LW_ROAD ! absorbed shortwave radiation over roads
 PABS_LW_SNOW_ROAD= DMT%XABS_LW_SNOW_ROAD ! absorbed longwave radiation over roads
 PABS_SW_WALL_A   = DMT%XABS_SW_WALL_A ! absorbed shortwave radiation over walls
 PABS_SW_WALL_B   = DMT%XABS_SW_WALL_B ! absorbed shortwave radiation over walls
 PABS_LW_WALL_A   = DMT%XABS_LW_WALL_A ! absorbed shortwave radiation over walls
 PABS_LW_WALL_B   = DMT%XABS_LW_WALL_B ! absorbed shortwave radiation over walls
 PABS_SW_GARDEN   = DMT%XABS_SW_GARDEN ! absorbed shortwave radiation over green areas
 PABS_LW_GARDEN   = DMT%XABS_LW_GARDEN ! absorbed shortwave radiation over green areas
 PABS_SW_GREENROOF= DMT%XABS_SW_GREENROOF ! absorbed shortwave radiation over green roofs
 PABS_LW_GREENROOF= DMT%XABS_LW_GREENROOF ! absorbed shortwave radiation over green roofs
 PG_GREENROOF_ROOF= DMT%XG_GREENROOF_ROOF ! Heat flux between green roof and structural roof
 PABS_SW_PANEL    = DMT%XABS_SW_PANEL ! absorbed shortwave radiation over solar panels
 PABS_LW_PANEL    = DMT%XABS_LW_PANEL ! absorbed longwave  radiation over solar panels
!
 PRN_PANEL        = DMT%XRN_PANEL ! net radiation           over solar panels (W/m2)
 PH_PANEL         = DMT%XH_PANEL ! sensible heat flux      over solar panels (W/m2)
 PTHER_PROD_PANEL = DMT%XTHER_PROD_PANEL ! thermal      production of   solar panels (W/m2 thermal panel)
 PPHOT_PROD_PANEL = DMT%XPHOT_PROD_PANEL ! photovoltaic production of   solar panels (W/m2 photovoltaic panel)
 PPROD_PANEL      = DMT%XPROD_PANEL !              production of   solar panels (W/m2 panel)
 PTHER_PROD_BLD   = DMT%XTHER_PROD_BLD ! thermal      production of   solar panels (W/m2 bld)
 PPHOT_PROD_BLD   = DMT%XPHOT_PROD_BLD ! photovoltaic production of   solar panels (W/m2 bld)

 PH_BLD_COOL      = DMT%XH_BLD_COOL ! Sensible cooling energy demand  
                                                     ! of the building [W m-2(bld)]
 PT_BLD_COOL      = DMT%XT_BLD_COOL ! Total cooling energy demand  
                                                     ! of the building [W m-2(bld)]
 PH_BLD_HEAT      = DMT%XH_BLD_HEAT ! Heating energy demand       
                                                     ! of the building [W m-2(bld)]
 PLE_BLD_COOL     = DMT%XLE_BLD_COOL ! Latent cooling energy demand 
                                                     ! of the building [W m-2(bld)]
 PLE_BLD_HEAT     = DMT%XLE_BLD_HEAT ! Latent heating energy demand 
                                                     ! of the building [W m-2(bld)]
 PH_WASTE         = DMT%XH_WASTE ! Sensible waste heat from HVAC system
                                                     ! [W m-2(tot)]
 PLE_WASTE        = DMT%XLE_WASTE ! Latent waste heat from HVAC system
                                                     ! [W m-2(tot)]
 PHVAC_COOL       = DMT%XHVAC_COOL ! Energy consumption of the cooling system
                                                     ! [W m-2(bld)]
 PHVAC_HEAT       = DMT%XHVAC_HEAT ! Energy consumption of the heating system
                                                     ! [W m-2(bld)]
 PCAP_SYS         = DMT%XCAP_SYS ! Actual capacity of the cooling system
                                                     ! [W m-2(bld)] 
 PM_SYS           = DMT%XM_SYS ! Actual HVAC mass flow rate 
                                                     ! [kg s-1 m-2(bld)]
 PCOP             = DMT%XCOP ! COP of the cooling system
 PQ_SYS           = DMT%XQ_SYS ! Supply air specific humidity [kg kg-1]
 PT_SYS           = DMT%XT_SYS ! Supply air temperature [K]
 PTR_SW_WIN       = DMT%XTR_SW_WIN ! Solar radiation transmitted throught
                                                     ! windows [W m-2(bld)]
 PFAN_POWER       = DMT%XFAN_POWER ! HVAC fan power
 PABS_SW_WIN      = DMT%XABS_SW_WIN ! window absorbed shortwave radiation [W m-2] 
 PABS_LW_WIN      = DMT%XABS_LW_WIN ! absorbed infrared rad. [W m-2]

 PCUR_TCOOL_TARGET = DMT%XTCOOL_TARGET ! Cooling target temperature at current time (K)
 PCUR_THEAT_TARGET = DMT%XTHEAT_TARGET ! Heating target temperature at current time (K)
 PCUR_QIN          = DMT%XQIN          ! Internal heat gains        at current time (W/m2 floor)
!
!-------------------------------------------------------------------------------
!
CALL DEALLOC_TEB_STRUCT(TOP,T,BOP,B,TPN,TIR,DMT) 
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE TEB_GARDEN_STRUCT
