!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!   ##########################################################################
    SUBROUTINE TEB  (HZ0H, HIMPLICIT_WIND, HWALL_OPT, HBEM, TPTIME, PTSUN,     &
                     PT_CANYON, PQ_CANYON, PU_CANYON,                          &
                     PT_LOWCAN, PQ_LOWCAN, PU_LOWCAN, PZ_LOWCAN, PTI_BLD,      &
                     PT_ROOF, PT_ROAD, PT_WALL_A, PT_WALL_B,PWS_ROOF, PWS_ROAD,&
                     HSNOW_ROOF, PWSNOW_ROOF, PTSNOW_ROOF, PRSNOW_ROOF,        &
                     PASNOW_ROOF, PTSSNOW_ROOF, PESNOW_ROOF,                   &
                     HSNOW_ROAD, PWSNOW_ROAD, PTSNOW_ROAD, PRSNOW_ROAD,        &
                     PASNOW_ROAD, PTSSNOW_ROAD, PESNOW_ROAD,                   &
                     PPEW_A_COEF, PPEW_B_COEF,                                 &
                     PPEW_A_COEF_LOWCAN, PPEW_B_COEF_LOWCAN,                   &
                     PPS, PPA, PEXNS, PEXNA, PTA, PQA, PRHOA, PLW_RAD,         &
                     PRR, PSR, PZREF, PUREF,                                   &
                     PVMOD, PH_TRAFFIC, PLE_TRAFFIC, PH_INDUSTRY, PLE_INDUSTRY,&
                     PTSTEP, PZ0_TOWN, PBLD, PGARDEN, PROAD, PFRAC_GR,         &
                     PBLD_HEIGHT, PWALL_O_HOR, PCAN_HW_RATIO, PWALL_O_GRND,    &
                     PDF_ROOF, PDN_ROOF, PDF_ROAD,                             &
                     PDN_ROAD, PQSAT_ROOF, PQSAT_ROAD, PDELT_ROOF, PDELT_ROAD, &
                     PEMIS_ROOF, PHC_ROOF, PTC_ROOF, PD_ROOF,                  &
                     PEMIS_ROAD, PHC_ROAD, PTC_ROAD,                           &
                     PD_ROAD, PEMIS_WALL, PTS_GARDEN,                          &
                     PHC_WALL, PTC_WALL, PD_WALL, PRN_ROOF, PH_ROOF, PLE_ROOF, &
                     PLEW_ROOF, PGFLUX_ROOF, PRUNOFF_ROOF,                     &
                     PRN_GREENROOF, PH_GREENROOF, PLE_GREENROOF,               &
                     PGFLUX_GREENROOF, PUW_GREENROOF,                          &
                     PRUNOFF_GREENROOF, PDRAIN_GREENROOF,                      &
                     PRN_STRLROOF, PH_STRLROOF, PLE_STRLROOF, PGFLUX_STRLROOF, &
                     PRUNOFF_STRLROOF,                                         &
                     PRN_ROAD, PH_ROAD,                                        &
                     PLE_ROAD, PLEW_ROAD, PGFLUX_ROAD, PRUNOFF_ROAD,           &
                     PRN_WALL_A, PH_WALL_A, PLE_WALL_A, PGFLUX_WALL_A,         &
                     PRN_WALL_B, PH_WALL_B, PLE_WALL_B, PGFLUX_WALL_B,         &
                     PRN_BLT, PH_BLT, PLE_BLT, PGFLUX_BLT,                     &
                     PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,     &
                     PMELT_ROOF,                                               &
                     PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,     &
                     PMELT_ROAD,                                               &
                     PG_GREENROOF_ROOF,                                        &
                     PUW_ROAD, PUW_ROOF, PDUWDU_ROAD, PDUWDU_ROOF,             &
                     PUSTAR_TOWN, PCD, PCDN, PCH_TOWN, PRI_TOWN,               &
                     PRESA_TOWN, PDQS_TOWN, PQF_TOWN, PQF_BLD, PFLX_BLD,       &
                     PAC_ROOF, PAC_ROAD, PAC_WALL, PAC_TOP, PAC_GARDEN,        &
                     PAC_ROOF_WAT, PAC_ROAD_WAT, PABS_SW_ROOF, PABS_LW_ROOF,   &
                     PABS_SW_SNOW_ROOF, PABS_LW_SNOW_ROOF, PABS_SW_ROAD,       &
                     PABS_LW_ROAD, PABS_SW_SNOW_ROAD, PABS_LW_SNOW_ROAD,       &
                     PABS_SW_WALL_A,PABS_LW_WALL_A,                            &
                     PABS_SW_WALL_B,PABS_LW_WALL_B,                            &
                     PLW_WA_TO_WB,                                             &
                     PLW_WA_TO_R, PLW_WB_TO_R,                                 &
                     PLW_WA_TO_NR, PLW_WB_TO_NR,                               &
                     PLW_R_TO_WA, PLW_R_TO_WB,                                 &
                     PLW_G_TO_WA, PLW_G_TO_WB,                                 &
                     PLW_S_TO_WA, PLW_S_TO_WB, PLW_S_TO_R,                     &
                     PLW_S_TO_NR, PLW_NR_TO_WA, PLW_NR_TO_WB,                  &
                     PLW_NR_TO_WIN, PLW_WA_TO_WIN, PLW_WB_TO_WIN, PLW_G_TO_WIN,&
                     PLW_R_TO_WIN, PLW_S_TO_WIN, PLW_WIN_TO_WA, PLW_WIN_TO_WB, &
                     PLW_WIN_TO_R, PLW_WIN_TO_NR,                              &
                     HNATVENT,                                                 &
                     HCOOL_COIL, PF_WATER_COND, HHEAT_COIL, OAUTOSIZE,         &
                     KDAY, PAUX_MAX, PT_FLOOR, PT_MASS, PH_BLD_COOL,           &
                     PT_BLD_COOL, PH_BLD_HEAT, PLE_BLD_COOL, PLE_BLD_HEAT,     &
                     PH_WASTE, PLE_WASTE, PF_WASTE_CAN, PHVAC_COOL, PHVAC_HEAT,&
                     PQIN, PQIN_FRAD, PQIN_FLAT, PGR, PEFF_HEAT,               &
                     PINF, PTCOOL_TARGET, PTHEAT_TARGET, PHR_TARGET, PT_WIN2,  &
                     PQI_BLD, PV_VENT, PCAP_SYS_HEAT, PCAP_SYS_RAT, PT_ADP,    &
                     PM_SYS_RAT, PCOP_RAT, PCAP_SYS, PM_SYS, PCOP,  PQ_SYS,    &
                     PT_SYS, PTR_SW_WIN, PFAN_POWER, PHC_FLOOR, PTC_FLOOR,     &
                     PD_FLOOR,                                                 &
                     PT_WIN1, PABS_SW_WIN, PABS_LW_WIN, PUGG_WIN, PEMIT_LW_FAC,&
                     PEMIT_LW_ROAD, PT_RAD_IND, PHU_BLD, PTIME, PE_SHADING,    &
                     ONATVENT_NIGHT, PN_FLOOR, PWALL_O_BLD, PGLAZ_O_BLD, PMASS_O_BLD, &
                     PFLOOR_HW_RATIO, PF_FLOOR_MASS, PF_FLOOR_WALL,            &
                     PF_FLOOR_WIN, PF_FLOOR_ROOF, PF_WALL_FLOOR, PF_WALL_MASS, &
                     PF_WALL_WIN, PF_WIN_FLOOR, PF_WIN_MASS, PF_WIN_WALL,      &
                     PF_MASS_FLOOR, PF_MASS_WALL, PF_MASS_WIN, OCANOPY,        &
                     HCH_BEM, PROUGH_ROOF, PROUGH_WALL, PF_WIN_WIN,            &
                     OPAR_RD_IRRIG, PRD_START_MONTH, PRD_END_MONTH,            &
                     PRD_START_HOUR, PRD_END_HOUR, PRD_24H_IRRIG, PIRRIG_ROAD  )
!   ##########################################################################
!
!!****  *TEB*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the evoultion of prognostic variables and the fluxes
!     over artificial surfaces as towns, taking into account the canyon like
!     geometry of urbanized areas.
!         
!     
!!**  METHOD
!     ------
!
!     The prognostic variables are:
!       - the surface temperature for roofs, roads, and walls
!       - the water reservoir, whose maximum value is 10mm
!
!
!    1 : Warning about snow
!        ******************
!
!     Except for snow mantel evolution, all other computation with snow
!   variables must be performed with these variables at previous time-step,
!   and NOT new time-step. This insure coherence between snow fractions
!   (computed at the begining) and other snow characteristics (albedo, Ts).
!
!
!    2 : computation of input solar radiation on each surface
!        ****************************************************
!
!      Those are now done in subroutine urban_solar_abs.F90
!
!    3 : drag coefficient for momentum 
!        *****************************
!
!
!    4 : aerodynamical resistance for heat transfers
!        *******************************************
!
!
!    5 : equation for evolution of Ts_roof
!        *********************************
!
!
!       Rn = (dir_Rg + sca_Rg) (1-a) + emis * ( Rat - sigma Ts**4 (t+dt) )
!
!       H  = rho Cp CH V ( Ts (t+dt) - Tas )
!
!       LE = rho Lv CH V ( qs (t+dt) - qas )
!
!      where the as subscript denotes atmospheric values at ground level
!      (and not at first half level)
!
!
!    6 : equations for evolution of Ts_road and Ts_wall simultaneously
!        *************************************************************
!
!
!
!   Rn_w = abs_Rg_w 
!  - sigma * emis_w                                                   * Ts_w**4 (t+dt)
!  +         emis_w                       *      SVF_w                * Rat
!  + sigma * emis_w * emis_r              *      SVF_w                * Ts_r**4 (t+dt)
!  + sigma * emis_w * emis_w              * (1-2*SVF_w)               * Ts_w**4 (t+dt)
!  + sigma * emis_w * emis_w * (1-emis_r) *      SVF_w  * (1-  SVF_r) * Ts_w**4 (t+dt)
!  + sigma * emis_w * emis_w * (1-emis_w) * (1-2*SVF_w) * (1-2*SVF_w) * Ts_w**4 (t+dt)
!  + sigma * emis_w * emis_r * (1-emis_w) *      SVF_w  * (1-2*SVF_w) * Ts_r**4 (t+dt)
!
!   Rn_r = abs_Rg_r
!  - sigma * emis_r                                                   * Ts_r**4 (t+dt)
!  +         emis_r                       *    SVF_r                  * Rat
!  + sigma * emis_r * emis_w              * (1-SVF_r)                 * Ts_w**4 (t+dt)
!  + sigma * emis_r * emis_w * (1-emis_w) * (1-SVF_r)   * (1-2*SVF_w) * Ts_w**4 (t+dt)
!  + sigma * emis_r * emis_r * (1-emis_w) * (1-SVF_r)   *      SVF_w  * Ts_r**4 (t+dt)
!
!  H_w  = rho Cp CH V ( Ts_w (t+dt) - Ta_canyon )
!
!  LE_w = rho Lv CH V ( qs_w (t+dt) - qa_canyon )
!
!  H_r  = rho Cp CH V ( Ts_r (t+dt) - Ta_canyon )
!
!  LE_r = rho Lv CH V ( qs_r (t+dt) - qa_canyon )
!
! with again
!                AC_can * Swall/Sroad * Twall + AC_can * Troad + AC_top * Ta + H_traffic/Cp/rho/Sroad
!   Ta_canyon = -------------------------------------------------------------------------------------
!                AC_can * Swall/Sroad         + AC_can         + AC_top
!
!
!                 AC_can * delt_road * Hu_road * qsat(Troad) + AC_top * qa + LE_traffic/Lv/rho/Sroad
!   qa_canyon = ------------------------------------------------------------------------------------
!                 AC_can * delt_road                        + AC_top
!
!
!
!
!    7 : computation of fluxes for each surface type
!        *******************************************
!
!
!    8 : averaging of the fluxes
!        ***********************
!
!   This is done on the total exchange surface (roof + wall + road),
!  which is bigger than the horizontal surface (roof+road), leading
!  to bigger fluxes.
!
!   The fluxes due to industrial activity are directly added into the 
!  atmosphere
!
!
!    9 : road reservoir evolution
!        ************************
!
!   The roof reservoir runoff goes directly into the road reservoir.
!
!   Runoff occurs for road reservoir (too much water), as well as drainage
!   (evacuation system, typical time scale: 1 day)
!
!
!------------------------
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    MODD_CST
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
!!      Original    23/01/98 
!!     21 / 10 / 2003   P. Tulet    output aerodynamical resistance
!!     01 / 07 / 2005   P.Le Moigne Exner functions as arguments to urban_fluxes
!!     17 / 10 / 2005   (G. Pigeon) computation of anthropogenic heat from domestic heating
!!          01 / 2012   V. Masson   Separates the 2 walls 
!!     25 / 09 / 2012   B. Decharme new wind implicitation
!!          07 / 2013   V. Masson   Adds road watering
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TYPE_DATE_SURF,ONLY: DATE_TIME
USE MODD_CSTS,         ONLY : XTT, XSTEFAN, XCPD, XLVTT
USE MODD_SURF_PAR,     ONLY : XUNDEF
USE MODD_SNOW_PAR,     ONLY : XEMISSN, XANSMAX_ROOF, &
                          XANSMAX_ROAD,XWCRN_ROOF,XWCRN_ROAD
!
USE MODE_THERMOS
USE MODE_SURF_SNOW_FRAC
!
USE MODI_SNOW_COVER_1LAYER
USE MODI_URBAN_DRAG
USE MODI_URBAN_SNOW_EVOL
USE MODI_ROOF_LAYER_E_BUDGET
USE MODI_ROAD_LAYER_E_BUDGET
USE MODI_FACADE_E_BUDGET
USE MODI_URBAN_FLUXES
USE MODI_URBAN_HYDRO
USE MODI_BLD_E_BUDGET
USE MODI_WIND_THRESHOLD
USE MODI_BEM
USE MODI_TEB_IRRIG
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    Declarations of arguments
!
!
 CHARACTER(LEN=6), INTENT(IN)      :: HZ0H          ! TEB option for z0h roof & road
!                                                  ! 'MASC95' : Mascart et al 1995
!                                                  ! 'BRUT82' : Brustaert     1982
!                                                  ! 'KAND07' : Kanda         2007
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
!                                                     ! 'OLD' = direct
!                                                     ! 'NEW' = Taylor serie, order 1
 CHARACTER(LEN=4), INTENT(IN)      :: HWALL_OPT     ! TEB option for walls
                                                   ! 'UNIF' : uniform walls
                                                   ! 'TWO ' : two different opposite walls
 CHARACTER(LEN=3), INTENT(IN)      :: HBEM          ! Building Energy model 'DEF' or 'BEM'
TYPE(DATE_TIME)     , INTENT(IN)    :: TPTIME             ! current date and time from teb
REAL, DIMENSION(:),   INTENT(IN)    :: PTSUN              ! solar time   (s from midnight)
REAL, DIMENSION(:), INTENT(INOUT) :: PT_CANYON     ! canyon air temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PQ_CANYON     ! canyon air specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PU_CANYON     ! canyon hor. wind
REAL, DIMENSION(:), INTENT(IN)    :: PU_LOWCAN     ! wind near the road
REAL, DIMENSION(:), INTENT(IN)    :: PT_LOWCAN     ! temp. near the road
REAL, DIMENSION(:), INTENT(IN)    :: PQ_LOWCAN     ! hum. near the road
REAL, DIMENSION(:), INTENT(IN)    :: PZ_LOWCAN     ! height of atm. var. near the road
REAL, DIMENSION(:), INTENT(INOUT) :: PTI_BLD       ! inside building temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_ROOF     ! roof layers temperatures
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_ROAD     ! road layers temperatures
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_WALL_A   ! wall A layers temperatures
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_WALL_B   ! wall B layers temperatures
REAL, DIMENSION(:), INTENT(INOUT) :: PWS_ROOF      ! roof water reservoir
REAL, DIMENSION(:), INTENT(INOUT) :: PWS_ROAD      ! road water reservoir
 CHARACTER(LEN=*),   INTENT(IN)    :: HSNOW_ROOF    ! snow roof scheme
!                                                  ! 'NONE'
!                                                  ! 'D95 '
!                                                  ! '1-L '
 CHARACTER(LEN=*),   INTENT(IN)    :: HSNOW_ROAD    ! snow road scheme
!                                                  ! 'NONE'
!                                                  ! 'D95 '
!                                                  ! '1-L '
REAL, DIMENSION(:,:), INTENT(INOUT) :: PWSNOW_ROOF ! snow layers reservoir
REAL, DIMENSION(:,:), INTENT(INOUT) :: PTSNOW_ROOF ! snow layers temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRSNOW_ROOF ! snow layers density
REAL, DIMENSION(:),   INTENT(INOUT) :: PASNOW_ROOF ! snow albedo
REAL, DIMENSION(:),   INTENT(INOUT) :: PESNOW_ROOF ! snow emissivity
REAL, DIMENSION(:),   INTENT(INOUT) :: PTSSNOW_ROOF! snow surface temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PWSNOW_ROAD ! snow layers reservoir
REAL, DIMENSION(:,:), INTENT(INOUT) :: PTSNOW_ROAD ! snow layers temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRSNOW_ROAD ! snow layers density
REAL, DIMENSION(:),   INTENT(INOUT) :: PASNOW_ROAD ! snow albedo
REAL, DIMENSION(:),   INTENT(INOUT) :: PESNOW_ROAD ! snow emissivity
REAL, DIMENSION(:),   INTENT(INOUT) :: PTSSNOW_ROAD! snow surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_A_COEF   ! implicit coefficients
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_B_COEF   ! for wind coupling
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_A_COEF_LOWCAN ! implicit coefficients for wind coupling
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_B_COEF_LOWCAN ! between low canyon wind and road
REAL, DIMENSION(:), INTENT(IN)    :: PPS           ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PPA           ! pressure at the first atmospheric level
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS         ! surface exner function
REAL, DIMENSION(:), INTENT(IN)    :: PTA           ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA           ! specific humidity
                                                   ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD         ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PEXNA         ! exner function
                                                   ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA         ! air density
                                                   ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD       ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PRR           ! rain rate
REAL, DIMENSION(:), INTENT(IN)    :: PSR           ! snow rate
REAL, DIMENSION(:), INTENT(IN)    :: PH_TRAFFIC    ! anthropogenic sensible
!                                                  ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PLE_TRAFFIC   ! anthropogenic latent
!                                                  ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PH_INDUSTRY   ! anthropogenic sensible
!                                                  ! heat fluxes due to factories
REAL, DIMENSION(:), INTENT(IN)    :: PLE_INDUSTRY  ! anthropogenic latent
!                                                  ! heat fluxes due to factories
REAL, DIMENSION(:), INTENT(IN)    :: PZREF         ! reference height of the first
                                                   ! atmospheric level (temperature)
REAL, DIMENSION(:), INTENT(IN)    :: PUREF         ! reference height of the first
                                                   ! atmospheric level (wind)
REAL,               INTENT(IN)    :: PTSTEP        ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PZ0_TOWN      ! town roughness length
                                                   ! for momentum
REAL, DIMENSION(:), INTENT(IN)    :: PBLD          ! fraction of buildings
REAL, DIMENSION(:), INTENT(IN)    :: PGARDEN       ! fraction of GARDEN areas
REAL, DIMENSION(:), INTENT(IN)    :: PROAD         ! fraction of roads
REAL, DIMENSION(:), INTENT(IN)    :: PFRAC_GR      ! fraction of green roofs
REAL, DIMENSION(:), INTENT(IN)    :: PBLD_HEIGHT   ! buildings h
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_O_HOR   ! wall surf. / hor. surf.
REAL, DIMENSION(:), INTENT(IN)    :: PCAN_HW_RATIO ! canyon    h/W
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_O_GRND  ! wall surf. / (road + gard.  surf.)
!
REAL, DIMENSION(:), INTENT(INOUT) :: PDF_ROOF      ! snow-free    fraction on roofs
REAL, DIMENSION(:), INTENT(INOUT) :: PDN_ROOF      ! snow-covered fraction on roofs
REAL, DIMENSION(:), INTENT(INOUT) :: PDF_ROAD      ! snow-free    fraction on roads
REAL, DIMENSION(:), INTENT(INOUT) :: PDN_ROAD      ! snow-covered fraction on roads
REAL, DIMENSION(:), INTENT(OUT)   :: PQSAT_ROOF    ! hum at saturation over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PQSAT_ROAD    ! hum at saturation over road
REAL, DIMENSION(:), INTENT(OUT)   :: PDELT_ROOF    ! water fraction on roof
REAL, DIMENSION(:), INTENT(OUT)   :: PDELT_ROAD    ! water fraction on road
!
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_ROOF    ! roof emissivity
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_ROOF      ! heat capacity for roof layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_ROOF      ! thermal conductivity for roof layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_ROOF       ! depth of roof layers
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_ROAD    ! road emissivity
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_ROAD      ! heat capacity for road layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_ROAD      ! thermal conductivity for road layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_ROAD       ! depth of road layers
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_WALL    ! wall emissivity
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_WALL      ! heat capacity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_WALL      ! thermal conductivity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_WALL       ! depth of wall layers
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GARDEN    ! GARDEN area surf temp.
! greenroof
REAL, DIMENSION(:), INTENT(IN)    :: PG_GREENROOF_ROOF ! snow melt
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_ROOF     ! net radiation over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PH_ROOF      ! sensible heat flux over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_ROOF     ! latent heat flux over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PLEW_ROOF    ! latent heat flux over roof (snow)
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_ROOF  ! flux through the roof
REAL, DIMENSION(:), INTENT(OUT)   :: PRUNOFF_STRLROOF ! runoff over the ground
REAL, DIMENSION(:), INTENT(OUT)   :: PRUNOFF_ROOF ! runoff over roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_ROAD     ! net radiation over road
REAL, DIMENSION(:), INTENT(OUT)   :: PH_ROAD      ! sensible heat flux over road
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_ROAD     ! latent heat flux over road
REAL, DIMENSION(:), INTENT(OUT)   :: PLEW_ROAD    ! latent heat flux over road (snow)
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_ROAD  ! flux through the road
REAL, DIMENSION(:), INTENT(OUT)   :: PRUNOFF_ROAD ! runoff over the ground
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_WALL_A   ! net radiation over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PH_WALL_A    ! sensible heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WALL_A   ! latent heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_WALL_A! flux through the wall
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_WALL_B   ! net radiation over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PH_WALL_B    ! sensible heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WALL_B   ! latent heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_WALL_B! flux through the wall
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_BLT       ! net radiation over built surf 
REAL, DIMENSION(:), INTENT(OUT)   :: PH_BLT        ! sensible heat flux over built surf 
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_BLT       ! latent heat flux over built surf 
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_BLT    ! flux through the built surf 
!
REAL, DIMENSION(:), INTENT(IN)    :: PRN_GREENROOF     ! net radiation over greenroof
REAL, DIMENSION(:), INTENT(IN)    :: PH_GREENROOF      ! sensible heat flux over greenroof
REAL, DIMENSION(:), INTENT(IN)    :: PLE_GREENROOF     ! latent heat flux over greenroof
REAL, DIMENSION(:), INTENT(IN)    :: PGFLUX_GREENROOF  ! flux through the greenroof
REAL, DIMENSION(:), INTENT(IN)    :: PUW_GREENROOF     ! Momentum flux for greenroofs
REAL, DIMENSION(:), INTENT(IN)    :: PRUNOFF_GREENROOF ! runoff over green roofs
REAL, DIMENSION(:), INTENT(IN)    :: PDRAIN_GREENROOF  ! outlet drainage at base of green roofs
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_STRLROOF      ! net radiation over structural roof
REAL, DIMENSION(:), INTENT(OUT)   :: PH_STRLROOF       ! sensible heat flux over structural roof
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_STRLROOF      ! latent heat flux over structural roof
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_STRLROOF   ! flux through the structural roof
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSNOW_ROOF ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSNOW_ROOF  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESNOW_ROOF ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSNOW_ROOF  ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_ROOF   ! snow melt
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSNOW_ROAD ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSNOW_ROAD  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESNOW_ROAD ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSNOW_ROAD  ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_ROAD   ! snow melt
!
REAL, DIMENSION(:), INTENT(OUT)   :: PUW_ROAD     ! Momentum flux for roads
REAL, DIMENSION(:), INTENT(OUT)   :: PUW_ROOF     ! Momentum flux for roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PDUWDU_ROAD  !
REAL, DIMENSION(:), INTENT(OUT)   :: PDUWDU_ROOF  !
REAL, DIMENSION(:), INTENT(OUT)   :: PUSTAR_TOWN  ! friciton velocity over town
REAL, DIMENSION(:), INTENT(OUT)   :: PCD          ! town averaged drag coefficient
REAL, DIMENSION(:), INTENT(OUT)   :: PCDN         ! town averaged neutral drag coefficient
REAL, DIMENSION(:), INTENT(OUT)   :: PCH_TOWN     ! town averaged heat transfer
!                                                 ! coefficient
REAL, DIMENSION(:), INTENT(OUT)   :: PRI_TOWN      ! town averaged Richardson number
REAL, DIMENSION(:), INTENT(OUT)   :: PRESA_TOWN    ! town aerodynamical resistance
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_TOWN     ! heat storage inside town
REAL, DIMENSION(:), INTENT(OUT)   :: PQF_TOWN      ! total anthropogenic heat
REAL, DIMENSION(:), INTENT(OUT)   :: PQF_BLD       ! anthropogenic heat flux of
                                                   ! domestic heating  
REAL, DIMENSION(:), INTENT(OUT)   :: PFLX_BLD      ! heat flx from inside bld through its structure
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROOF      ! roof conductance
REAL, DIMENSION(:), INTENT(INOUT) :: PAC_ROAD      ! road conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_WALL      ! wall conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_TOP       ! top conductance
REAL, DIMENSION(:), INTENT(IN)    :: PAC_GARDEN    ! garden conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROOF_WAT  ! roof water conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROAD_WAT  ! roof water conductance
!
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROOF      ! absorbed solar rad by roof
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_SNOW_ROOF ! absorbed solar rad by snow on roof
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROAD      ! absorbed solar rad by road
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_SNOW_ROAD ! absorbed solar rad by snow on road
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WALL_A    ! absorbed solar rad by wall
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WALL_B    ! absorbed solar rad by wall
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_ROOF      ! absorbed IR rad by roof
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_SNOW_ROOF ! absorbed IR rad by snow on roof
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_ROAD      ! absorbed IR rad by road
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_SNOW_ROAD ! absorbed IR rad by snow on road
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_WALL_A    ! absorbed IR rad by wall
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_WALL_B    ! absorbed IR rad by wall
!
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_WB      ! LW contrib. wall A (orB) -> wall B (or A)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_R         ! LW contrib. wall       -> road
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WB_TO_R         ! LW contrib. wall       -> road
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_NR        ! LW contrib. wall       -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WB_TO_NR        ! LW contrib. wall       -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_WA         ! LW contrib. road       -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_WB         ! LW contrib. road       -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_WA         ! LW contrib. GARDEN     -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_WB         ! LW contrib. GARDEN     -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_WA        ! LW contrib. road(snow) -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_WB        ! LW contrib. road(snow) -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_WA         ! LW contrib. sky        -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_WB         ! LW contrib. sky        -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_R          ! LW contrib. sky        -> road
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_NR         ! LW contrib. sky        -> road(snow)
!
! new arguments after BEM
!
 CHARACTER(LEN=4), DIMENSION(:), INTENT(IN) :: HNATVENT
 CHARACTER(LEN=6), INTENT(IN)       :: HCOOL_COIL
REAL, DIMENSION(:), INTENT(IN)     :: PF_WATER_COND!fraction of evaporation for the condensers
 CHARACTER(LEN=6), INTENT(IN)       :: HHEAT_COIL
LOGICAL,            INTENT(IN)     :: OAUTOSIZE    ! Flag to activate autosize calculations
INTEGER,            INTENT(IN)     :: KDAY         ! Simulation day
REAL, DIMENSION(:), INTENT(INOUT)  :: PAUX_MAX     ! Auxiliar variable for autosize calcs
REAL, DIMENSION(:,:), INTENT(INOUT):: PT_FLOOR     ! Floor layers temperatures [K]
REAL, DIMENSION(:,:), INTENT(INOUT):: PT_MASS      ! Internal mass layers temperatures [K]
REAL, DIMENSION(:), INTENT(OUT)    :: PH_BLD_COOL  ! Sensible cooling energy demand  
                                                   ! of the building [W m-2(bld)]
REAL, DIMENSION(:), INTENT(OUT)   :: PT_BLD_COOL   ! Total cooling energy demand  
                                                   ! of the building [W m-2(bld)]
REAL, DIMENSION(:), INTENT(OUT)   :: PH_BLD_HEAT   ! Heating energy demand       
                                                   ! of the building [W m-2(bld)]
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_BLD_COOL  ! Latent cooling energy demand 
                                                   ! of the building [W m-2(bld)]
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_BLD_HEAT  ! Latent heating energy demand 
                                                   ! of the building [W m-2(bld)]
REAL, DIMENSION(:), INTENT(OUT)   :: PH_WASTE      ! Sensible waste heat from HVAC system
                                                   ! [W m-2(tot)]
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WASTE     ! Latent waste heat from HVAC system
                                                   ! [W m-2(tot)]
REAL, DIMENSION(:), INTENT(IN)    :: PF_WASTE_CAN  ! fraction of waste heat released into the canyon
REAL, DIMENSION(:), INTENT(OUT)   :: PHVAC_COOL    ! Energy consumption of the cooling system
                                                   ! [W m-2(bld)]
REAL, DIMENSION(:), INTENT(OUT)  :: PHVAC_HEAT     ! Energy consumption of the heating system
                                                   ! [W m-2(bld)]
REAL, DIMENSION(:), INTENT(IN)   :: PQIN          ! Internal heat gains [W m-2(floor)]
REAL, DIMENSION(:), INTENT(IN)   :: PQIN_FRAD     ! Radiant fraction of internal heat gains
REAL, DIMENSION(:), INTENT(IN)   :: PQIN_FLAT     ! Latent franction of internal heat gains
REAL, DIMENSION(:), INTENT(IN)   :: PGR           ! Glazing ratio
REAL, DIMENSION(:), INTENT(IN)   :: PEFF_HEAT     ! Efficiency of the heating system
REAL, DIMENSION(:), INTENT(IN)   :: PINF          ! Infiltration flow rate [AC/H]
REAL, DIMENSION(:), INTENT(IN)   :: PTCOOL_TARGET ! Cooling setpoint of HVAC system [K]
REAL, DIMENSION(:), INTENT(IN)   :: PTHEAT_TARGET ! Heating setpoint of HVAC system [K]
REAL, DIMENSION(:), INTENT(IN)   :: PHR_TARGET    ! Relative humidity setpoint
REAL, DIMENSION(:), INTENT(INOUT):: PT_WIN2       ! Indoor window temperature [K]
REAL, DIMENSION(:), INTENT(INOUT):: PQI_BLD       ! Indoor air specific humidity [kg kg-1]
REAL, DIMENSION(:), INTENT(IN)   :: PV_VENT       ! Ventilation flow rate [AC/H]
REAL, DIMENSION(:), INTENT(IN)   :: PCAP_SYS_HEAT ! Capacity of the heating system 
                                                  ! [W m-2(bld)]
REAL, DIMENSION(:), INTENT(INOUT):: PCAP_SYS_RAT  ! Rated capacity of the cooling system
                                                  ! [W m-2(bld)]
REAL, DIMENSION(:), INTENT(IN)   :: PT_ADP        ! Apparatus dewpoint temperature of the
                                                  ! cooling coil [K]
REAL, DIMENSION(:), INTENT(INOUT):: PM_SYS_RAT    ! Rated HVAC mass flow rate 
                                                  ! [kg s-1 m-2(bld)]
REAL, DIMENSION(:), INTENT(IN)   :: PCOP_RAT      ! Rated COP of the cooling system
REAL, DIMENSION(:), INTENT(OUT)  :: PCAP_SYS      ! Actual capacity of the cooling system
                                                  ! [W m-2(bld)] 
REAL, DIMENSION(:), INTENT(OUT)  :: PM_SYS        ! Actual HVAC mass flow rate 
                                                  ! [kg s-1 m-2(bld)]
REAL, DIMENSION(:), INTENT(OUT)  :: PCOP          ! COP of the cooling system
REAL, DIMENSION(:), INTENT(OUT)  :: PQ_SYS        ! Supply air specific humidity [kg kg-1]
REAL, DIMENSION(:), INTENT(OUT)  :: PT_SYS        ! Supply air temperature [K]
REAL, DIMENSION(:), INTENT(IN)   :: PTR_SW_WIN    ! Solar radiation transmitted throught
                                                  ! windows [W m-2(bld)]
REAL, DIMENSION(:), INTENT(OUT)  :: PFAN_POWER    ! HVAC fan power
REAL, DIMENSION(:,:),INTENT(IN)  :: PHC_FLOOR     ! heat capacity for road layers
REAL, DIMENSION(:,:),INTENT(IN)  :: PTC_FLOOR     ! thermal conductivity for 
                                                  ! road layers
REAL, DIMENSION(:,:),INTENT(IN)  :: PD_FLOOR      ! depth of road layers
REAL, DIMENSION(:), INTENT(IN)   :: PLW_WA_TO_WIN ! Radiative heat trasfer coeff wall-window 
                                                  ! [W K-1 m-2] 
REAL, DIMENSION(:), INTENT(IN)   :: PLW_WB_TO_WIN ! Radiative heat trasfer coeff wall-window 
                                                  ! [W K-1 m-2] 
REAL, DIMENSION(:), INTENT(IN)   :: PLW_G_TO_WIN  ! Radiative heat trasfer coeff garden-window 
                                                  ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(IN)   :: PLW_R_TO_WIN  ! Radiative heat trasfer coeff road-window 
                                                  ! [W K-1 m-2] 
REAL, DIMENSION(:), INTENT(IN)   :: PLW_S_TO_WIN ! Radiative heat trasfer coeff window-sky 
                                                 ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(IN)   :: PLW_WIN_TO_WA! Radiative heat trasfer coeff window-wall
                                                 ! [W K-1 m-2] 
REAL, DIMENSION(:), INTENT(IN)   :: PLW_WIN_TO_WB! Radiative heat trasfer coeff window-wall
                                                 ! [W K-1 m-2] 
REAL, DIMENSION(:), INTENT(IN)   :: PLW_WIN_TO_R ! Radiative heat trasfer coeff window-road 
                                                 ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(IN)   :: PLW_NR_TO_WIN! Radiative heat trasfer coeff road(snow)-win 
                                                 ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(IN)   :: PLW_WIN_TO_NR! Radiative heat trasfer coeff win-road(snow) 
                                                 ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(INOUT):: PT_WIN1       ! outdoor window temperature [K]
REAL, DIMENSION(:), INTENT(IN)   :: PABS_SW_WIN   ! window absorbed shortwave radiation [W m-2] 
REAL, DIMENSION(:), INTENT(OUT)  :: PABS_LW_WIN   ! absorbed infrared rad. [W m-2]
REAL, DIMENSION(:), INTENT(IN)   :: PUGG_WIN      ! window glass-to-glass U-factro [W m-2 K-1]
 !new argument for PET calculation
REAL, DIMENSION(:), INTENT(OUT) :: PEMIT_LW_ROAD ! LW fluxes emitted by road (W/m2 surf road)
REAL, DIMENSION(:), INTENT(OUT) :: PEMIT_LW_FAC  ! LW fluxes emitted by wall (W/m2 surf wall)
REAL, DIMENSION(:), INTENT(OUT) :: PT_RAD_IND    ! Indoor mean radiant temperature [K]
REAL, DIMENSION(:), INTENT(OUT) :: PHU_BLD       ! Indoor relative humidity 0 < (-) < 1
REAL,                INTENT(IN)  :: PTIME        ! current time since midnight (UTC, s)
REAL, DIMENSION(:), INTENT(IN)  :: PE_SHADING    !energy not ref., nor absorbed, nor
                                                 !trans. by glazing [Wm-2(win)]
LOGICAL, DIMENSION(:),INTENT(INOUT) :: ONATVENT_NIGHT ! has natural ventilation
                                                 !been necessary/possible during the night
!
REAL, DIMENSION(:)  , INTENT(IN)  :: PN_FLOOR        ! Number of floors     
REAL, DIMENSION(:)  , INTENT(IN)  :: PWALL_O_BLD     ! Wall area [m2_wall/m2_bld]
REAL, DIMENSION(:)  , INTENT(IN)  :: PGLAZ_O_BLD     ! Window area [m2_win/m2_bld]
REAL, DIMENSION(:)  , INTENT(IN)  :: PMASS_O_BLD     ! Mass area [m2_mass/m2_bld]
REAL, DIMENSION(:)  , INTENT(IN)  :: PFLOOR_HW_RATIO ! H/W ratio of 1 floor level
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_FLOOR_MASS   ! View factor floor-mass
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_FLOOR_WALL   ! View factor floor-wall
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_FLOOR_WIN    ! View factor floor-window
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_FLOOR_ROOF   ! View factor floor-roof
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_WALL_FLOOR   ! View factor wall-floor
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_WALL_MASS    ! View factor wall-mass
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_WALL_WIN     ! View factor wall-win
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_WIN_FLOOR    ! View factor win-floor
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_WIN_MASS     ! View factor win-mass
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_WIN_WALL     ! View factor win-wall
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_MASS_FLOOR   ! View factor mass-floor
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_MASS_WALL    ! View factor mass-wall
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_MASS_WIN     ! View factor mass-window
LOGICAL,              INTENT(IN)  :: OCANOPY         ! is canopy active ?
 CHARACTER(LEN=5), INTENT(IN)      :: HCH_BEM         ! TEB option for outside building coef
                                                     ! "DOE-2"
REAL, DIMENSION(:)  , INTENT(IN)  :: PROUGH_ROOF     ! roof roughness coef
REAL, DIMENSION(:)  , INTENT(IN)  :: PROUGH_WALL     ! wall roughness coef
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_WIN_WIN      ! indoor win to win view factor
!
!* arguments for road watering
!
LOGICAL,              INTENT(IN)  :: OPAR_RD_IRRIG  ! True if there is potentially irrigation
REAL, DIMENSION(:),   INTENT(IN)  :: PRD_START_MONTH! First month of irrigation
REAL, DIMENSION(:),   INTENT(IN)  :: PRD_END_MONTH  ! Last month of irrigation
REAL, DIMENSION(:),   INTENT(IN)  :: PRD_START_HOUR ! First hour of irrigation (included)
REAL, DIMENSION(:),   INTENT(IN)  :: PRD_END_HOUR   ! Last hour of irrigation  (excluded)
REAL, DIMENSION(:),   INTENT(IN)  :: PRD_24H_IRRIG  ! diurnal averaged irrigation flux (kg/m2/s)
REAL, DIMENSION(:),   INTENT(OUT) :: PIRRIG_ROAD    ! man-made watering for road (kg/m2/s)
!
!*      0.2    Declarations of local variables
!
REAL, DIMENSION(SIZE(PTA)) :: ZVMOD          ! wind
REAL, DIMENSION(SIZE(PTA)) :: ZWS_ROOF_MAX   ! maximum deepness of roof
REAL, DIMENSION(SIZE(PTA)) :: ZWS_ROAD_MAX   ! and road water reservoirs
!
REAL, DIMENSION(SIZE(PTA)) :: ZAC_BLD        ! surface conductance inside the building itself in DEF building model
REAL, DIMENSION(SIZE(PTA)) :: ZTA            ! air temperature extrapolated at roof level
REAL, DIMENSION(SIZE(PTA)) :: ZQA            ! air humidity extrapolated at roof level
!
REAL, DIMENSION(SIZE(PTA)) :: ZDQS_ROAD      ! heat storage inside road
REAL, DIMENSION(SIZE(PTA)) :: ZDQS_ROOF      ! heat storage inside roof
REAL, DIMENSION(SIZE(PTA)) :: ZDQS_WALL_A    ! heat storage inside wall
REAL, DIMENSION(SIZE(PTA)) :: ZDQS_WALL_B    ! heat storage inside wall
REAL, DIMENSION(SIZE(PTA)) :: ZFLX_BLD_ROOF  !heat flux from inside through roof
REAL, DIMENSION(SIZE(PTA)) :: ZFLX_BLD_WALL_A!heat flux from inside through wall
REAL, DIMENSION(SIZE(PTA)) :: ZFLX_BLD_WALL_B!heat flux from inside through wall
REAL, DIMENSION(SIZE(PTA)) :: ZFLX_BLD_FLOOR !heat flux from inside through floor
REAL, DIMENSION(SIZE(PTA)) :: ZFLX_BLD_MASS  !heat flux from inside through mass
!
REAL, DIMENSION(SIZE(PTA)) :: ZDQS_SNOW_ROOF ! heat storage inside roof snowpack
REAL, DIMENSION(SIZE(PTA)) :: ZDQS_SNOW_ROAD ! heat storage inside road snowpack
REAL, DIMENSION(SIZE(PTA)) :: ZMELT_BLT      ! Snow melt for built & impervious part
!
! coefficients for LW computations over snow (from previous time-step)
!
REAL, DIMENSION(SIZE(PTA)) :: ZTSSNOW_ROAD ! road snow temperature
!                                          ! at previous time-step
! new local variables after BEM
!
REAL, DIMENSION(SIZE(PTA)) :: ZIMB_ROOF      ! residual energy imbalance
                                             ! of the roof for
                                             ! verification
REAL, DIMENSION(SIZE(PTA)) :: ZIMB_ROAD      ! road residual energy imbalance 
                                             ! for verification [W m-2]
REAL, DIMENSION(SIZE(PTA)) :: ZIMB_WALL      ! wall residual energy imbalance 
                                             ! for verification [W m-2]
REAL, DIMENSION(SIZE(PTA)) :: ZTS_ROAD       ! road surface temperature 
!                                            ! at previous time-step
REAL, DIMENSION(SIZE(PTA)) :: ZTS_WALL_A     ! wall A surface temperature 
!                                            ! at previous time-step
REAL, DIMENSION(SIZE(PTA)) :: ZTS_WALL_B     ! wall B surface temperature 
!                                            ! at previous time-step
REAL, DIMENSION(SIZE(PTA)) :: ZTS_WALL       ! averaged wall surface temperature 
!                                            ! at previous time-step
REAL, DIMENSION(SIZE(PTA)) :: ZTS_ROOF       ! roof surface temperature 
!                                            ! at previous time-step
REAL, DIMENSION(SIZE(PTA),SIZE(PT_WALL_A,2)) :: ZT_WALL ! averaged wall surface temperature 
!
INTEGER :: IWALL, IROOF                      ! number of wall, roof layer
REAL, DIMENSION(SIZE(PTA)) :: ZRADHT_IN     ! Indoor radiant heat transfer coefficient
                                                    ! [W K-1 m-2]
REAL, DIMENSION(SIZE(PTA)) :: ZTS_FLOOR       ! floor surface temperature [K]
REAL, DIMENSION(SIZE(PTA)) :: ZRAD_ROOF_WALL  ! rad. flux from roof to averaged wall [W m-2(roof)]
REAL, DIMENSION(SIZE(PTA)) :: ZRAD_ROOF_WIN   ! rad. flux from roof to window [W m-2(roof)]
REAL, DIMENSION(SIZE(PTA)) :: ZRAD_ROOF_FLOOR ! rad. flux from roof to floor [W m-2(roof)]
REAL, DIMENSION(SIZE(PTA)) :: ZRAD_ROOF_MASS  ! rad. flux from roof to mass [W m-2(roof)]
REAL, DIMENSION(SIZE(PTA)) :: ZCONV_ROOF_BLD  ! rad. flux from roof to bld [W m-2(roof)]
REAL, DIMENSION(SIZE(PTA)) :: ZRAD_WALL_FLOOR ! rad. flux from averaged wall to floor [W m-2(wall)]
REAL, DIMENSION(SIZE(PTA)) :: ZRAD_WALL_MASS  ! rad. flux from averaged wall to mass [W m-2(wall)]
REAL, DIMENSION(SIZE(PTA)) :: ZRAD_WIN_FLOOR  ! rad. flux from averaged wall to floor [W m-2(win)]
REAL, DIMENSION(SIZE(PTA)) :: ZRAD_WIN_MASS   ! rad. flux from averaged wall to mass [W m-2(win)]
REAL, DIMENSION(SIZE(PTA)) :: ZCONV_WALL_BLD  ! rad. flux from roof to bld [W m-2(wall)]
REAL, DIMENSION(SIZE(PTA)) :: ZCONV_WIN_BLD   ! rad. flux from roof to bld [W m-2(win)]
REAL, DIMENSION(SIZE(PTA)) :: ZAC_WIN         ! window aerodynamic conductance

REAL, DIMENSION(SIZE(PTA)) :: ZLOAD_IN_ROOF   ! indoor load on roof W/m2[roof]
REAL, DIMENSION(SIZE(PTA)) :: ZLOAD_IN_FLOOR   ! indoor load on floor W/m2[floor]
REAL, DIMENSION(SIZE(PTA)) :: ZLOAD_IN_WALL   ! indoor load on wall W/m2[wall]
REAL, DIMENSION(SIZE(PTA)) :: ZLOAD_IN_WIN   ! indoor load on win W/m2[win]
REAL, DIMENSION(SIZE(PTA)) :: ZLOAD_IN_MASS   ! indoor load on mass W/m2[mass]
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE                                             
!-------------------------------------------------------------------------------
!
!*      1.     Initializations
!              ---------------
!
!*      1.1    Water reservoirs
!              ----------------
!
IF (LHOOK) CALL DR_HOOK('TEB',0,ZHOOK_HANDLE)
ZWS_ROOF_MAX =  1. ! (1mm) maximum deepness of roof water reservoir
ZWS_ROAD_MAX =  1. ! (1mm) maximum deepness of road water reservoir
!
!*      1.2    radiative snow variables at previous time-step
!              ----------------------------------------------
!
ZTSSNOW_ROAD(:)=PTSSNOW_ROAD(:)
!
!
!*      1.3    indoor aerodynamique conductance for DEF case
!              ----------------------------------------------
!
ZAC_BLD(:) = XUNDEF
IF (HBEM=='DEF') ZAC_BLD=1. / 0.123 / (XCPD * PRHOA(:)) !* (normalized by rho Cp for convenience)
!-------------------------------------------------------------------------------
!
!*      1.3    number of roof/wall layer
!              -------------------------
!
IWALL = SIZE(PT_WALL_A,2)
IROOF = SIZE(PT_ROOF,2)
!
ZTS_WALL_A  (:)=PT_WALL_A   (:,1)
ZTS_WALL_B  (:)=PT_WALL_B   (:,1)
ZTS_WALL    (:)=0.5 * (ZTS_WALL_A(:)+ZTS_WALL_B(:))
ZTS_ROAD    (:)=PT_ROAD     (:,1)
ZTS_ROOF    (:)=PT_ROOF     (:,1)
!
!
!*      1.4    load on indoor walls
!              -------------------------
!
IF (HBEM=='BEM') THEN
  !
  ZLOAD_IN_ROOF = PF_FLOOR_WIN * PTR_SW_WIN + PQIN * PN_FLOOR * (1-PQIN_FLAT) * PQIN_FRAD  &
           / (2 + PWALL_O_BLD + PGLAZ_O_BLD + PMASS_O_BLD ) ! W/m2 [ROOF]
  ZLOAD_IN_FLOOR = PF_FLOOR_WIN * PTR_SW_WIN + PQIN * PN_FLOOR * (1-PQIN_FLAT) * PQIN_FRAD  &
           / (2 + PWALL_O_BLD + PGLAZ_O_BLD + PMASS_O_BLD )
  ZLOAD_IN_MASS = PF_MASS_WIN * PTR_SW_WIN + PQIN * PN_FLOOR * (1-PQIN_FLAT) * PQIN_FRAD  &
           / (2 + PWALL_O_BLD + PGLAZ_O_BLD + PMASS_O_BLD )
  ZLOAD_IN_WALL = PF_WALL_WIN * PTR_SW_WIN + PQIN * PN_FLOOR * (1-PQIN_FLAT) * PQIN_FRAD  &
           / (2 + PWALL_O_BLD + PGLAZ_O_BLD + PMASS_O_BLD )
  ZLOAD_IN_WIN = PF_WIN_WIN * PTR_SW_WIN + PQIN * PN_FLOOR * (1-PQIN_FLAT) * PQIN_FRAD  &
           / (2 + PWALL_O_BLD + PGLAZ_O_BLD + PMASS_O_BLD )
ELSE
  ZLOAD_IN_ROOF = 0.
  ZLOAD_IN_FLOOR = 0.
  ZLOAD_IN_MASS = 0.
  ZLOAD_IN_WALL = 0.
  ZLOAD_IN_WIN = 0.
ENDIF
!-------------------------------------------------------------------------------
!
!*      2.     Snow-covered surfaces relative effects
!              --------------------------------------
!
!*      2.1    Effects on water reservoirs
!              ---------------------------
!
ZWS_ROOF_MAX(:) = ZWS_ROOF_MAX(:) * PDF_ROOF(:)
ZWS_ROAD_MAX(:) = ZWS_ROAD_MAX(:) * PDF_ROAD(:)
!
!-------------------------------------------------------------------------------
!
!*      3.     Surface drag
!              ------------
!
 CALL URBAN_DRAG(HZ0H, HIMPLICIT_WIND, PTSTEP, PT_CANYON, PQ_CANYON,  &
                PU_CANYON, PT_LOWCAN, PQ_LOWCAN, PU_LOWCAN,          &
                PZ_LOWCAN, ZTS_ROOF, ZTS_ROAD, ZTS_WALL,             &
                PTS_GARDEN, PDN_ROOF, PDN_ROAD,                      &
                PEXNS, PEXNA, PTA, PQA, PPS, PRHOA,                  &
                PZREF, PUREF, PVMOD,                                 &
                PZ0_TOWN,                                            &
                PBLD, PGARDEN, PROAD,                                &
                PBLD_HEIGHT, PCAN_HW_RATIO,                          &
                PWALL_O_GRND,                                        &
                PWS_ROOF, PWS_ROAD,                                  &
                ZWS_ROOF_MAX, ZWS_ROAD_MAX,                          &
                PPEW_A_COEF, PPEW_B_COEF,                            &
                PPEW_A_COEF_LOWCAN, PPEW_B_COEF_LOWCAN,              &
                PQSAT_ROOF, PQSAT_ROAD, PDELT_ROOF, PDELT_ROAD,      &
                PCD, PCDN, PAC_ROOF, PAC_ROOF_WAT,                   &
                PAC_WALL, PAC_ROAD, PAC_ROAD_WAT, PAC_TOP,           &
                PAC_GARDEN, PRI_TOWN,                                &
                PUW_ROAD, PUW_ROOF, PDUWDU_ROAD, PDUWDU_ROOF,        &
                PUSTAR_TOWN, OCANOPY, PT_WIN1,  ZAC_WIN, HCH_BEM,    &
                PROUGH_ROOF, PROUGH_WALL                             )
!
!* area-averaged heat transfer coefficient
!
ZVMOD(:) = WIND_THRESHOLD(PVMOD(:),PUREF(:))
!
PCH_TOWN(:) = (PBLD(:) * PAC_ROOF(:) + (1.-PBLD(:)) * PAC_TOP (:)) / ZVMOD(:)
!
!* aggregation of momentum fluxes for roofs (=> derivate of flux also recalculated)
!
PUW_ROOF (:) = (1-PFRAC_GR(:)) * PUW_ROOF(:) + PFRAC_GR(:) * PUW_GREENROOF(:)
WHERE (PVMOD(:)/=0.) PDUWDU_ROOF(:) = 2. * PUW_ROOF(:) / PVMOD(:)
!
!-------------------------------------------------------------------------------
!
!*      4.     Extrapolation of atmospheric T and q at roof level (for fluxes computation)
!              --------------------------------------------------
!
ZTA(:) = PTA(:) * PEXNS(:) / PEXNA(:)
ZQA(:) = PQA(:) * QSAT(PTA(:),PPS(:)) / QSAT(ZTA(:),PPA(:))
!
!-------------------------------------------------------------------------------
!
!*      5.     Snow mantel model
!              -----------------
!
 CALL URBAN_SNOW_EVOL(PT_LOWCAN, PQ_LOWCAN, PU_LOWCAN,                         &
                     ZTS_ROOF,ZTS_ROAD,ZTS_WALL_A, ZTS_WALL_B,                &
                     PT_ROOF, PD_ROOF, PTC_ROOF, PHC_ROOF,                    &
                     HSNOW_ROOF,                                              &
                     PWSNOW_ROOF, PTSNOW_ROOF, PRSNOW_ROOF, PASNOW_ROOF,      &
                     PTSSNOW_ROOF, PESNOW_ROOF,                               &
                     HSNOW_ROAD,                                              &
                     PWSNOW_ROAD, PTSNOW_ROAD, PRSNOW_ROAD, PASNOW_ROAD,      &
                     PTSSNOW_ROAD, PESNOW_ROAD,                               &
                     PPS, ZTA, ZQA, PRHOA,                                    &
                     PLW_RAD,                                                 &
                     PSR, PZREF, PUREF, PVMOD,                                &
                     PTSTEP,                                                  &
                     PZ_LOWCAN,                                               &
                     PDN_ROOF, PABS_SW_SNOW_ROOF, PABS_LW_SNOW_ROOF,          &
                     PDN_ROAD, PABS_SW_SNOW_ROAD, PABS_LW_SNOW_ROAD,          &
                     PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,    &
                     PMELT_ROOF,                                              &
                     PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,    &
                     PMELT_ROAD,                                              &
                     PLW_WA_TO_NR, PLW_WB_TO_NR, PLW_S_TO_NR, PLW_WIN_TO_NR,  &
                     ZDQS_SNOW_ROOF, ZDQS_SNOW_ROAD, PT_WIN1                  )
!
!-------------------------------------------------------------------------------
!
!*      6.    LW properties
!              -------------
!
PDF_ROAD (:) = 1. - PDN_ROAD (:)
!
!-------------------------------------------------------------------------------
!
!*      7.    Indoor radiative temperature
!              ---------------------------
!
! uses the averaged temperature of both walls for the building energy balance
ZT_WALL   (:,:)=0.5 * (PT_WALL_A(:,:)+PT_WALL_B(:,:))
!
SELECT CASE(HBEM)
   CASE("DEF")
      ZTS_FLOOR(:) = 19. + XTT
      PT_RAD_IND(:) = (PWALL_O_HOR(:) / PBLD(:) * ZT_WALL(:,IWALL) + &
                    PT_ROOF(:,IROOF) + ZTS_FLOOR(:)) &
                   / (PWALL_O_HOR(:) / PBLD(:) + 1. + 1.) 
      ZRADHT_IN(:) = XUNDEF
   CASE("BEM")
      ZTS_FLOOR(:) = PT_FLOOR(:,1)
      PT_RAD_IND(:)  = (PT_MASS(:,1)*PMASS_O_BLD(:) + ZT_WALL(:,IWALL)*PWALL_O_BLD(:)     &
                 + ZTS_FLOOR(:) + PT_ROOF(:,IROOF) + PT_WIN2(:) * PGLAZ_O_BLD(:)) &
                 /(PMASS_O_BLD(:) + PWALL_O_BLD(:) + 1. + 1. + PGLAZ_O_BLD(:))
      !             Assuming indoor surface emissivities of 0.9
      ZRADHT_IN(:)   = 0.9 * 0.9 * 4 * XSTEFAN * PT_RAD_IND(:)**3          
END SELECT
!
!
!*      7.    Roof Ts computation
!              -------------------
!
!* ts_roof and qsat_roof are updated
!

 CALL ROOF_LAYER_E_BUDGET(PT_ROOF, PQSAT_ROOF, PTI_BLD, ZAC_BLD, PTSTEP, &
                         HBEM, PHC_ROOF, PTC_ROOF, PD_ROOF, PDN_ROOF,   &
                         PRHOA, PAC_ROOF, PAC_ROOF_WAT, PLW_RAD, PPS,   &
                         PDELT_ROOF, ZTA, ZQA, PEXNA, PEXNS,            &
                         PABS_SW_ROOF, PGSNOW_ROOF, PEMIS_ROOF,         &
                         ZFLX_BLD_ROOF, ZDQS_ROOF, PABS_LW_ROOF,        &
                         PH_ROOF, PLEW_ROOF, ZIMB_ROOF,                 &
                         PFRAC_GR, PG_GREENROOF_ROOF,                   &
                         PF_FLOOR_MASS, PF_FLOOR_WALL, PF_FLOOR_WIN,    &
                         PF_FLOOR_ROOF, ZRADHT_IN,         &
                         PT_MASS(:,1), PT_WIN2, ZTS_FLOOR(:),           &
                         ZT_WALL(:,IWALL), ZRAD_ROOF_WALL,              &
                         ZRAD_ROOF_WIN, ZRAD_ROOF_FLOOR, ZRAD_ROOF_MASS,&
                         ZCONV_ROOF_BLD,                                &
                         PRR, & !modif to add heating/cooling of rain
                         ZLOAD_IN_ROOF )
!
!-------------------------------------------------------------------------------
!
!*      8.    Road Ts computations
!              -----------------------------
!
!* Road watering

 CALL TEB_IRRIG(OPAR_RD_IRRIG, PTSTEP, TPTIME%TDATE%MONTH, PTSUN, &
               PRD_START_MONTH, PRD_END_MONTH, PRD_START_HOUR,   &
               PRD_END_HOUR, PRD_24H_IRRIG, PIRRIG_ROAD          )

!* ts_road, ts_wall, qsat_road, t_canyon and q_canyon are updated
!
 CALL ROAD_LAYER_E_BUDGET(PT_ROAD, PTSTEP, PHC_ROAD, PTC_ROAD, PD_ROAD,       &
                         PDN_ROAD, PRHOA, PAC_ROAD, PAC_ROAD_WAT,            &
                         PLW_RAD, PPS, PQSAT_ROAD, PDELT_ROAD, PEXNS,        &
                         PABS_SW_ROAD, PGSNOW_ROAD, PQ_LOWCAN, PT_LOWCAN,    &
                         ZTS_WALL_A, ZTS_WALL_B, ZTSSNOW_ROAD,               &
                         PTS_GARDEN, PT_WIN1,                                &
                         PLW_WA_TO_R, PLW_WB_TO_R, PLW_S_TO_R, PLW_WIN_TO_R, &
                         PEMIS_ROAD, PEMIT_LW_ROAD, ZDQS_ROAD, PABS_LW_ROAD, &
                         PH_ROAD, PLEW_ROAD, ZIMB_ROAD, PRR+PIRRIG_ROAD      )
!
!-------------------------------------------------------------------------------
!
!*      8.     Wall Ts computations
!              -----------------------------
!
 CALL FACADE_E_BUDGET(HWALL_OPT, HBEM,                                    &
                     PT_WALL_A, PT_WALL_B, PTSTEP,                       &
                     PHC_WALL, PTC_WALL, PD_WALL,                        &
                     PDN_ROAD, PRHOA, PAC_WALL,                          &
                     ZAC_BLD, PTI_BLD, PLW_RAD, PPS, PEXNS,              &
                     PABS_SW_WALL_A, PABS_SW_WALL_B, PT_CANYON,          &
                     ZTS_ROAD, ZTSSNOW_ROAD, PTS_GARDEN, PUGG_WIN, PGR,  &
                     PABS_SW_WIN, PF_WIN_MASS, PF_WIN_WALL,              &
                     PF_WIN_FLOOR, PT_MASS(:,1), ZTS_FLOOR(:),           &
                     PEMIS_WALL,                                         &
                     PLW_WA_TO_WB,                                       &
                     PLW_R_TO_WA, PLW_R_TO_WB, PLW_G_TO_WA, PLW_G_TO_WB, &
                     PLW_S_TO_WA, PLW_S_TO_WB, PLW_NR_TO_WA,PLW_NR_TO_WB,&
                     PLW_WIN_TO_WA, PLW_WIN_TO_WB,                       &
                     PLW_S_TO_WIN, PLW_WA_TO_WIN, PLW_WB_TO_WIN,         &
                     PLW_R_TO_WIN, PLW_G_TO_WIN, PLW_NR_TO_WIN, PT_WIN1, &
                     ZFLX_BLD_WALL_A, ZDQS_WALL_A,                       &
                     ZFLX_BLD_WALL_B, ZDQS_WALL_B,                       &
                     PT_WIN2, PEMIT_LW_FAC,                              &
                     PABS_LW_WALL_A, PABS_LW_WALL_B, PABS_LW_WIN,        &
                     PH_WALL_A, PH_WALL_B,                               &
                     ZIMB_WALL, PF_WALL_MASS, PF_WALL_FLOOR,             &
                     PF_WALL_WIN, ZRADHT_IN, ZRAD_ROOF_WALL,             &
                     ZRAD_ROOF_WIN, ZRAD_WALL_FLOOR, ZRAD_WALL_MASS,     &
                     ZRAD_WIN_FLOOR, ZRAD_WIN_MASS, ZCONV_WALL_BLD,      &
                     ZCONV_WIN_BLD, ZAC_WIN,                             &
                     ZLOAD_IN_WALL, ZLOAD_IN_WIN                         )
!
!-------------------------------------------------------------------------------
!
!*      9.     Evolution of interior building air temperature
!              ----------------------------------------------
!
! uses the averaged temperature of both walls for the building energy balance
ZT_WALL   (:,:)=0.5 * (PT_WALL_A(:,:)+PT_WALL_B(:,:))
!
SELECT CASE(HBEM)
 CASE("DEF")
!
   CALL BLD_E_BUDGET(.TRUE., PTSTEP, PBLD, PWALL_O_HOR,        &
                     PRHOA, PT_ROOF, ZT_WALL, PTI_BLD, ZTS_FLOOR(:) )

   !variables that needs to be computed apart
   PQI_BLD = 0.5 * QSAT(PTI_BLD, PPS)
   !variables that need to be set 0 for calculation
   PH_WASTE      (:) = 0.
   PLE_WASTE     (:) = 0.
   ZFLX_BLD_FLOOR(:) = 0.
   ZFLX_BLD_MASS (:) = 0.
   !other variables
   PH_BLD_COOL(:) = XUNDEF
   PLE_BLD_COOL(:)= XUNDEF
   PT_BLD_COOL(:) = XUNDEF
   PH_BLD_HEAT(:) = XUNDEF
   PLE_BLD_HEAT(:)= XUNDEF
   PHVAC_COOL(:)  = XUNDEF
   PHVAC_HEAT(:)  = XUNDEF
   PCAP_SYS(:)    = XUNDEF
   PM_SYS(:)      = XUNDEF
   PCOP(:)        = XUNDEF
   PQ_SYS(:)      = XUNDEF
   PFAN_POWER(:)  = XUNDEF
   PHU_BLD(:)     = XUNDEF

 CASE("BEM")
  CALL BEM(PTSTEP, PTSUN, HCOOL_COIL, HHEAT_COIL,                 &
    OAUTOSIZE, KDAY, HNATVENT,                                    &
    PPS, PRHOA, PT_CANYON, PQ_CANYON, PU_CANYON,                  &
    PT_ROOF, PT_WALL_A, PT_WALL_B, PBLD, PBLD_HEIGHT, PWALL_O_HOR,&
    PHC_FLOOR, PTC_FLOOR, PD_FLOOR,                               &
    PT_WIN2, PGR, PQIN, PQIN_FRAD, PQIN_FLAT,                     &
    PEFF_HEAT, PINF, PTCOOL_TARGET, PTHEAT_TARGET, PHR_TARGET,    &
    PF_WATER_COND, PV_VENT, PCAP_SYS_HEAT, PT_ADP, PCOP_RAT,      &
    PAUX_MAX, PT_FLOOR, PT_MASS, PTI_BLD, PQI_BLD,                &
    PCAP_SYS_RAT, PM_SYS_RAT, ONATVENT_NIGHT,                     &
    PH_BLD_COOL, PH_BLD_HEAT, PLE_BLD_COOL, PLE_BLD_HEAT,         &
    PT_BLD_COOL, PHVAC_COOL, PT_SYS, PQ_SYS, PH_WASTE, PLE_WASTE, & 
    PFAN_POWER, PHVAC_HEAT, PM_SYS, PCOP, PCAP_SYS, PHU_BLD,      &
    PTR_SW_WIN, PT_RAD_IND, ZFLX_BLD_FLOOR, ZFLX_BLD_MASS,        &
    PF_FLOOR_MASS, PF_FLOOR_WALL, PF_FLOOR_WIN,                   &
    ZRADHT_IN, PN_FLOOR, PWALL_O_BLD, PGLAZ_O_BLD, PMASS_O_BLD,   &
    PFLOOR_HW_RATIO, PF_FLOOR_ROOF, PF_MASS_FLOOR, PF_MASS_WALL,  &
    PF_MASS_WIN, ZRAD_ROOF_MASS, ZRAD_ROOF_FLOOR, ZRAD_WALL_MASS, &
    ZRAD_WALL_FLOOR, ZRAD_WIN_MASS, ZRAD_WIN_FLOOR,               &
    ZCONV_ROOF_BLD, ZCONV_WALL_BLD, ZCONV_WIN_BLD,                &
    ZLOAD_IN_FLOOR, ZLOAD_IN_MASS                                 )

   PH_WASTE  = PH_WASTE  * PBLD
   PLE_WASTE = PLE_WASTE * PBLD
END SELECT
!
!-------------------------------------------------------------------------------
!
!*      10.    Fluxes over built surfaces
!              --------------------------
!
 CALL URBAN_FLUXES   (HIMPLICIT_WIND, OCANOPY, PT_CANYON,                      &
                     PPEW_A_COEF, PPEW_B_COEF,                                &
                     PEXNS, PRHOA, PVMOD,                                     &
                     PH_TRAFFIC, PLE_TRAFFIC, PH_INDUSTRY, PLE_INDUSTRY,      &
                     PBLD, PROAD, PWALL_O_HOR, PFRAC_GR,                      &
                     PABS_SW_ROOF, PABS_LW_ROOF,                              &
                     PABS_SW_WALL_A, PABS_LW_WALL_A,                          &
                     PABS_SW_WALL_B, PABS_LW_WALL_B,                          &
                     PABS_SW_ROAD, PABS_LW_ROAD,                              &
                     PAC_WALL, PCD,                                           &
                     PDF_ROOF, PDN_ROOF, PDF_ROAD, PDN_ROAD,                  &
                     PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,    &
                     PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,    &
                     PMELT_ROOF, ZDQS_ROOF,                                   &
                     PMELT_ROAD, ZDQS_ROAD,                                   &
                     ZDQS_WALL_A, ZDQS_WALL_B,                                &
                     ZFLX_BLD_ROOF, ZFLX_BLD_WALL_A, ZFLX_BLD_WALL_B,         &
                     ZFLX_BLD_FLOOR,                                          &
                     ZFLX_BLD_MASS, PE_SHADING, HBEM, PGR, PT_WIN1, PH_WASTE, &
                     PLE_WASTE, PF_WASTE_CAN, PABS_SW_WIN, PABS_LW_WIN, PQIN, &
                     PN_FLOOR, PLEW_ROOF, PLEW_ROAD,                          &
                     PRN_GREENROOF, PH_GREENROOF, PLE_GREENROOF,              &
                     PGFLUX_GREENROOF,                                        &
                     PH_ROOF, PH_ROAD, PH_WALL_A, PH_WALL_B,                  &
                     PRN_ROOF,PLE_ROOF,                                       &
                     PGFLUX_ROOF, PRN_ROAD, PLE_ROAD, PGFLUX_ROAD,            &
                     PRN_WALL_A, PLE_WALL_A, PGFLUX_WALL_A,                   &
                     PRN_WALL_B, PLE_WALL_B, PGFLUX_WALL_B,                   &
                     PRN_STRLROOF, PH_STRLROOF, PLE_STRLROOF,                 & 
                     PGFLUX_STRLROOF,                                         &
                     PRN_BLT, PH_BLT, PLE_BLT,                                &
                     PGFLUX_BLT, ZMELT_BLT, PQF_BLD, PFLX_BLD, PDQS_TOWN,     &
                     PQF_TOWN, PUSTAR_TOWN, PHVAC_COOL, PHVAC_HEAT            )
!
!
! Water transfer from snow reservoir to water reservoir in case of snow melt
!
WHERE (PMELT_ROOF(:) .GT. 0.)
 PWS_ROOF(:) = MIN(ZWS_ROOF_MAX,PWS_ROOF(:) + PMELT_ROOF(:)*PTSTEP)
ENDWHERE
!
WHERE (PMELT_ROAD(:) .GT. 0.)
 PWS_ROAD(:) = MIN(ZWS_ROAD_MAX,PWS_ROAD(:) + PMELT_ROAD(:)*PTSTEP)
ENDWHERE
!
!-------------------------------------------------------------------------------
!
!*      11.    Roof ans road reservoirs evolution
!              ----------------------------------
!
 CALL URBAN_HYDRO(ZWS_ROOF_MAX,ZWS_ROAD_MAX, PWS_ROOF, PWS_ROAD,      &
                 PRR, PIRRIG_ROAD, PTSTEP, PBLD, PLE_ROOF, PLE_ROAD,  &
                 PRUNOFF_STRLROOF,                                    &
                 PRUNOFF_ROAD                                         )
!
PRUNOFF_ROOF(:) =   ( 1. - PFRAC_GR(:) ) *   PRUNOFF_STRLROOF(:)                        &
                  +        PFRAC_GR(:)   * ( PRUNOFF_GREENROOF(:) + PDRAIN_GREENROOF(:) )
                                                           
!
!-------------------------------------------------------------------------------
!
!*      19.    Compute aerodynamical resistance 
!              --------------------------------
!
PRESA_TOWN(:) = 1. / ( PBLD(:) * PAC_ROOF(:)  + ( 1. - PBLD(:)) * PAC_TOP (:))
!
IF (LHOOK) CALL DR_HOOK('TEB',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE TEB
