!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!   ##########################################################################
    SUBROUTINE TEB  (TOP, T, BOP, B, TIR, DMT, HIMPLICIT_WIND, PTSUN,                   &
                     PT_CANYON, PQ_CANYON, PU_CANYON, PT_LOWCAN, PQ_LOWCAN, PU_LOWCAN,  &
                     PZ_LOWCAN, PPEW_A_COEF, PPEW_B_COEF, PPEW_A_COEF_LOWCAN,           &
                     PPEW_B_COEF_LOWCAN, PPS, PPA, PEXNS, PEXNA, PTA, PQA, PRHOA,       &
                     PLW_RAD, PRR, PSR, PZREF, PUREF, PVMOD, PH_TRAFFIC, PLE_TRAFFIC,   &
                     PTSTEP, PDF_RF, PDN_RF, PDF_RD, PDN_RD, PQSAT_RF, PQSAT_RD,        &
                     PDELT_RF, PDELT_RD, PTS_GARDEN, PLEW_RF, PUW_GR, PLEW_RD, PLE_WL_A,&
                     PLE_WL_B, PRNSN_RF, PHSN_RF, PLESN_RF, PGSN_RF, PMELT_RF, PRN_GR, &
                     PH_GR, PLE_GR, PGFLUX_GR, PDRAIN_GR, PRUNOFF_GR, PRNSN_RD,    &
                     PHSN_RD, PLESN_RD, PGSN_RD, PMELT_RD, PUW_RD, PUW_RF, PDUWDU_RD,   &
                     PDUWDU_RF, PUSTAR_TWN, PCD, PCDN, PCH_TWN, PRI_TWN, PRESA_TWN,     &
                     PAC_RF, PAC_RD, PAC_WL, PAC_TOP, PAC_GARDEN, PAC_RF_WAT,           &
                     PAC_RD_WAT, PLW_WA_TO_WB, PLW_WA_TO_R, PLW_WB_TO_R, PLW_WA_TO_NR,  &
                     PLW_WB_TO_NR, PLW_R_TO_WA, PLW_R_TO_WB, PLW_G_TO_WA, PLW_G_TO_WB,  &
                     PLW_S_TO_WA, PLW_S_TO_WB, PLW_S_TO_R, PLW_S_TO_NR, PLW_NR_TO_WA,   &
                     PLW_NR_TO_WB, PLW_NR_TO_WIN, PLW_WA_TO_WIN, PLW_WB_TO_WIN,         &
                     PLW_G_TO_WIN, PLW_R_TO_WIN, PLW_S_TO_WIN, PLW_WIN_TO_WA,           &
                     PLW_WIN_TO_WB, PLW_WIN_TO_R, PLW_WIN_TO_NR, KDAY, PEMIT_LW_FAC,    &
                     PEMIT_LW_RD, PT_RAD_IND, PHU_BLD, PTIME, PE_SHADING )
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
!     ------------
!
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_TEB_IRRIG_n, ONLY : TEB_IRRIG_t
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_t
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
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(BEM_OPTIONS_t), INTENT(INOUT) :: BOP
TYPE(BEM_t), INTENT(INOUT) :: B
TYPE(TEB_IRRIG_t), INTENT(INOUT) :: TIR
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: DMT
!
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
!                                                     ! 'OLD' = direct
!                                                     ! 'NEW' = Taylor serie, order 1
REAL, DIMENSION(:),   INTENT(IN)    :: PTSUN              ! solar time   (s from midnight)
REAL, DIMENSION(:), INTENT(INOUT) :: PT_CANYON     ! canyon air temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PQ_CANYON     ! canyon air specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PU_CANYON     ! canyon hor. wind
REAL, DIMENSION(:), INTENT(IN)    :: PU_LOWCAN     ! wind near the road
REAL, DIMENSION(:), INTENT(IN)    :: PT_LOWCAN     ! temp. near the road
REAL, DIMENSION(:), INTENT(IN)    :: PQ_LOWCAN     ! hum. near the road
REAL, DIMENSION(:), INTENT(IN)    :: PZ_LOWCAN     ! height of atm. var. near the road
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
REAL, DIMENSION(:), INTENT(IN)    :: PH_TRAFFIC    ! anthropogenic sensible
!                                                  ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PLE_TRAFFIC   ! anthropogenic latent
!                                                  ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PEXNA         ! exner function
                                                   ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA         ! air density
                                                   ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD       ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PRR           ! rain rate
REAL, DIMENSION(:), INTENT(IN)    :: PSR           ! snow rate
REAL, DIMENSION(:), INTENT(IN)    :: PZREF         ! reference height of the first
                                                   ! atmospheric level (temperature)
REAL, DIMENSION(:), INTENT(IN)    :: PUREF         ! reference height of the first
                                                   ! atmospheric level (wind)
REAL,               INTENT(IN)    :: PTSTEP        ! time step
!
REAL, DIMENSION(:), INTENT(INOUT) :: PDF_RF      ! snow-free    fraction on roofs
REAL, DIMENSION(:), INTENT(INOUT) :: PDN_RF      ! snow-covered fraction on roofs
REAL, DIMENSION(:), INTENT(INOUT) :: PDF_RD      ! snow-free    fraction on roads
REAL, DIMENSION(:), INTENT(INOUT) :: PDN_RD      ! snow-covered fraction on roads
REAL, DIMENSION(:), INTENT(OUT)   :: PQSAT_RF    ! hum at saturation over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PQSAT_RD    ! hum at saturation over road
REAL, DIMENSION(:), INTENT(OUT)   :: PDELT_RF    ! water fraction on roof
REAL, DIMENSION(:), INTENT(OUT)   :: PDELT_RD    ! water fraction on road
!
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GARDEN    ! GARDEN area surf temp.
! greenroof
REAL, DIMENSION(:), INTENT(OUT)   :: PLEW_RF    ! latent heat flux over roof (snow)
REAL, DIMENSION(:), INTENT(OUT)   :: PLEW_RD    ! latent heat flux over road (snow)
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WL_A   ! latent heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WL_B   ! latent heat flux over wall
!
REAL, DIMENSION(:), INTENT(IN)    :: PUW_GR     ! Momentum flux for greenroofs
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSN_RF ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSN_RF  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESN_RF ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSN_RF  ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_RF   ! snow melt
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSN_RD ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSN_RD  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESN_RD ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSN_RD  ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_RD   ! snow melt
!
REAL, DIMENSION(:), INTENT(OUT)   :: PUW_RD     ! Momentum flux for roads
REAL, DIMENSION(:), INTENT(OUT)   :: PUW_RF     ! Momentum flux for roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PDUWDU_RD  !
REAL, DIMENSION(:), INTENT(OUT)   :: PDUWDU_RF  !
REAL, DIMENSION(:), INTENT(OUT)   :: PUSTAR_TWN ! friciton velocity over town
!
REAL, DIMENSION(:), INTENT(IN)    :: PRN_GR     ! net radiation over greenroof
REAL, DIMENSION(:), INTENT(IN)    :: PH_GR      ! sensible heat flux over greenroof
REAL, DIMENSION(:), INTENT(IN)    :: PLE_GR     ! latent heat flux over greenroof
REAL, DIMENSION(:), INTENT(IN)    :: PGFLUX_GR  ! flux through the greenroof
REAL, DIMENSION(:), INTENT(IN)    :: PRUNOFF_GR ! runoff over green roofs
REAL, DIMENSION(:), INTENT(IN)    :: PDRAIN_GR  ! outlet drainage at base of green roofs
!
REAL, DIMENSION(:), INTENT(OUT)   :: PCD          ! town averaged drag coefficient
REAL, DIMENSION(:), INTENT(OUT)   :: PCDN         ! town averaged neutral drag coefficient
REAL, DIMENSION(:), INTENT(OUT)   :: PCH_TWN     ! town averaged heat transfer
!                                                 ! coefficient
REAL, DIMENSION(:), INTENT(OUT)   :: PRI_TWN      ! town averaged Richardson number
REAL, DIMENSION(:), INTENT(OUT)   :: PRESA_TWN    ! town aerodynamical resistance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_RF      ! roof conductance
REAL, DIMENSION(:), INTENT(INOUT) :: PAC_RD      ! road conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_WL      ! wall conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_TOP       ! top conductance
REAL, DIMENSION(:), INTENT(IN)    :: PAC_GARDEN    ! garden conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_RF_WAT  ! roof water conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_RD_WAT  ! roof water conductance
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
INTEGER,            INTENT(IN)     :: KDAY         ! Simulation day
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
 !new argument for PET calculation
REAL, DIMENSION(:), INTENT(OUT) :: PEMIT_LW_RD ! LW fluxes emitted by road (W/m2 surf road)
REAL, DIMENSION(:), INTENT(OUT) :: PEMIT_LW_FAC  ! LW fluxes emitted by wall (W/m2 surf wall)
REAL, DIMENSION(:), INTENT(OUT) :: PT_RAD_IND    ! Indoor mean radiant temperature [K]
REAL, DIMENSION(:), INTENT(OUT) :: PHU_BLD       ! Indoor relative humidity 0 < (-) < 1
REAL,                INTENT(IN)  :: PTIME        ! current time since midnight (UTC, s)
REAL, DIMENSION(:), INTENT(IN)  :: PE_SHADING    !energy not ref., nor absorbed, nor
                                                 !trans. by glazing [Wm-2(win)]
!
!*      0.2    Declarations of local variables
!
REAL, DIMENSION(SIZE(PTA)) :: ZVMOD          ! wind
REAL, DIMENSION(SIZE(PTA)) :: ZWS_RF_MAX   ! maximum deepness of roof
REAL, DIMENSION(SIZE(PTA)) :: ZWS_RD_MAX   ! and road water reservoirs
!
REAL, DIMENSION(SIZE(PTA)) :: ZAC_BLD        ! surface conductance inside the building itself in DEF building model
REAL, DIMENSION(SIZE(PTA)) :: ZTA            ! air temperature extrapolated at roof level
REAL, DIMENSION(SIZE(PTA)) :: ZQA            ! air humidity extrapolated at roof level
!
REAL, DIMENSION(SIZE(PTA)) :: ZDQS_RD      ! heat storage inside road
REAL, DIMENSION(SIZE(PTA)) :: ZDQS_RF      ! heat storage inside roof
REAL, DIMENSION(SIZE(PTA)) :: ZDQS_WL_A    ! heat storage inside wall
REAL, DIMENSION(SIZE(PTA)) :: ZDQS_WL_B    ! heat storage inside wall
REAL, DIMENSION(SIZE(PTA)) :: ZFLX_BLD_RF  !heat flux from inside through roof
REAL, DIMENSION(SIZE(PTA)) :: ZFLX_BLD_WL_A!heat flux from inside through wall
REAL, DIMENSION(SIZE(PTA)) :: ZFLX_BLD_WL_B!heat flux from inside through wall
REAL, DIMENSION(SIZE(PTA)) :: ZFLX_BLD_FL !heat flux from inside through floor
REAL, DIMENSION(SIZE(PTA)) :: ZFLX_BLD_MA  !heat flux from inside through mass
!
REAL, DIMENSION(SIZE(PTA)) :: ZDQS_SN_RF ! heat storage inside roof snowpack
REAL, DIMENSION(SIZE(PTA)) :: ZDQS_SN_RD ! heat storage inside road snowpack
REAL, DIMENSION(SIZE(PTA)) :: ZMELT_BLT      ! Snow melt for built & impervious part
!
! coefficients for LW computations over snow (from previous time-step)
!
REAL, DIMENSION(SIZE(PTA)) :: ZTSSN_RD ! road snow temperature
!                                          ! at previous time-step
! new local variables after BEM
!
REAL, DIMENSION(SIZE(PTA)) :: ZIMB_RF      ! residual energy imbalance
                                             ! of the roof for
                                             ! verification
REAL, DIMENSION(SIZE(PTA)) :: ZIMB_RD      ! road residual energy imbalance 
                                             ! for verification [W m-2]
REAL, DIMENSION(SIZE(PTA)) :: ZIMB_WL      ! wall residual energy imbalance 
                                             ! for verification [W m-2]
REAL, DIMENSION(SIZE(PTA)) :: ZTS_RD       ! road surface temperature 
!                                            ! at previous time-step
REAL, DIMENSION(SIZE(PTA)) :: ZTS_WL_A     ! wall A surface temperature 
!                                            ! at previous time-step
REAL, DIMENSION(SIZE(PTA)) :: ZTS_WL_B     ! wall B surface temperature 
!                                            ! at previous time-step
REAL, DIMENSION(SIZE(PTA)) :: ZTS_WL       ! averaged wall surface temperature 
!                                            ! at previous time-step
REAL, DIMENSION(SIZE(PTA)) :: ZTS_RF       ! roof surface temperature 
!                                            ! at previous time-step
REAL, DIMENSION(SIZE(PTA),SIZE(T%XT_WALL_A,2)) :: ZT_WL ! averaged wall surface temperature 
!
INTEGER :: IWL, IRF                      ! number of wall, roof layer
REAL, DIMENSION(SIZE(PTA)) :: ZRADHT_IN     ! Indoor radiant heat transfer coefficient
                                                    ! [W K-1 m-2]
REAL, DIMENSION(SIZE(PTA)) :: ZTS_FL       ! floor surface temperature [K]
REAL, DIMENSION(SIZE(PTA)) :: ZRAD_RF_WL  ! rad. flux from roof to averaged wall [W m-2(roof)]
REAL, DIMENSION(SIZE(PTA)) :: ZRAD_RF_WIN   ! rad. flux from roof to window [W m-2(roof)]
REAL, DIMENSION(SIZE(PTA)) :: ZRAD_RF_FL ! rad. flux from roof to floor [W m-2(roof)]
REAL, DIMENSION(SIZE(PTA)) :: ZRAD_RF_MA  ! rad. flux from roof to mass [W m-2(roof)]
REAL, DIMENSION(SIZE(PTA)) :: ZCONV_RF_BLD  ! rad. flux from roof to bld [W m-2(roof)]
REAL, DIMENSION(SIZE(PTA)) :: ZRAD_WL_FL ! rad. flux from averaged wall to floor [W m-2(wall)]
REAL, DIMENSION(SIZE(PTA)) :: ZRAD_WL_MA  ! rad. flux from averaged wall to mass [W m-2(wall)]
REAL, DIMENSION(SIZE(PTA)) :: ZRAD_WIN_FL  ! rad. flux from averaged wall to floor [W m-2(win)]
REAL, DIMENSION(SIZE(PTA)) :: ZRAD_WIN_MA   ! rad. flux from averaged wall to mass [W m-2(win)]
REAL, DIMENSION(SIZE(PTA)) :: ZCONV_WL_BLD  ! rad. flux from roof to bld [W m-2(wall)]
REAL, DIMENSION(SIZE(PTA)) :: ZCONV_WIN_BLD   ! rad. flux from roof to bld [W m-2(win)]
REAL, DIMENSION(SIZE(PTA)) :: ZAC_WIN         ! window aerodynamic conductance

REAL, DIMENSION(SIZE(PTA)) :: ZLOAD_IN_RF   ! indoor load on roof W/m2[roof]
REAL, DIMENSION(SIZE(PTA)) :: ZLOAD_IN_FL   ! indoor load on floor W/m2[floor]
REAL, DIMENSION(SIZE(PTA)) :: ZLOAD_IN_WL   ! indoor load on wall W/m2[wall]
REAL, DIMENSION(SIZE(PTA)) :: ZLOAD_IN_WIN   ! indoor load on win W/m2[win]
REAL, DIMENSION(SIZE(PTA)) :: ZLOAD_IN_MA   ! indoor load on mass W/m2[mass]
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE                                             
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TEB',0,ZHOOK_HANDLE)
!
!*      1.     Initializations
!              ---------------
!
!*      1.1    Water reservoirs
!              ----------------
!
ZWS_RF_MAX =  1. ! (1mm) maximum deepness of roof water reservoir
ZWS_RD_MAX =  1. ! (1mm) maximum deepness of road water reservoir
!
!*      1.2    radiative snow variables at previous time-step
!              ----------------------------------------------
!
ZTSSN_RD(:) = T%TSNOW_ROAD%TS(:)
!
!
!*      1.3    indoor aerodynamique conductance for DEF case
!              ----------------------------------------------
!
ZAC_BLD(:) = XUNDEF
IF (TOP%CBEM=='DEF') ZAC_BLD=1. / 0.123 / (XCPD * PRHOA(:)) !* (normalized by rho Cp for convenience)
!-------------------------------------------------------------------------------
!
!*      1.3    number of roof/wall layer
!              -------------------------
!
IWL = SIZE(T%XT_WALL_A,2)
IRF = SIZE(T%XT_ROOF,2)
!
ZTS_WL_A  (:)=T%XT_WALL_A   (:,1)
ZTS_WL_B  (:)=T%XT_WALL_B   (:,1)
ZTS_WL    (:)=0.5 * (ZTS_WL_A(:)+ZTS_WL_B(:))
ZTS_RD    (:)=T%XT_ROAD     (:,1)
ZTS_RF    (:)=T%XT_ROOF     (:,1)
!
!
!*      1.4    load on indoor walls
!              -------------------------
!
IF (TOP%CBEM=='BEM') THEN
  !
  ZLOAD_IN_RF = B%XF_FLOOR_WIN * DMT%XTR_SW_WIN + DMT%XQIN * B%XN_FLOOR * (1-B%XQIN_FLAT) * B%XQIN_FRAD  &
           / (2 + T%XWALL_O_BLD + B%XGLAZ_O_BLD + B%XMASS_O_BLD ) ! W/m2 [ROOF]
  ZLOAD_IN_FL = B%XF_FLOOR_WIN * DMT%XTR_SW_WIN + DMT%XQIN * B%XN_FLOOR * (1-B%XQIN_FLAT) * B%XQIN_FRAD  &
           / (2 + T%XWALL_O_BLD + B%XGLAZ_O_BLD + B%XMASS_O_BLD )
  ZLOAD_IN_MA = B%XF_MASS_WIN * DMT%XTR_SW_WIN + DMT%XQIN * B%XN_FLOOR * (1-B%XQIN_FLAT) * B%XQIN_FRAD  &
           / (2 + T%XWALL_O_BLD + B%XGLAZ_O_BLD + B%XMASS_O_BLD )
  ZLOAD_IN_WL = B%XF_WALL_WIN * DMT%XTR_SW_WIN + DMT%XQIN * B%XN_FLOOR * (1-B%XQIN_FLAT) * B%XQIN_FRAD  &
           / (2 + T%XWALL_O_BLD + B%XGLAZ_O_BLD + B%XMASS_O_BLD )
  ZLOAD_IN_WIN = B%XF_WIN_WIN * DMT%XTR_SW_WIN + DMT%XQIN * B%XN_FLOOR * (1-B%XQIN_FLAT) * B%XQIN_FRAD  &
           / (2 + T%XWALL_O_BLD + B%XGLAZ_O_BLD + B%XMASS_O_BLD )
ELSE
  ZLOAD_IN_RF = 0.
  ZLOAD_IN_FL = 0.
  ZLOAD_IN_MA = 0.
  ZLOAD_IN_WL = 0.
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
ZWS_RF_MAX(:) = ZWS_RF_MAX(:) * PDF_RF(:)
ZWS_RD_MAX(:) = ZWS_RD_MAX(:) * PDF_RD(:)
!
!-------------------------------------------------------------------------------
!
!*      3.     Surface drag
!              ------------
!
 CALL URBAN_DRAG(TOP, T, B, HIMPLICIT_WIND, PTSTEP, PT_CANYON, PQ_CANYON, &
                 PU_CANYON, PT_LOWCAN, PQ_LOWCAN, PU_LOWCAN, PZ_LOWCAN, &
                 ZTS_RF, ZTS_RD, ZTS_WL, PTS_GARDEN, PDN_RF, PDN_RD,    &
                 PEXNS, PEXNA, PTA, PQA, PPS, PRHOA, PZREF, PUREF,      &
                 PVMOD, ZWS_RF_MAX, ZWS_RD_MAX, PPEW_A_COEF,            &
                 PPEW_B_COEF, PPEW_A_COEF_LOWCAN, PPEW_B_COEF_LOWCAN,   &
                 PQSAT_RF, PQSAT_RD, PDELT_RF, PDELT_RD, PCD, PCDN,     &
                 PAC_RF, PAC_RF_WAT, PAC_WL, PAC_RD, PAC_RD_WAT,        &
                 PAC_TOP, PAC_GARDEN, PRI_TWN, PUW_RD, PUW_RF,          &
                 PDUWDU_RD, PDUWDU_RF, PUSTAR_TWN, ZAC_WIN    )
!
!* area-averaged heat transfer coefficient
!
ZVMOD(:) = WIND_THRESHOLD(PVMOD(:),PUREF(:))
!
PCH_TWN(:) = (T%XBLD(:) * PAC_RF(:) + (1.-T%XBLD(:)) * PAC_TOP (:)) / ZVMOD(:)
!
!* aggregation of momentum fluxes for roofs (=> derivate of flux also recalculated)
!
PUW_RF (:) = (1-T%XGREENROOF(:)) * PUW_RF(:) + T%XGREENROOF(:) * PUW_GR(:)
WHERE (PVMOD(:)/=0.) PDUWDU_RF(:) = 2. * PUW_RF(:) / PVMOD(:)
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
 CALL URBAN_SNOW_EVOL(T, B, PT_LOWCAN, PQ_LOWCAN, PU_LOWCAN, ZTS_RF, ZTS_RD, ZTS_WL_A,    &
                      ZTS_WL_B, PPS, ZTA, ZQA, PRHOA, PLW_RAD, PSR, PZREF, PUREF, PVMOD,  &
                      PTSTEP,  PZ_LOWCAN, PDN_RF, DMT%XABS_SW_SNOW_ROOF,                 &
                      DMT%XABS_LW_SNOW_ROOF, PDN_RD, DMT%XABS_SW_SNOW_ROAD,             &
                      DMT%XABS_LW_SNOW_ROAD, PRNSN_RF, PHSN_RF, PLESN_RF, PGSN_RF,       &
                      PMELT_RF, PRNSN_RD, PHSN_RD, PLESN_RD, PGSN_RD, PMELT_RD,           &
                      PLW_WA_TO_NR, PLW_WB_TO_NR, PLW_S_TO_NR, PLW_WIN_TO_NR, ZDQS_SN_RF, &
                      ZDQS_SN_RD      )
!
!-------------------------------------------------------------------------------
!
!*      6.    LW properties
!              -------------
!
PDF_RD (:) = 1. - PDN_RD (:)
!
!-------------------------------------------------------------------------------
!
!*      7.    Indoor radiative temperature
!              ---------------------------
!
! uses the averaged temperature of both walls for the building energy balance
ZT_WL   (:,:)=0.5 * (T%XT_WALL_A(:,:)+T%XT_WALL_B(:,:))
!
SELECT CASE(TOP%CBEM)
   CASE("DEF")
      ZTS_FL(:) = 19. + XTT
      PT_RAD_IND(:) = ( T%XWALL_O_HOR(:) / T%XBLD(:) * ZT_WL(:,IWL) + &
                  T%XT_ROOF(:,IRF) + ZTS_FL(:) ) / (T%XWALL_O_HOR(:) / T%XBLD(:) + 1. + 1.) 
      ZRADHT_IN(:) = XUNDEF
   CASE("BEM")
      ZTS_FL(:) = B%XT_FLOOR(:,1)
      PT_RAD_IND(:)  = (B%XT_MASS(:,1)*B%XMASS_O_BLD(:) + ZT_WL(:,IWL)*T%XWALL_O_BLD(:)     &
                   + ZTS_FL(:) + T%XT_ROOF(:,IRF) + B%XT_WIN2(:) * B%XGLAZ_O_BLD(:)) &
                   /(B%XMASS_O_BLD(:) + T%XWALL_O_BLD(:) + 1. + 1. + B%XGLAZ_O_BLD(:))
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

 CALL ROOF_LAYER_E_BUDGET(TOP, T, B, PQSAT_RF, ZAC_BLD, PTSTEP, PDN_RF, PRHOA,    &
                          PAC_RF, PAC_RF_WAT, PLW_RAD, PPS, PDELT_RF, ZTA, ZQA,   &
                          PEXNA, PEXNS, DMT%XABS_SW_ROOF, PGSN_RF,  ZFLX_BLD_RF, &
                          ZDQS_RF, DMT%XABS_LW_ROOF, DMT%XH_ROOF, PLEW_RF, ZIMB_RF, &
                          DMT%XG_GREENROOF_ROOF, ZRADHT_IN, ZTS_FL, ZT_WL(:,IWL),&
                          ZRAD_RF_WL, ZRAD_RF_WIN, ZRAD_RF_FL, ZRAD_RF_MA, ZCONV_RF_BLD, &
                          PRR, & !modif to add heating/cooling of rain
                          ZLOAD_IN_RF )
!
!-------------------------------------------------------------------------------
!
!*      8.    Road Ts computations
!              -----------------------------
!
!* Road watering

 CALL TEB_IRRIG(TIR%LPAR_RD_IRRIG, PTSTEP, TOP%TTIME%TDATE%MONTH, PTSUN,   &
               TIR%XRD_START_MONTH, TIR%XRD_END_MONTH, TIR%XRD_START_HOUR, &
               TIR%XRD_END_HOUR, TIR%XRD_24H_IRRIG, DMT%XIRRIG_ROAD      )

!* ts_road, ts_wall, qsat_road, t_canyon and q_canyon are updated
!
 CALL ROAD_LAYER_E_BUDGET(T, B, PTSTEP, PDN_RD, PRHOA, PAC_RD, PAC_RD_WAT, &
                          PLW_RAD, PPS, PQSAT_RD, PDELT_RD, PEXNS,         &
                          DMT%XABS_SW_ROAD, PGSN_RD, PQ_LOWCAN, PT_LOWCAN,&
                          ZTS_WL_A, ZTS_WL_B, ZTSSN_RD,  PTS_GARDEN,       &
                          PLW_WA_TO_R, PLW_WB_TO_R, PLW_S_TO_R,            &
                          PLW_WIN_TO_R, PEMIT_LW_RD, ZDQS_RD, DMT%XABS_LW_ROAD,  &
                          DMT%XH_ROAD, PLEW_RD, ZIMB_RD, PRR+DMT%XIRRIG_ROAD    )
!
!-------------------------------------------------------------------------------
!
!*      8.     Wall Ts computations
!              -----------------------------
!
 CALL FACADE_E_BUDGET(TOP, T, B, DMT, PTSTEP, PDN_RD, PRHOA, PAC_WL, ZAC_BLD,   &
                      PLW_RAD, PPS, PEXNS, PT_CANYON, ZTS_RD, ZTSSN_RD, PTS_GARDEN, &
                      ZTS_FL, PLW_WA_TO_WB, PLW_R_TO_WA, PLW_R_TO_WB,      &
                      PLW_G_TO_WA, PLW_G_TO_WB, PLW_S_TO_WA, PLW_S_TO_WB,  &
                      PLW_NR_TO_WA, PLW_NR_TO_WB, PLW_WIN_TO_WA,           &
                      PLW_WIN_TO_WB, PLW_S_TO_WIN, PLW_WA_TO_WIN,          &
                      PLW_WB_TO_WIN, PLW_R_TO_WIN, PLW_G_TO_WIN,           &
                      PLW_NR_TO_WIN, ZFLX_BLD_WL_A, ZDQS_WL_A,             &
                      ZFLX_BLD_WL_B, ZDQS_WL_B, PEMIT_LW_FAC, ZIMB_WL,     &
                      ZRADHT_IN, ZRAD_RF_WL, ZRAD_RF_WIN, ZRAD_WL_FL,      &
                      ZRAD_WL_MA, ZRAD_WIN_FL, ZRAD_WIN_MA, ZCONV_WL_BLD,  &
                      ZCONV_WIN_BLD, ZAC_WIN, ZLOAD_IN_WL, ZLOAD_IN_WIN   )
!
!-------------------------------------------------------------------------------
!
!*      9.     Evolution of interior building air temperature
!              ----------------------------------------------
!
! uses the averaged temperature of both walls for the building energy balance
ZT_WL   (:,:)=0.5 * (T%XT_WALL_A(:,:)+T%XT_WALL_B(:,:))
!
SELECT CASE(TOP%CBEM)
CASE("DEF")
!
   CALL BLD_E_BUDGET(.TRUE., PTSTEP, T%XBLD, T%XWALL_O_HOR,        &
                     PRHOA, T%XT_ROOF, ZT_WL, B%XTI_BLD, ZTS_FL(:) )

   !variables that needs to be computed apart
   B%XQI_BLD = 0.5 * QSAT(B%XTI_BLD, PPS)
   !variables that need to be set 0 for calculation
   ZFLX_BLD_FL(:) = 0.
   ZFLX_BLD_MA (:) = 0.
   !other variables
   PHU_BLD(:)     = XUNDEF

CASE("BEM")
  CALL BEM(BOP, T, B, DMT, PTSTEP, PTSUN, KDAY, PPS, PRHOA, PT_CANYON, &
           PQ_CANYON, PU_CANYON, PHU_BLD, PT_RAD_IND, ZFLX_BLD_FL,&
           ZFLX_BLD_MA, ZRADHT_IN, ZRAD_RF_MA, ZRAD_RF_FL,        &
           ZRAD_WL_MA, ZRAD_WL_FL, ZRAD_WIN_MA, ZRAD_WIN_FL,      &
           ZCONV_RF_BLD, ZCONV_WL_BLD, ZCONV_WIN_BLD, ZLOAD_IN_FL,&
           ZLOAD_IN_MA                                 )

   DMT%XH_WASTE  = DMT%XH_WASTE  * T%XBLD
   DMT%XLE_WASTE = DMT%XLE_WASTE * T%XBLD
END SELECT
!
!-------------------------------------------------------------------------------
!
!*      10.    Fluxes over built surfaces
!              --------------------------
!
 CALL URBAN_FLUXES   (TOP, T, B, DMT, HIMPLICIT_WIND, PT_CANYON, PPEW_A_COEF, PPEW_B_COEF,      &
                      PEXNS, PRHOA, PVMOD, PH_TRAFFIC, PLE_TRAFFIC,PAC_WL, PCD, PDF_RF,         &
                      PDN_RF, PDF_RD, PDN_RD, PRNSN_RF, PHSN_RF, PLESN_RF, PGSN_RF,             &
                      PRNSN_RD, PHSN_RD, PLESN_RD, PGSN_RD, PMELT_RF, ZDQS_RF, PMELT_RD,        &
                      ZDQS_RD, ZDQS_WL_A, ZDQS_WL_B, ZFLX_BLD_RF, ZFLX_BLD_WL_A,                &
                      ZFLX_BLD_WL_B, ZFLX_BLD_FL, ZFLX_BLD_MA, PE_SHADING, PLEW_RF,             &
                      PRN_GR, PH_GR, PLE_GR, PGFLUX_GR,                                         &
                      PLEW_RD, PLE_WL_A, PLE_WL_B, ZMELT_BLT, PUSTAR_TWN                        )
!
!
! Water transfer from snow reservoir to water reservoir in case of snow melt
!
WHERE (PMELT_RF(:) .GT. 0.)
  T%XWS_ROOF(:) = MIN(ZWS_RF_MAX,T%XWS_ROOF(:) + PMELT_RF(:)*PTSTEP)
ENDWHERE
!
WHERE (PMELT_RD(:) .GT. 0.)
  T%XWS_ROAD(:) = MIN(ZWS_RD_MAX,T%XWS_ROAD(:) + PMELT_RD(:)*PTSTEP)
ENDWHERE
!
!-------------------------------------------------------------------------------
!
!*      11.    Roof ans road reservoirs evolution
!              ----------------------------------
!
 CALL URBAN_HYDRO(ZWS_RF_MAX, ZWS_RD_MAX, T%XWS_ROOF, T%XWS_ROAD, PRR,          &
                  DMT%XIRRIG_ROAD, PTSTEP, T%XBLD, DMT%XLE_ROOF, DMT%XLE_ROAD,  &
                  DMT%XRUNOFF_STRLROOF, DMT%XRUNOFF_ROAD   )
!
IF (TOP%LGREENROOF) THEN
  DMT%XRUNOFF_ROOF(:) =  (1.-T%XGREENROOF(:)) * DMT%XRUNOFF_STRLROOF(:) &
                        + T%XGREENROOF(:) * (PRUNOFF_GR(:) + PDRAIN_GR(:))
ELSE
  DMT%XRUNOFF_ROOF(:) =  DMT%XRUNOFF_STRLROOF(:)
ENDIF                                                      
!
!-------------------------------------------------------------------------------
!
!*      19.    Compute aerodynamical resistance 
!              --------------------------------
!
PRESA_TWN(:) = 1. / ( T%XBLD(:) * PAC_RF(:)  + ( 1. - T%XBLD(:)) * PAC_TOP (:))
!
IF (LHOOK) CALL DR_HOOK('TEB',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE TEB
