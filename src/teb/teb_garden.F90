!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL licence
!SFX_LIC version 2.1. See Licence_CeCILL_V2.1-en.txt and Licence_CeCILL_V2.1-fr.txt  
!SFX_LIC for details.
!     #########
    SUBROUTINE TEB_GARDEN (TOP, T, BOP, B, TPN, TIR, DMT,                                         &
                           HIMPLICIT_WIND, PTSUN, PT_CAN, PQ_CAN, PU_CAN, PT_LOWCAN, PQ_LOWCAN,   &
                           PU_LOWCAN, PZ_LOWCAN, PPEW_A_COEF, PPEW_B_COEF, PPEW_A_COEF_LOWCAN,    &
                           PPEW_B_COEF_LOWCAN, PPS, PPA, PEXNS, PEXNA, PTA, PQA, PRHOA, PCO2,     &
                           PLW_RAD, PDIR_SW, PSCA_SW, PSW_BANDS, KSW, PZENITH, PAZIM, PRR, PSR,   &
                           PZREF, PUREF, PVMOD, PH_TRAFFIC, PLE_TRAFFIC, PTSTEP, PLEW_RF, PLEW_RD,&
                           PLE_WL_A, PLE_WL_B, PRNSN_RF, PHSN_RF, PLESN_RF, PGSN_RF, PMELT_RF,    &
                           PRNSN_RD, PHSN_RD, PLESN_RD, PGSN_RD, PMELT_RD, PRN_GRND, PH_GRND,     &
                           PLE_GRND, PGFLX_GRND, PRN_TWN, PH_TWN, PLE_TWN, PGFLX_TWN, PEVAP_TWN,  &
                           PSFCO2, PUW_GRND, PUW_RF, PDUWDU_GRND, PDUWDU_RF,                      &
                           PUSTAR_TWN, PCD, PCDN, PCH_TWN, PRI_TWN, PTS_TWN, PEMIS_TWN,           &
                           PDIR_ALB_TWN, PSCA_ALB_TWN, PRESA_TWN, PAC_RD, PAC_GD, PAC_GR,         &
                           PAC_RD_WAT, PAC_GD_WAT, PAC_GR_WAT, KDAY, PEMIT_LW_FAC, PEMIT_LW_GRND, &
                           PT_RAD_IND, PREF_SW_GRND, PREF_SW_FAC, PHU_BLD, PTIME, PPROD_BLD  )
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
!!                04/2012 add B%XTRAN_WIN
!!    modified    08/2012 TOP%CCH_BEM/ ROUGH_WALL and ROUGH_ROOF for buildind conv coef.
!!    modified    10/2012 add B%XF_WIN_WIN as arg
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_TEB_PANEL_n, ONLY : TEB_PANEL_t
USE MODD_TEB_IRRIG_n, ONLY : TEB_IRRIG_t
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_t
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
USE MODI_TEB_VEG_PROPERTIES
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
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(BEM_OPTIONS_t), INTENT(INOUT) :: BOP
TYPE(BEM_t), INTENT(INOUT) :: B
TYPE(TEB_PANEL_t), INTENT(INOUT) :: TPN
TYPE(TEB_IRRIG_t), INTENT(INOUT) :: TIR
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: DMT
!
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND      ! wind implicitation option
!                                                         ! 'OLD' = direct
!                                                         ! 'NEW' = Taylor serie, order 1
REAL, DIMENSION(:),   INTENT(IN)    :: PTSUN              ! solar time   (s from midnight)
!                                                         
REAL, DIMENSION(:)  , INTENT(INOUT) :: PT_CAN             ! canyon air temperature
REAL, DIMENSION(:)  , INTENT(INOUT) :: PQ_CAN             ! canyon air specific humidity
REAL, DIMENSION(:)  , INTENT(IN)    :: PU_CAN             ! canyon hor. wind
REAL, DIMENSION(:)  , INTENT(IN)    :: PU_LOWCAN          ! wind near the road
REAL, DIMENSION(:)  , INTENT(IN)    :: PT_LOWCAN          ! temp. near the road
REAL, DIMENSION(:)  , INTENT(IN)    :: PQ_LOWCAN          ! hum. near the road
REAL, DIMENSION(:)  , INTENT(IN)    :: PZ_LOWCAN          ! height of atm. var. near the road
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
REAL, DIMENSION(:)  , INTENT(IN)    :: PZREF              ! reference height of the first atm level (temperature)
REAL, DIMENSION(:)  , INTENT(IN)    :: PUREF              ! reference height of the first atm level (wind)
REAL                , INTENT(IN)    :: PTSTEP             ! time step
!
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLEW_RF          ! latent heat flux over roof (snow)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLEW_RD          ! latent heat flux over road (snow)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_WL_A         ! latent heat flux over wall
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_WL_B         ! latent heat flux over wall

!
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRNSN_RF       ! net radiation over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PHSN_RF        ! sensible heat flux over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLESN_RF       ! latent heat flux over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGSN_RF        ! flux under the snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PMELT_RF         ! snow melt
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRNSN_RD       ! net radiation over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PHSN_RD        ! sensible heat flux over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLESN_RD       ! latent heat flux over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGSN_RD        ! flux under the snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PMELT_RD       ! snow melt
!
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_GRND          ! net radiation over ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_GRND           ! sensible heat flux over ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_GRND          ! latent heat flux over ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLX_GRND        ! flux through the ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_TWN           ! net radiation over town
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_TWN            ! sensible heat flux over town
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_TWN           ! latent heat flux over town
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLX_TWN         ! flux through the ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PEVAP_TWN         ! evaporation flux (kg/m2/s)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PSFCO2            ! flux of CO2       (m/s*kg_CO2/kg_air)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PUW_GRND          ! momentum flux for ground built surf
REAL, DIMENSION(:)  , INTENT(OUT)   :: PUW_RF            ! momentum flux for roofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDUWDU_GRND       !
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDUWDU_RF         !
REAL, DIMENSION(:)  , INTENT(OUT)   :: PUSTAR_TWN        ! friciton velocity over town
REAL, DIMENSION(:)  , INTENT(OUT)   :: PCD               ! town averaged drag coefficient
REAL, DIMENSION(:)  , INTENT(OUT)   :: PCDN              ! town averaged neutral drag coefficient
REAL, DIMENSION(:)  , INTENT(OUT)   :: PCH_TWN           ! town averaged heat transfer coefficient
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRI_TWN           ! town averaged Richardson number
REAL, DIMENSION(:)  , INTENT(OUT)   :: PTS_TWN           ! town surface temperature
REAL, DIMENSION(:)  , INTENT(OUT)   :: PEMIS_TWN         ! town equivalent emissivity
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDIR_ALB_TWN      ! town equivalent direct albedo
REAL, DIMENSION(:)  , INTENT(OUT)   :: PSCA_ALB_TWN      ! town equivalent diffuse albedo
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRESA_TWN         ! town aerodynamical resistance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_RD            ! road conductance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_GD            ! green area conductance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_GR            ! green roof conductance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_RD_WAT        ! road conductance for latent heat
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_GD_WAT        ! green area conductance for latent heat
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_GR_WAT        ! green roof conductance for latent heat
!
! new arguments created after BEM
!
INTEGER             , INTENT(IN)    :: KDAY               ! Simulation day
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
!
!* arguments for solar panels
!
REAL, DIMENSION(:),   INTENT(OUT)     :: PPROD_BLD        ! Averaged     Energy production of solar panel on roofs (W/m2 bld  )
!
!*      0.2    Declarations of local variables
!
REAL, DIMENSION(SIZE(PTA)) :: ZTA            ! air temperature extrapolated at roof level
REAL, DIMENSION(SIZE(PTA)) :: ZQA            ! air humidity extrapolated at roof level
!
REAL, DIMENSION(SIZE(PTA)) :: ZDN_RF       ! snow fraction on roofs
REAL, DIMENSION(SIZE(PTA)) :: ZDN_RD       ! snow fraction on roads
REAL, DIMENSION(SIZE(PTA)) :: ZDF_RF       ! free-snow fraction on roofs
REAL, DIMENSION(SIZE(PTA)) :: ZDF_RD       ! free-snow fraction on roads
REAL, DIMENSION(SIZE(PTA)) :: ZDELT_RD     ! fraction of water on roads
REAL, DIMENSION(SIZE(PTA)) :: ZDELT_RF     ! fraction of water on roofs
REAL, DIMENSION(SIZE(PTA)) :: ZAC_RF       ! roof conductance
REAL, DIMENSION(SIZE(PTA)) :: ZAC_RF_WAT   ! roof water conductance
REAL, DIMENSION(SIZE(PTA)) :: ZAC_WL       ! wall conductance
REAL, DIMENSION(SIZE(PTA)) :: ZAC_TOP      ! top conductance
REAL, DIMENSION(SIZE(PTA)) :: ZQSAT_RD     ! hum of saturation for roads
REAL, DIMENSION(SIZE(PTA)) :: ZQSAT_GD     ! hum of saturation for green areas
REAL, DIMENSION(SIZE(PTA)) :: ZQSAT_RF     ! hum of saturation for roofs
REAL, DIMENSION(SIZE(PTA)) :: ZQSAT_GR     ! hum of saturation for green roofs
!
! coefficients for LW computations over snow (from previous time-step)
!
REAL, DIMENSION(SIZE(PTA)) :: ZTSSN_RF   ! roof snow temp at previous time-step
REAL, DIMENSION(SIZE(PTA)) :: ZTSSN_RD   ! road snow temp at previous time-step
REAL, DIMENSION(SIZE(PTA)) :: ZESN_RF    ! snow emissivity at previous time-step
REAL, DIMENSION(SIZE(PTA)) :: ZESN_RD    ! snow emissivity at previous time-step
!
! incoming shortwave radiation
!
REAL, DIMENSION(SIZE(PTA)) :: ZDIR_SW             ! direct  solar rad
REAL, DIMENSION(SIZE(PTA)) :: ZSCA_SW             ! diffuse solar rad
INTEGER                    :: JSWB
!
! albedo & emissivity
!
REAL, DIMENSION(SIZE(PTA)) :: ZALB_GD    ! albedo     for green areas
REAL, DIMENSION(SIZE(PTA)) :: ZEMIS_GD   ! emissivity for green areas
REAL, DIMENSION(SIZE(PTA)) :: ZALB_GR    ! albedo     for green roofs
REAL, DIMENSION(SIZE(PTA)) :: ZEMIS_GR   ! emissivity for green roofs
!
! radiation received by surfaces
!
REAL, DIMENSION(SIZE(PTA)) :: ZREC_SW_RD      ! solar rad received by roads
REAL, DIMENSION(SIZE(PTA)) :: ZREC_SW_WL_A    ! solar rad received by walls
REAL, DIMENSION(SIZE(PTA)) :: ZREC_SW_WL_B    ! solar rad received by walls
REAL, DIMENSION(SIZE(PTA)) :: ZREC_SW_GD      ! solar rad received by gardens
REAL, DIMENSION(SIZE(PTA)) :: ZREC_SW_SN_RD   ! solar rad received by snow on roads
!
REAL, DIMENSION(SIZE(PTA)) :: ZREC_LW_GD      ! IR rad received by gardens
!
REAL, DIMENSION(SIZE(PTA)) :: ZSW_RAD_GD      ! solar radiation reaching urban green areas
REAL, DIMENSION(SIZE(PTA)) :: ZREC_SW_RF      ! solar rad received by roofs in presence of solar panels
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
REAL, DIMENSION(SIZE(PTA)) :: ZREC_LW_RF          ! Incoming LW on roofs in presence of solar panels
!
! local variable at previous time-step
!
REAL, DIMENSION(SIZE(PTA)) :: ZPET_A_COEF          
REAL, DIMENSION(SIZE(PTA)) :: ZPET_B_COEF          
REAL, DIMENSION(SIZE(PTA)) :: ZPEQ_A_COEF          
REAL, DIMENSION(SIZE(PTA)) :: ZPEQ_B_COEF          
!
REAL, DIMENSION(SIZE(PTA)) :: ZUW_RD       ! momentum flux for roads
REAL, DIMENSION(SIZE(PTA)) :: ZUW_GD       ! momentum flux for green areas
REAL, DIMENSION(SIZE(PTA)) :: ZUW_GR       ! momentum flux for green roofs
REAL, DIMENSION(SIZE(PTA)) :: ZDUWDU_RD    !
!
REAL, DIMENSION(SIZE(PTA)) :: ZAC_AGG_GD   ! aggreg. aeodynamic resistance for green areas
REAL, DIMENSION(SIZE(PTA)) :: ZHU_AGG_GD   ! aggreg. relative humidity for green areas
REAL, DIMENSION(SIZE(PTA)) :: ZAC_AGG_GR   ! aggreg. aeodynamic resistance for green roofs
REAL, DIMENSION(SIZE(PTA)) :: ZHU_AGG_GR   ! aggreg. relative humidity for green roofs
!
!  surfaces relative fractions
!
REAL, DIMENSION(SIZE(PTA)) :: ZRF_FRAC        ! roof, wall and
REAL, DIMENSION(SIZE(PTA)) :: ZWL_FRAC        ! road fractions
REAL, DIMENSION(SIZE(PTA)) :: ZRD_FRAC        ! of exchange surf.
REAL, DIMENSION(SIZE(PTA)) :: ZGD_FRAC        !                  
REAL, DIMENSION(SIZE(PTA)) :: ZTOTS_O_HORS    ! total canyon+roof surface
!                                             ! over horizontal surface
REAL, DIMENSION(SIZE(PTA)) :: ZWL_O_RD        ! wall surface over road surface
REAL, DIMENSION(SIZE(PTA)) :: ZWL_O_GRND      ! wall surface over (road+green area) surface
!
! surface temperatures
!
REAL, DIMENSION(SIZE(PTA)) :: ZMTC_O_GR_R1 ! mean thermal conductivity over distance 
!                                          ! between two layers (bottom GR & roof)
!
! fluxes from green surfaces
!
REAL, DIMENSION(SIZE(PTA)) :: ZSFCO2_GD      ! CO2 fluxes (m/s*kg_CO2/kg_air)
REAL, DIMENSION(SIZE(PTA)) :: ZEMIT_LW_GD    ! LW flux emitted by the garden (W/m2 garden)
REAL, DIMENSION(SIZE(PTA)) :: ZSFCO2_GR      ! CO2 fluxes over greenroofs (m/s*kg_CO2/kg_air)
!
! fluxes from built surfaces
REAL, DIMENSION(SIZE(PTA)) :: ZEMIT_LW_RD    ! LW flux emitted by the road (W/m2 road)
!
! fluxes from/to solar panel
REAL, DIMENSION(SIZE(PTA)) :: ZEMIT_LWDN_PANEL  ! LW flux emitted DOWNWARDS by the solar panel (W/m2 panel)
REAL, DIMENSION(SIZE(PTA)) :: ZEMIT_LWUP_PANEL  ! LW flux emitted UPWARDS   by the solar panel (W/m2 panel)
REAL, DIMENSION(SIZE(PTA)) :: ZEMIT_LW_RF       ! LW flux emitted UPWARDS   by the roof        (W/m2 roof )
!
REAL, DIMENSION(SIZE(PTA)) :: ZRN_GD, ZH_GD, ZLE_GD, ZGFLUX_GD, ZEVAP_GD, ZTSRAD_GD, ZRUNOFF_GD, ZDRAIN_GD, ZIRRIG_GD
REAL, DIMENSIOn(SIZE(PTA)) :: ZRN_GR, ZH_GR, ZLE_GR, ZGFLUX_GR
REAL, DIMENSION(SIZE(PTA)) :: ZEVAP_GR, ZTSRAD_GR, ZRUNOFF_GR, ZDRAIN_GR, ZIRRIG_GR
!
!new local variables for shading
REAL, DIMENSION(SIZE(PTA)) :: ZE_SHADING          ! energy not ref., nor absorbed, nor
                                                  ! trans. by glazing [Wm-2(win)]
LOGICAL, DIMENSION(SIZE(PTA)) :: GSHADE           ! describes if one encounters the
!                                                 ! conditions to close windows

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
!
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
DO JJ=1,SIZE(T%XROAD)
  IF (T%XROAD(JJ) .GT. 0.) THEN
    ZTOTS_O_HORS(JJ) = 1. + T%XWALL_O_HOR(JJ)
    ZRF_FRAC  (JJ) = T%XBLD       (JJ) / ZTOTS_O_HORS(JJ)
    ZWL_FRAC  (JJ) = T%XWALL_O_HOR(JJ) / ZTOTS_O_HORS(JJ)
    ZRD_FRAC  (JJ) = T%XROAD      (JJ) / ZTOTS_O_HORS(JJ)
    ZGD_FRAC  (JJ) = T%XGARDEN    (JJ) / ZTOTS_O_HORS(JJ)
    ZWL_O_RD  (JJ) = ZWL_FRAC     (JJ) / ZRD_FRAC(JJ)
    ZWL_O_GRND(JJ) = ZWL_FRAC     (JJ) / (ZRD_FRAC(JJ)+ZGD_FRAC(JJ))
  ELSEIF (T%XBLD(JJ) .EQ. 0. .AND. T%XROAD(JJ) .EQ. 0.) THEN
    ZTOTS_O_HORS(JJ) = 1. + T%XWALL_O_HOR(JJ)
    ZRF_FRAC  (JJ) = 0.
    ZWL_FRAC  (JJ) = T%XWALL_O_HOR(JJ) / ZTOTS_O_HORS(JJ)
    ZRD_FRAC  (JJ) = 0.
    ZGD_FRAC  (JJ) = T%XGARDEN    (JJ) / ZTOTS_O_HORS(JJ)
    ZWL_O_RD  (JJ) = 0.
    ZWL_O_GRND(JJ) = ZWL_FRAC     (JJ) / ZGD_FRAC(JJ)
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
 CALL SNOW_FRAC_ROAD(T%TSNOW_ROAD%WSNOW(:,1),PSR(:)>0.,ZDN_RD,ZDF_RD)
 CALL SNOW_FRAC_ROOF(T%TSNOW_ROOF%WSNOW(:,1),PSR(:)>0.,ZDN_RF,ZDF_RF)
!
!* new snow albedo
!
WHERE (T%TSNOW_ROAD%WSNOW(:,1)==0. .AND. PSR(:)>0.) T%TSNOW_ROAD%ALB(:) = XANSMAX
WHERE (T%TSNOW_ROOF%WSNOW(:,1)==0. .AND. PSR(:)>0.) T%TSNOW_ROOF%ALB(:) = XANSMAX
!
!*      2.2    If snow was not present at previous time-step but is falling
!              ------------------------------------------------------------
!
WHERE (T%TSNOW_ROAD%WSNOW(:,1)==0. .AND. PSR(:)>0.)
  T%TSNOW_ROAD%ALB (:) = XANSMAX
  T%TSNOW_ROAD%EMIS(:) = XEMISSN
  T%TSNOW_ROAD%TS  (:) = MIN(T%XT_ROAD(:,1), XTT)
END WHERE
WHERE (T%TSNOW_ROOF%WSNOW(:,1)==0. .AND. PSR(:)>0.)
  T%TSNOW_ROOF%ALB (:) = XANSMAX
  T%TSNOW_ROOF%EMIS(:) = XEMISSN
  T%TSNOW_ROOF%TS  (:) = MIN(T%XT_ROOF(:,1), XTT)
END WHERE
!
!*      2.3    Radiative snow variables at previous time-step
!              ----------------------------------------------
!
ZESN_RF  (:) = T%TSNOW_ROOF%EMIS(:)
ZESN_RD  (:) = T%TSNOW_ROAD%EMIS(:)
ZTSSN_RF (:) = T%TSNOW_ROOF%TS  (:)
ZTSSN_RD (:) = T%TSNOW_ROAD%TS  (:)
!
!-------------------------------------------------------------------------------
!
!*      3.     Extrapolation of atmospheric T and q at roof level (for fluxes computation)
!              --------------------------------------------------
!
ZTA(:) = PTA(:) * PEXNS(:) / PEXNA(:)
ZQA(:) = PQA(:) * QSAT(ZTA(:),PPS(:)) / QSAT(PTA(:),PPA(:))
!
!-------------------------------------------------------------------------------
!

!*      4.     Grid-averaged albedo and emissivity of green areas
!              --------------------------------------------------
!
ZALB_GD   = XUNDEF
ZEMIS_GD  = XUNDEF
ZTSRAD_GD = XUNDEF
!
IF (TOP%LGARDEN) THEN
 CALL TEB_VEG_PROPERTIES(PDIR_SW, PSCA_SW, PSW_BANDS, KSW,   &
                        ZTSRAD_GD, ZEMIS_GD, ZALB_GD,       &
                        PTA=PT_LOWCAN                       )
ENDIF
!
! for greenroofs :
!
ZALB_GR   = XUNDEF
ZEMIS_GR  = XUNDEF
ZTSRAD_GR = XUNDEF
!
IF (TOP%LGREENROOF) THEN
 CALL TEB_VEG_PROPERTIES(  PDIR_SW, PSCA_SW, PSW_BANDS, KSW,    &
                           ZTSRAD_GR, ZEMIS_GR, ZALB_GR,        &
                           PTA=PTA                              )
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      5.     Occupation of buildings
!              -----------------------
!
!* when building in unoccupied, target temperature is modified
!
IF (TOP%CBEM=="BEM") THEN
  CALL BLD_OCC_CALENDAR(TOP%TTIME, PTSUN, T, B, 1., DMT%XTCOOL_TARGET, DMT%XTHEAT_TARGET, DMT%XQIN )
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
IF (TOP%CBEM=="BEM") THEN
 CALL WINDOW_SHADING_AVAILABILITY(B%LSHADE, B%XTI_BLD, DMT%XTCOOL_TARGET, GSHADE)
!
!* Reinitialize shading of windows when changing day
!
  WHERE (PTSUN .LT. PTSTEP + 1E-3) B%LSHAD_DAY(:) = .FALSE.
END IF
!
!
!*       5.2   computes solar radiation exchanges
!              ----------------------------------
!
 CALL URBAN_SOLAR_ABS(TOP, T, B, DMT, ZDIR_SW, ZSCA_SW, PZENITH, PAZIM,  &
                      TPN%XFRAC_PANEL, TPN%XALB_PANEL, ZALB_GD,     &
                      T%XSVF_GARDEN, ZALB_GR, ZDN_RF, ZDF_RF,       &
                      ZDN_RD, ZDF_RD, ZREC_SW_RD, ZREC_SW_SN_RD,    &
                      ZREC_SW_WL_A, ZREC_SW_WL_B, ZREC_SW_GD,       &
                      ZREC_SW_RF, PDIR_ALB_TWN, PSCA_ALB_TWN,       &
                      ZSW_RAD_GD, ZREC_SW_WIN,PREF_SW_GRND,         &
                      PREF_SW_FAC, ZE_SHADING, B%LSHAD_DAY, GSHADE )
!
!-------------------------------------------------------------------------------
!
!*      6.     LW properties
!              -------------
!
  CALL URBAN_LW_COEF(B, T, PLW_RAD, ZEMIS_GD, T%TSNOW_ROAD%TS, ZTSRAD_GD,  &
                     ZLW_WA_TO_WB, ZLW_WA_TO_R, ZLW_WB_TO_R,            &
                     ZLW_WA_TO_NR, ZLW_WB_TO_NR, ZLW_WA_TO_G,           &
                     ZLW_WB_TO_G, ZLW_WA_TO_WIN, ZLW_WB_TO_WIN,         &
                     ZLW_R_TO_WA, ZLW_R_TO_WB, ZLW_R_TO_WIN,            &
                     ZLW_G_TO_WA, ZLW_G_TO_WB, ZLW_G_TO_WIN,            &
                     ZLW_S_TO_WA, ZLW_S_TO_WB, ZLW_S_TO_R,              &
                     ZLW_S_TO_NR, ZLW_S_TO_G,ZLW_S_TO_WIN,              &
                     ZLW_WIN_TO_WA, ZLW_WIN_TO_WB, ZLW_WIN_TO_R,        &
                     ZLW_WIN_TO_NR, ZLW_WIN_TO_G, ZLW_NR_TO_WA,         &
                     ZLW_NR_TO_WB, ZLW_NR_TO_WIN          )

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
IF (TOP%LGARDEN) THEN
  ZT_SKY  (:) = (PLW_RAD(:)/XSTEFAN)**0.25
  ZREC_LW_GD(:) = (ZLW_S_TO_G  (:) * (ZT_SKY(:)       - ZTSRAD_GD(:))   &
                 + ZLW_WA_TO_G (:) * (T%XT_WALL_A(:,1)- ZTSRAD_GD(:))   &
                 + ZLW_WB_TO_G (:) * (T%XT_WALL_B(:,1)- ZTSRAD_GD(:))   &
                 + ZLW_WIN_TO_G(:) * (B%XT_WIN1(:)    - ZTSRAD_GD(:)) ) &
                   / ZEMIS_GD(:) + XSTEFAN * ZTSRAD_GD(:)**4
ELSE
  ZREC_LW_GD      (:) = XUNDEF
END IF
!
!*     7.2     Effect of solar panels on incoming LW on roofs
!              ----------------------------------------------
!
IF (TOP%LSOLAR_PANEL) THEN
  !
  ! solar panels downwards surface is supposed to be at air temperature
  ! and to be have an emissivity of 1.
  !
  ZEMIT_LWDN_PANEL = XSTEFAN * PTA**4
  !
  ! note that, for the time being, one considers that the solar panel 
  ! intercept radiation both above roof and greenroofs (if any)
  ZREC_LW_RF = (1.-TPN%XFRAC_PANEL(:)) * PLW_RAD + TPN%XFRAC_PANEL(:) * ZEMIT_LWDN_PANEL
ELSE
  ZEMIT_LWDN_PANEL = XUNDEF
  ZREC_LW_RF       = PLW_RAD
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
IF (TOP%LGARDEN) THEN
!
  CALL GARDEN(HIMPLICIT_WIND, TOP%TTIME, PTSUN, PPEW_A_COEF_LOWCAN, PPEW_B_COEF_LOWCAN, &
              ZPET_A_COEF, ZPEQ_A_COEF, ZPET_B_COEF, ZPEQ_B_COEF, PTSTEP, PZ_LOWCAN,    &
              PT_LOWCAN, PQ_LOWCAN, PEXNS, PRHOA, PCO2, PPS, PRR, PSR, PZENITH,         &
              ZREC_SW_GD, ZREC_LW_GD, PU_LOWCAN,  ZRN_GD, ZH_GD, ZLE_GD, ZGFLUX_GD,     &
              ZSFCO2_GD, ZEVAP_GD, ZUW_GD, ZRUNOFF_GD, PAC_GD, ZQSAT_GD, ZTSRAD_GD,     &
              ZAC_AGG_GD, ZHU_AGG_GD, ZDRAIN_GD, ZIRRIG_GD )

  PAC_GD_WAT(:) = PAC_GD(:)
  DMT%XABS_SW_GARDEN(:) = (1.-ZALB_GD(:)) * ZREC_SW_GD
  DMT%XABS_LW_GARDEN(:) = ZEMIS_GD(:) * ZREC_LW_GD(:) - XSTEFAN * ZEMIS_GD(:) * ZTSRAD_GD(:)**4 
  ZEMIT_LW_GD(:) = XSTEFAN * ZTSRAD_GD(:)**4 + (1 - ZEMIS_GD(:)) / ZEMIS_GD(:) * DMT%XABS_LW_GARDEN(:)

ELSE
  !
  ZRN_GD    (:) = 0.
  ZH_GD     (:) = 0.
  ZLE_GD    (:) = 0.
  ZGFLUX_GD (:) = 0.
  ZEVAP_GD  (:) = 0.
  ZRUNOFF_GD(:) = 0. 
  ! 
  ZTSRAD_GD (:) = XUNDEF
  !
  ZUW_GD     (:) = 0.
  PAC_GD     (:) = 0.
  ZSFCO2_GD  (:) = 0.
  ZQSAT_GD   (:) = XUNDEF
  ZAC_AGG_GD (:) = XUNDEF
  ZHU_AGG_GD (:) = XUNDEF
  PAC_GD_WAT (:) = XUNDEF 
  ZEMIT_LW_GD(:) = 0.
  !
  DMT%XABS_SW_GARDEN (:) = XUNDEF
  DMT%XABS_LW_GARDEN (:) = XUNDEF
  !
ENDIF
!
!*      8.3    Call ISBA for greenroofs
!              -------------------------
!
IF (TOP%LGREENROOF) THEN
  !
  CALL GREENROOF(HIMPLICIT_WIND, TOP%TTIME, PTSUN, PPEW_A_COEF, PPEW_B_COEF,         &
                 ZPET_A_COEF, ZPEQ_A_COEF, ZPET_B_COEF, ZPEQ_B_COEF, PTSTEP, PZREF,  &
                 PUREF, PTA, PQA, PEXNS, PEXNA,PRHOA, PCO2, PPS, PRR, PSR, PZENITH,  &
                 ZREC_SW_RF, ZREC_LW_RF, PVMOD, ZRN_GR, ZH_GR, ZLE_GR,               &
                 ZGFLUX_GR, ZSFCO2_GR, ZEVAP_GR, ZUW_GR,                             &
                 PAC_GR, ZQSAT_GR, ZTSRAD_GR, ZAC_AGG_GR, ZHU_AGG_GR,                &
                 DMT%XG_GREENROOF_ROOF, ZRUNOFF_GR, ZDRAIN_GR, ZIRRIG_GR ) 
  !
  PAC_GR_WAT(:) = PAC_GR(:)
  DMT%XABS_SW_GREENROOF(:) = (1.-ZALB_GR(:)) * ZREC_SW_RF
  DMT%XABS_LW_GREENROOF(:) = ZEMIS_GR * ZREC_LW_RF - XSTEFAN * ZEMIS_GR * ZTSRAD_GR**4

ELSE
  !
  ZRN_GR    (:) = 0.
  ZH_GR     (:) = 0.
  ZLE_GR    (:) = 0.
  ZGFLUX_GR (:) = 0.
  ZEVAP_GR  (:) = 0.
  ZRUNOFF_GR(:) = 0.
  ZDRAIN_GR (:) = 0.  
  !
  ZTSRAD_GR (:) = XUNDEF
  !
  ZUW_GR    (:) = 0.
  PAC_GR    (:) = 0.
  ZSFCO2_GR (:) = 0.
  ZQSAT_GR  (:) = XUNDEF
  ZAC_AGG_GR(:) = XUNDEF
  ZHU_AGG_GR(:) = XUNDEF 
  ZMTC_O_GR_R1(:) = XUNDEF 
  !
  DMT%XIRRIG_GREENROOF (:) = 0.
  DMT%XABS_SW_GREENROOF(:) = XUNDEF
  DMT%XABS_LW_GREENROOF(:) = XUNDEF
  DMT%XG_GREENROOF_ROOF(:) = XUNDEF
  !
ENDIF

END SUBROUTINE TEB_GARDEN2
!-------------------------------------------------------------------------------
SUBROUTINE TEB_GARDEN3
!
!*     9.      Treatment of built covers
!              -------------------------
!
  CALL TEB  (TOP, T, BOP, B, TIR, DMT, HIMPLICIT_WIND,                             &
             PTSUN, PT_CAN, PQ_CAN, PU_CAN, PT_LOWCAN, PQ_LOWCAN, PU_LOWCAN,       &
             PZ_LOWCAN, PPEW_A_COEF, PPEW_B_COEF, PPEW_A_COEF_LOWCAN,              &
             PPEW_B_COEF_LOWCAN, PPS, PPA, PEXNS, PEXNA, PTA, PQA, PRHOA, PLW_RAD, &
             PRR, PSR, PZREF, PUREF, PVMOD, PH_TRAFFIC, PLE_TRAFFIC, PTSTEP,       &
             ZDF_RF, ZDN_RF, ZDF_RD, ZDN_RD, ZQSAT_RF, ZQSAT_RD, ZDELT_RF,         &
             ZDELT_RD, ZTSRAD_GD, PLEW_RF, ZUW_GR, PLEW_RD, PLE_WL_A,              &
             PLE_WL_B, PRNSN_RF, PHSN_RF, PLESN_RF, PGSN_RF, PMELT_RF,             &
             ZRN_GR, ZH_GR, ZLE_GR, ZGFLUX_GR, ZDRAIN_GR, ZRUNOFF_GR,              &
             PRNSN_RD, PHSN_RD, PLESN_RD, PGSN_RD, PMELT_RD, ZUW_RD, PUW_RF,       &
             ZDUWDU_RD, PDUWDU_RF, PUSTAR_TWN, PCD, PCDN, PCH_TWN, PRI_TWN, PRESA_TWN, &
             ZAC_RF, PAC_RD, ZAC_WL, ZAC_TOP, PAC_GD, ZAC_RF_WAT, PAC_RD_WAT,      &
             ZLW_WA_TO_WB, ZLW_WA_TO_R, ZLW_WB_TO_R, ZLW_WA_TO_NR, ZLW_WB_TO_NR,   &
             ZLW_R_TO_WA, ZLW_R_TO_WB, ZLW_G_TO_WA, ZLW_G_TO_WB, ZLW_S_TO_WA,      &
             ZLW_S_TO_WB, ZLW_S_TO_R, ZLW_S_TO_NR, ZLW_NR_TO_WA, ZLW_NR_TO_WB,     &
             ZLW_NR_TO_WIN, ZLW_WA_TO_WIN, ZLW_WB_TO_WIN, ZLW_G_TO_WIN,            &
             ZLW_R_TO_WIN, ZLW_S_TO_WIN, ZLW_WIN_TO_WA, ZLW_WIN_TO_WB,             &
             ZLW_WIN_TO_R, ZLW_WIN_TO_NR, KDAY, PEMIT_LW_FAC, ZEMIT_LW_RD,         &
             PT_RAD_IND, PHU_BLD, PTIME, ZE_SHADING )
!
!-------------------------------------------------------------------------------
!
!*    10.      Treatment of solar panels
!              -------------------------
!
IF (TOP%LSOLAR_PANEL) THEN
  ! 
  !* LW radiation coming upwards from roofs
  !
  ZEMIT_LW_RF =  ZREC_LW_RF   &
     - (        T%XGREENROOF(:)  *             DMT%XABS_LW_GREENROOF(:)  &
          + (1.-T%XGREENROOF(:)) * ZDF_RF(:) * DMT%XABS_LW_ROOF(:)       &
          + (1.-T%XGREENROOF(:)) * ZDN_RF(:) * DMT%XABS_LW_SNOW_ROOF(:) )

  !
  ! note that, for the time being, one considers that the solar panel 
  ! intercept radiation both above roof and greenroofs (if any)
  CALL SOLAR_PANEL(TPN, DMT, PTSTEP, PTSUN, T%XRESIDENTIAL, ZEMIT_LW_RF, ZEMIT_LWDN_PANEL, &
                   PLW_RAD, PTA, B%XN_FLOOR, PPROD_BLD   )
ELSE
  PPROD_BLD        = XUNDEF  
END IF
!-------------------------------------------------------------------------------
!
!*     11.     Aggregation
!              -----------
!
 CALL AVG_URBAN_FLUXES(TOP, T, B, TPN, DMT,                                               &
                       PTS_TWN, PEMIS_TWN, PT_CAN, PQ_CAN, PT_LOWCAN, PQ_LOWCAN,          &
                       ZTA, ZQA, PRHOA, PPS, PH_TRAFFIC,  PLE_TRAFFIC, ZWL_O_GRND,        &
                       ZESN_RF, ZEMIS_GR, PLW_RAD,  ZAC_RF, ZAC_RF_WAT, ZAC_WL, PAC_RD,  &
                       PAC_RD_WAT, ZAC_TOP, PAC_GD, ZQSAT_GD, ZAC_AGG_GD, ZHU_AGG_GD,     &
                       ZQSAT_RF, ZQSAT_RD, ZDELT_RF, ZDELT_RD, ZRF_FRAC, ZWL_FRAC,        &
                       ZRD_FRAC, ZGD_FRAC, ZTOTS_O_HORS, ZDF_RF, ZDN_RF, ZDF_RD, ZDN_RD,  &
                       PLE_WL_A, PLE_WL_B, PLEW_RF, PLESN_RF, PLEW_RD, PLESN_RD, PHSN_RD, &
                       ZTSRAD_GD, ZRN_GD, ZH_GD, ZLE_GD, ZGFLUX_GD, ZEVAP_GD,             &
                       ZRUNOFF_GD, ZEVAP_GR, ZRUNOFF_GR, ZDRAIN_GR,                       &
                       PRN_GRND, PH_GRND, PLE_GRND, PGFLX_GRND, PRN_TWN, PH_TWN, PLE_TWN, &
                       PGFLX_TWN, PEVAP_TWN, ZEMIT_LW_RD,ZEMIT_LW_GD, PEMIT_LW_GRND, ZEMIS_GD )
!
PSFCO2(:) = T%XGARDEN(:) * ZSFCO2_GD(:) + T%XBLD(:) * T%XGREENROOF(:) * ZSFCO2_GR(:) ! no CO2 flux from built and road yet.
!
!-------------------------------------------------------------------------------
!
!*     12.     Momentum flux for ground built surfaces
!              ---------------------------------------
!
PUW_GRND (:)     = (T%XROAD(:)*ZUW_RD(:) + T%XGARDEN(:)*ZUW_GD(:)) / (T%XROAD(:)+T%XGARDEN(:))
!
PDUWDU_GRND (:)  = 0.
!
END SUBROUTINE TEB_GARDEN3
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE TEB_GARDEN
