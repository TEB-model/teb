!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL licence
!SFX_LIC version 2.1. See Licence_CeCILL_V2.1-en.txt and Licence_CeCILL_V2.1-fr.txt  
!SFX_LIC for details.
!   ##########################################################################
        SUBROUTINE BEM(BOP, T, B, DMT, PTSTEP, PSUNTIME, KDAY, PPS, PRHOA, PT_CAN,  &
                       PQ_CAN, PU_CAN, PHU_BLD, PT_RAD_IND, PFLX_BLD_FL, PFLX_BLD_MA,&
                       PRADHT_IN, PRAD_RF_MA, PRAD_RF_FL, PRAD_WL_MA, PRAD_WL_FL,&
                       PRAD_WIN_MA, PRAD_WIN_FL, PCONV_RF_BLD, PCONV_WL_BLD,     &
                       PCONV_WIN_BLD, PLOAD_IN_FL, PLOAD_IN_MA                )
!   ##########################################################################
!
!!****  *BEM*
!!
!!    PURPOSE
!!    -------
!
!     Computes the temperature and humidity evolution of indoor air, 
!     building energy demand, HVAC energy consumption, 
!     waste heat from HVAC systems, and heat fluxes from indoor to building surfaces.
!
!
!!**  METHOD
!     ------
!
!              NOMENCLATURE: bld  - refers to building plant area; 
!                            floor- refers to building plant area multiplied 
!                                   by the number of floors;
!                            wall - refers to wall area (excluding windows).
!                            win  - refers to window area. 
!                            mass - refers to internal mass area. 
!
!
!        solar radiation transmitted through windows
!        *******************************************
!
!     Qsol_tr_win = Qsol_facade * tr_win * GR 
!
!
!        indoor wall conv/rad heat transfer coefficients
!        ***********************************************
!
!     The calculation of CHTC accounts for favorable or unfavorable convection 
!     depending on the relative position between the hot layer and cold layer
!
! 
!        building energy demand
!        **********************
!
!     Calculation of the cooling and heating, sensible and latent building energy demand.
!     The sensible demand includes the convective heat transfer from indoor surfaces, the 
!     convective fraction of internal heat gains, and sensible infiltration/ventilation heat
!     gains. The latent demand includes the latent fraction of internal heat gains and latent
!     infiltration/ventilation heat gains.  
!
!        surface areas and volummes (referred to m2_bld)
!        ***********************************************
!
!     Awall   =  WALL_O_BLD [m2_wall/m2_bld]
!     Awin    =  GLAZ_O_BLD [m2_win/m2_bld]   
!     Amass   =  2 * N_FLOOR                 [m2_mass/m2_bld]  
!     N_FLOOR  =  BLD_HEIGHT / FLOOR_HEIGHT  [#]
!     Aroof   =  1                           [m2_roof/m2_bld]  
!     Afloor  =  1                           [m2_floor/m2_bld]   
!     Vol_air =  BLD_HEIGHT                  [m3_bld/m2_bld]
!
!
!        evolution of the internal temperature
!        *************************************
!

!                                  dTin  
!     Vol_air * ro_air * cp_air * ---- = h_wall * Awall * (Twall - Tin)
!                                   dt    + h_roof * Aroof * (Troof -Tin)
!                                         + h_floor * Afloor *(Tfloor - Tin)
!                                         + h_mass * Amass * (Tmass - Tin)  
!                                         + h_win * Awin * (Twin - Tin)
!                                         + Qig * (1 - fig_rad) * (1-fig_lat)
!                                         + Vinf * ro_air * cp_air * (Tout - Tin) 
!                                         + Vsys * ro_air * cp_air * (Tsys - Tin) 
!
!
!        evolution of the internal specific humidity
!        *******************************************
!
!                                  dQin  
!      Vol_air * ro_air * lv_air * ---- = Qig * fig_lat
!                                   dt    + Vinf * ro_air * lv_air * (Qout - Qin) 
!                                         + Vsys * ro_air * lv_air * (Qsys - Qin) 
!
!
!        heat fluxes from indoor to surfaces
!        ***********************************
!
!      Qin_wall  = h_wall  * (Tin - Twall)  [W/m2_wall]
!      Qin_roof  = h_roof  * (Tin - Troof)  [W/m2_roof] 
!      Qin_floor = h_floor * (Tin - Tfloor) [W/m2_floor] 
!      Qin_mass  = h_wall  * (Tin - Tmass)  
!                + Qig * fig_rad * (1-fig_lat)/ 2  
!                + Qsol_tr_win              [W/m2_mass]
!
!
!        energy consumption and waste heat from cooling system
!        *****************************************************
!
!      Qhvac  = Qbld / COP
!      Qwaste = Qbld + Qhvac
!
!
!        energy consumption and waste heat from heating system
!        *****************************************************
!
!      Qhvac  = Qbld / Eff
!      Qwaste = Qhvac - Qbld
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
!!      B. Bueno           * Meteo-France *
!!
!!!    MODIFICATIONS
!!    -------------
!!     Original 2010
!!     G. Pigeon nov. 2011: inclusion floor/mass budget inside
!!                          add automatic/manual ventilation
!!                          conserve exchanges with the different surfaces inside 1 time step
!!    G. Pigeon sept. 2012: use of TARP/DOE coef for indoor convection
!!                          use of both T%XT_WALL_A and T%XT_WALL_B for calculations
!!                          the internal mass depth is 1/2 of the floor depth
!!                          add the option of no atmospheric heat releases by HVAC system (B%XF_WATER_COND < 0)
!!    G. Pigeon oct. 2012:  use indoor air density + new solar heat gain distribution
!!    V. Masson May  2013   implicitation of internal building temperature evolution
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_t
!
USE MODD_CSTS,ONLY : XCPD,XSTEFAN,XLVTT,XG, XRV, XRD
 
USE MODE_THERMOS
USE MODE_PSYCHRO
USE MODI_DX_AIR_COOLING_COIL_CV
USE MODI_FLOOR_LAYER_E_BUDGET
USE MODI_MASS_LAYER_E_BUDGET
USE MODE_CONV_DOE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODD_CSTS, ONLY : XTT
!
use MinimalDXCoolingDriver, only: SimMinimalDXCooling
use mode_psychrolib, only: SetUnitSystem, SI, GetHumRatioFromSpecificHum, GetSpecificHumFromHumRatio
!
IMPLICIT NONE
!
!*      0.1    Declarations of arguments
!
TYPE(BEM_OPTIONS_t), INTENT(INOUT) :: BOP
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(BEM_t), INTENT(INOUT) :: B
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: DMT
!
REAL,                INTENT(IN)   :: PTSTEP        ! Time step
REAL, DIMENSION(:),  INTENT(IN)   :: PSUNTIME       ! current solar time since midnight (solar time, s)
INTEGER,             INTENT(IN)   :: KDAY          ! Simulation day
!
REAL, DIMENSION(:),   INTENT(IN)  :: PPS          ! Canyon air pressure [Pa]
REAL, DIMENSION(:),   INTENT(IN)  :: PRHOA        ! Air density at the lowest level [kg m-3]
REAL, DIMENSION(:),   INTENT(IN)  :: PT_CAN    ! Canyon air temperature [K]
REAL, DIMENSION(:),   INTENT(IN)  :: PQ_CAN    ! Canyon air specific humidity [kg kg-1]
REAL, DIMENSION(:),   INTENT(IN)  :: PU_CAN    ! Canyon wind speed (m s-1)
!
REAL, DIMENSION(:),   INTENT(OUT)  :: PHU_BLD       ! Indoor relative humidity 0 < (-) < 1
REAL, DIMENSION(:),   INTENT(IN)  :: PT_RAD_IND    ! Indoor mean radiant temperature [K]
!
REAL, DIMENSION(:),   INTENT(OUT)  :: PFLX_BLD_FL! Heat flux from indoor air to floor 
                                                    ! [W m-2(bld)]
REAL, DIMENSION(:),   INTENT(OUT)  :: PFLX_BLD_MA ! Heat flux from indoor air to mass 
                                                    ! [W m-2(bld)]
REAL, DIMENSION(:),   INTENT(IN) :: PRADHT_IN     ! Indoor radiant heat transfer coefficient
                                                    ! [W K-1 m-2]
REAL, DIMENSION(:)  , INTENT(IN)  :: PRAD_RF_MA  ! Rad. fluxes between roof and mass
REAL, DIMENSION(:)  , INTENT(IN)  :: PRAD_RF_FL ! Rad. fluxes between roof and floor
REAL, DIMENSION(:)  , INTENT(IN)  :: PRAD_WL_MA  ! Rad. fluxes between wall and mass
REAL, DIMENSION(:)  , INTENT(IN)  :: PRAD_WL_FL ! Rad. fluxes between wall and floor
REAL, DIMENSION(:)  , INTENT(IN)  :: PRAD_WIN_MA   ! Rad. fluxes between wind. and mass
REAL, DIMENSION(:)  , INTENT(IN)  :: PRAD_WIN_FL  ! Rad. fluxes between wind. and floor
REAL, DIMENSION(:)  , INTENT(IN)  :: PCONV_RF_BLD  ! Conv. fluxes between roof and indoor air
REAL, DIMENSION(:)  , INTENT(IN)  :: PCONV_WL_BLD  ! Conv. fluxes between wall and indoor air
REAL, DIMENSION(:)  , INTENT(IN)  :: PCONV_WIN_BLD   ! Conv. fluxes between wind. and indoor air
REAL, DIMENSION(:)  , INTENT(IN)  :: PLOAD_IN_FL  ! solar + int heat gain on floor W/m2 [floor]
REAL, DIMENSION(:)  , INTENT(IN)  :: PLOAD_IN_MA   ! solar + int heat gain on floor W/m2 [mass]
!
!*      0.2    Declarations of local variables 
!
INTEGER                        :: IRF        ! Number of roof layers
INTEGER                        :: IWL        ! Number of wall layers
!REAL                           :: ZTCOMF_MAX   ! Maximum comfort temperature for nat.vent [K]
!
!
LOGICAL, DIMENSION(SIZE(B%XTI_BLD)):: GSCHED     ! Day-night schedule flag 
                                               ! *to be transported to inputs*
!
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZF_NIGHT     ! Reduction factor of int.gains at night
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZF_DAY       ! Amplification factor of int.gains at daytime
!
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZAC_IN_MA_COOL, ZAC_IN_FL_COOL, &
                                 ZAC_IN_RF_COOL, ZAC_IN_WL_A_COOL, &
                                 ZAC_IN_WL_B_COOL, ZAC_IN_WIN_COOL   
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZAC_IN_MA_HEAT, ZAC_IN_FL_HEAT, &
                                 ZAC_IN_RF_HEAT, ZAC_IN_WL_A_HEAT, &
                                 ZAC_IN_WL_B_HEAT, ZAC_IN_WIN_HEAT   
!
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZQIN          ! Internal heat gains [W m-2(bld)]
!
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZV_VENT       ! Ventilation flow rate [m3 s-1 m-2(bld)]
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZINF          ! Infiltration flow rate [m3 s-1 m-2(bld)]
!
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZNAT_VENT     ! Nat.vent airflow rate [m3 s-1 m-2(bld)]
REAL,DIMENSION(SIZE(B%XTI_BLD)) :: ZTI_BLD       ! Indoor air temperature at time step t + dt [K]
REAL,DIMENSION(SIZE(B%XTI_BLD)) :: ZTI_ZTI_BLD  ! Indoor air temperature if windows opened
REAL,DIMENSION(SIZE(B%XTI_BLD)) :: ZTI_BLD_CLOSED! Indoor air temperature if windows closed
!
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZQCOOL_TRGT ! Specific humidity cooling setpoing [kg kg-1]
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZQHEAT_TRGT ! Specific humidity heating setpoing [kg kg-1]
!
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZSHR          ! Rated sensible heat rate
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZM_SYS_RAT    ! Auxiliar mass flow rate [kg s-1 m-2(bld)]
!
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZXMIX         ! Outdoor mixing fraction
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZT_MIX        ! Mixing air temperature [K]
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZQ_MIX        ! Mixing air specific humidity [kg kg-1]
!
REAL,DIMENSION(SIZE(B%XTI_BLD)) :: ZQI_BLD       ! Indoor air humidity at time step t + dt [K}
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZWASTE
!
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZDQS_FL
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZIMB_FL
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZDQS_MA
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZIMB_MA
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZLOAD_FL   ! sum of solar and internal loads on floor
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZLOAD_MA    ! sum of solar and internal loads on mass
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZRAD_FL_MA ! Rad. fluxes from floor to mass
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZCONV_FL_BLD ! Conv. fluxes from floor to indoor air
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZCONV_MA_BLD  ! Conv. fluxes from mass to indoor air
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZRHOI  ! indoor air density
!
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZFRAC_VENT
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: PINFCALC       ! Infiltration calculated [vol/h]

REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZPSATCOOL
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZPSATHEAT

REAL :: ZINFI
REAL :: ZMECH

INTEGER :: JJ                                  ! Loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Local variable used for the coupling with MinimalDX
REAL :: OutdoorHumRatio, InletHumRatio, OutletHumRatio, OutdoorTDryBulb, InletTDryBulb, OutletTemperature
REAL :: TotalHeatingEnergyRate ! Unused by BEM/TEB
INTEGER :: FanMode = 1 ! ON
LOGICAL :: PrintWarnings = .FALSE. ! Do not print warnings to stdout.
!
! Psychrolib: use International System of Units
call SetUnitSystem(SI)
!
!!REAL :: ZEXPL = 0.5 !explicit coefficient for internal temperature evol.
!!REAL :: ZIMPL = 0.5 !implicit coef..
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('BEM',0,ZHOOK_HANDLE)
!
  ! 1.   Initializations
  !      ---------------
  !
  IRF  = SIZE(T%XT_ROOF,2)
  IWL  = SIZE(T%XT_WALL_A,2)
  !
  ! Robert: Calculation of infiltration rate as a function of
  !         building characteristics and meteorological conditions
  !
  CALL GET_INFILTRATION(PRHOA,B%XN50,B%XTI_BLD,PT_CAN,PU_CAN,T%XBLD_HEIGHT,PINFCALC)
  !
  ! Conversion AC/H -> [m3 s-1 m-2(bld)]
  !
  ZINF=PINFCALC*T%XBLD_HEIGHT/3600
  !
  ! Robert:
  ! For the case where a permanent mechanical ventilation
  ! system is present, the ventilation due to this
  ! mechanical ventilation is added to the infiltration
  ! rate using the sum of squares. 
  ! CAUTION: this might not always be justified.
  !
  DO JJ=1,SIZE(B%XISMECH,1)
    IF (B%XISMECH(JJ).GT.0.5) THEN
      !
      ZMECH=B%XMECHRATE(JJ)*T%XBLD_HEIGHT(JJ)/3600.
      ZINFI=ZINF(JJ)
      !
      ZINF(JJ)=SQRT(ZMECH**2+ZINFI**2)
      PINFCALC(JJ)=3600.0*ZINF(JJ)/T%XBLD_HEIGHT(JJ)
      !
    ENDIF
  ENDDO
  !
  ! Conversion W/m²/floor -> W/m²(bld)
  ! We cacluate the number of floors directly
  ! as B%XN_FLOOR may have been rounded thus leading to
  ! under/over estimation of interla heat gains.
  ZQIN = B%XQIN * B%XN_FLOOR
  !
  ! Indoor air density
  !
  ZRHOI(:) = PPS(:) / (XRD * B%XTI_BLD(:) * ( 1.+((XRV/XRD)-1.)*B%XQI_BLD(:) ) )
  !
  !###########################################################
  ! 2.   heat balance for building floor and mass
  !###########################################################
  !
  ! 2.1 FLOOR HEAT BALANCE
  !
  CALL FLOOR_LAYER_E_BUDGET(B, PTSTEP, PFLX_BLD_FL, ZDQS_FL, ZIMB_FL, PRADHT_IN, &
                            PRAD_WL_FL, PRAD_RF_FL, PRAD_WIN_FL, PLOAD_IN_FL,    &
                            ZRAD_FL_MA, ZCONV_FL_BLD)
  !
  ! 2.2 MASS HEAT BALANCE
  !
  CALL MASS_LAYER_E_BUDGET(B, PTSTEP, PFLX_BLD_MA, ZDQS_MA, ZIMB_MA, PRADHT_IN, &
                          PRAD_WL_MA, PRAD_RF_MA, PRAD_WIN_MA, PLOAD_IN_MA,    &
                          ZRAD_FL_MA, ZCONV_MA_BLD  )
  !
  ZAC_IN_WL_A_COOL = CHTC_VERT_DOE(T%XT_WALL_A(:,IWL), B%XTCOOL_TARGET)
  ZAC_IN_WL_B_COOL = CHTC_VERT_DOE(T%XT_WALL_B(:,IWL), B%XTCOOL_TARGET)
  ZAC_IN_WIN_COOL  = CHTC_VERT_DOE(B%XT_WIN2         , B%XTCOOL_TARGET)
  ZAC_IN_MA_COOL   = CHTC_VERT_DOE(B%XT_MASS(:,1)    , B%XTCOOL_TARGET)
  ZAC_IN_RF_COOL   = CHTC_DOWN_DOE(T%XT_ROOF(:,IRF)  , B%XTCOOL_TARGET)
  ZAC_IN_FL_COOL   = CHTC_UP_DOE  (B%XT_FLOOR(:,1)   , B%XTCOOL_TARGET)
  !
  ZAC_IN_WL_A_HEAT = CHTC_VERT_DOE(T%XT_WALL_A(:,IWL), B%XTHEAT_TARGET)
  ZAC_IN_WL_B_HEAT = CHTC_VERT_DOE(T%XT_WALL_B(:,IWL), B%XTHEAT_TARGET)
  ZAC_IN_WIN_HEAT  = CHTC_VERT_DOE(B%XT_WIN2         , B%XTHEAT_TARGET)
  ZAC_IN_MA_HEAT   = CHTC_VERT_DOE(B%XT_MASS(:,1)    , B%XTHEAT_TARGET)
  ZAC_IN_RF_HEAT   = CHTC_DOWN_DOE(T%XT_ROOF(:,IRF)  , B%XTHEAT_TARGET)
  ZAC_IN_FL_HEAT   = CHTC_UP_DOE  (B%XT_FLOOR(:,1)   , B%XTHEAT_TARGET)
  !
  ! Calculation of ventilation exchange rate before loop
  !
  CALL GET_NAT_VENT(B%XTI_BLD, PT_CAN, PU_CAN, B%XGR, B%XFLOOR_HW_RATIO, T%XBLD_HEIGHT(:), ZNAT_VENT,B%XFOPEN)
  !
  ! Impose lower limit for coefficients
  !
  ZAC_IN_WL_A_COOL(:) = MAX(1., ZAC_IN_WL_A_COOL(:))
  ZAC_IN_WL_B_COOL(:) = MAX(1., ZAC_IN_WL_B_COOL(:))
  ZAC_IN_WIN_COOL (:) = MAX(1., ZAC_IN_WIN_COOL(:))
  ZAC_IN_MA_COOL  (:) = MAX(1., ZAC_IN_MA_COOL (:))
  ZAC_IN_RF_COOL  (:) = MAX(1., ZAC_IN_RF_COOL (:))
  ZAC_IN_FL_COOL  (:) = MAX(1., ZAC_IN_FL_COOL (:))
  !
  ZAC_IN_WL_A_HEAT(:) = MAX(1., ZAC_IN_WL_A_HEAT(:))
  ZAC_IN_WL_B_HEAT(:) = MAX(1., ZAC_IN_WL_B_HEAT(:))
  ZAC_IN_WIN_HEAT (:) = MAX(1., ZAC_IN_WIN_HEAT(:))
  ZAC_IN_MA_HEAT  (:) = MAX(1., ZAC_IN_MA_HEAT (:))
  ZAC_IN_RF_HEAT  (:) = MAX(1., ZAC_IN_RF_HEAT (:))
  ZAC_IN_FL_HEAT  (:) = MAX(1., ZAC_IN_FL_HEAT (:))
  !
  ! 4. Indoor energy balance calculation
  !
  ! Calculation of ventilation depending on choice of NATVENT
  !
  DO JJ=1,SIZE(PT_CAN)
    IF (B%CNATVENT(JJ)=='NONE') THEN
      !
      ! No ventilation
      !
      B%XDIAGVENT(JJ)=0.0
      !
    ELSEIF (B%CNATVENT(JJ)=='AUTO') THEN
      !
      STOP "Ventilation case not curretnly implemented"
      ! Note that in 8.2 version this is somewhat implemented but
      ! not availabe due to a stop condition placed on entry
      ! `STOP "Implement constraints due to user needs"` therefore
      ! here we simply omit this case atogether.
      !
    ELSEIF (B%CNATVENT(JJ)=='MANU') THEN
      !
      ZFRAC_VENT(JJ) = (1.0 / (1.0 + exp(-1.0*( B%XTI_BLD(JJ) - B%XTDESV(JJ)-2.0)))) * &
                        (1.0 / (1.0 + exp(-1.0*( B%XTI_BLD(JJ) - PT_CAN(JJ)   ))))
      !
      IF ( (ZFRAC_VENT(JJ).LT.-1.0E-10).OR.(ZFRAC_VENT(JJ).GT.(1.0+1.0E-10)) ) STOP "Wrong ventilation fraction"
      !
      B%XDIAGVENT(JJ) = B%XVENT_BEHAV_ANYWAY(JJ) + ZFRAC_VENT(JJ) * B%XVENT_BEHAV_ADAPTI(JJ)
      !
      ! Ventilation fractions below 0.2 are set to 0.0 since ventilation
      ! is very non-linear and interpheres with heating
      !
      IF (B%XDIAGVENT(JJ).LT.0.2) B%XDIAGVENT(JJ)=0.0
      !
      IF ( (B%XDIAGVENT(JJ).LT.-1.0E-10).OR.(B%XDIAGVENT(JJ).GT.(1.0+1.0E-10)) ) STOP "Wrong ventilation fraction"
      !
    ELSE
      STOP "This ventilation case is not implemented"
    ENDIF
    !
    ! Robert: Ventilation for air refreshment independ of indoor/outdoor conditions
    ! This shall consider that people open the windows after getting up,
    ! and/or when they arive home from work and/or before they go to bed, etc. 
    ! At the moment only for manual ventilation.
    ! However, one might consider that people open the window also
    ! for the cases of automatic and mechanical ventilation.
    ! If the house is not occupied the LFVSHO input variable needs to be
    ! set to .FALSE., since in that case the windows remain closed.
    !
    IF ((B%XFVSHO(JJ).GT.0.5).AND.(B%CNATVENT(JJ)=='MANU')) THEN
      !
      STOP "Solution for short-term ventilation needs to be implemented"
      !
    ENDIF
    !
  END DO
  !
  !################################################################################
  ! Calculation of heating and cooling energy demand
  !################################################################################
  !
  ! Re-calculation of ventilation rate
  CALL GET_NAT_VENT(B%XTI_BLD, PT_CAN, PU_CAN, B%XGR, B%XFLOOR_HW_RATIO, T%XBLD_HEIGHT(:), ZNAT_VENT,B%XFOPEN)
  !
  ! In case of manual ventilation, the ventilation flow rate is multiplied
  ! by the fraction of windows that have been opened
  !
  DO JJ=1,SIZE(B%XDIAGVENT)
    IF (B%CNATVENT(JJ)=='MANU') ZNAT_VENT(JJ) = ZNAT_VENT(JJ) * B%XDIAGVENT(JJ)
  ENDDO
  !
  ! If natural or mechanical surventilation ACTIVE
  !
  DO JJ=1,SIZE(B%XDIAGVENT)
    IF (B%XDIAGVENT(JJ).GT.0.0) THEN
      !
      ! The following lines would need to be changed in case 
      ! the HVAC system is not switched off during ventilation.
      !
      DMT%XH_BLD_COOL(JJ) = 0.0
      DMT%XH_BLD_HEAT (JJ) = 0.0
      DMT%XLE_BLD_COOL(JJ) = 0.0
      DMT%XLE_BLD_HEAT(JJ) = 0.0
      !
      DMT%XT_BLD_COOL(JJ) = 0.0
      DMT%XHVAC_COOL (JJ) = 0.0
      DMT%XT_SYS(JJ) = B%XTI_BLD(JJ)
      DMT%XQ_SYS(JJ) = B%XQI_BLD(JJ)
      !
      DMT%XH_WASTE(JJ) = 0.0
      DMT%XLE_WASTE(JJ) = 0.0
      !
      DMT%XHVAC_HEAT = 0.0
      DMT%XM_SYS(JJ) = 0.0
      DMT%XCOP(JJ) = 0.0
      DMT%XCAP_SYS(JJ) = 0.0
      !
    ENDIF
  ENDDO
  !
  ! *If natural surventilation INACTIVE
  !
  WHERE (B%XDIAGVENT(:).LE.0.0)
    !
    ZNAT_VENT(:) = 0.
    !
    ! ------------------------------------------------
    ! Building energy demand for heating and cooling
    ! CAUTION: ventilation = 0.0 assumed
    ! ------------------------------------------------
    !
    DMT%XH_BLD_COOL(:)                                                                            &
      = 0.5* T%XWALL_O_BLD(:) * (ZAC_IN_WL_A_COOL(:) * (T%XT_WALL_A(:,IWL) - B%XTCOOL_TARGET(:))  &
      +                          ZAC_IN_WL_B_COOL(:) * (T%XT_WALL_B(:,IWL) - B%XTCOOL_TARGET(:))) &
      +         B%XGLAZ_O_BLD(:) * ZAC_IN_WIN_COOL(:) * (B%XT_WIN2(:)      - B%XTCOOL_TARGET(:))  &
      +          ZAC_IN_MA_COOL (:)* B%XMASS_O_BLD(:) * (B%XT_MASS(:,1)    - B%XTCOOL_TARGET(:))  &
      +                             ZAC_IN_RF_COOL(:) * (T%XT_ROOF(:,IRF)  - B%XTCOOL_TARGET(:))  &
      +                             ZAC_IN_FL_COOL(:) * (B%XT_FLOOR(:,1)   - B%XTCOOL_TARGET(:))  &
      + ZQIN(:) * (1 - B%XQIN_FRAD(:)) * (1 - B%XQIN_FLAT(:))                                     &
      + ZINF(:)*ZRHOI(:) * XCPD * (PT_CAN(:) - B%XTCOOL_TARGET(:))
    !
    DMT%XH_BLD_HEAT(:)                                                                                &
      = - ( 0.5* T%XWALL_O_BLD(:) * (ZAC_IN_WL_A_HEAT(:) * (T%XT_WALL_A(:,IWL) - B%XTHEAT_TARGET(:))  &
      +                              ZAC_IN_WL_B_HEAT(:) * (T%XT_WALL_B(:,IWL) - B%XTHEAT_TARGET(:))) &
      +             B%XGLAZ_O_BLD(:) * ZAC_IN_WIN_HEAT(:) * (B%XT_WIN2(:)      - B%XTHEAT_TARGET(:))  &
      +               ZAC_IN_MA_HEAT(:)* B%XMASS_O_BLD(:) * (B%XT_MASS(:,1)    - B%XTHEAT_TARGET(:))  &
      +                                 ZAC_IN_RF_HEAT(:) * (T%XT_ROOF(:,IRF)  - B%XTHEAT_TARGET(:))  &
      +                                 ZAC_IN_FL_HEAT(:) * (B%XT_FLOOR(:,1)   - B%XTHEAT_TARGET(:))  &
      + ZQIN(:) * (1 - B%XQIN_FRAD(:))* (1 - B%XQIN_FLAT(:))                                          &
      + ZINF(:) * ZRHOI(:) * XCPD * (PT_CAN(:) - B%XTHEAT_TARGET(:)))
    !
    ! The design specific humidity is calculated based on
    ! the desing temperature and the design relative humidity
    !
    ZPSATCOOL(:) = PSAT(B%XTCOOL_TARGET(:))
    ZPSATHEAT(:) = PSAT(B%XTHEAT_TARGET(:))
    !
    ZQCOOL_TRGT(:) = 0.62198 * B%XHR_TARGET(:) * ZPSATCOOL(:) / &
                      (PPS(:)- B%XHR_TARGET(:) * ZPSATCOOL(:))
    !
    ZQHEAT_TRGT(:) = 0.62198 * B%XHR_TARGET(:) * ZPSATHEAT(:) / &
        (PPS(:)- B%XHR_TARGET(:) * ZPSATHEAT(:))
    !
    ! The latent heat required to reach the target specific humidity is calculated
    !
    DMT%XLE_BLD_COOL(:) = ZQIN(:) * B%XQIN_FLAT(:)+ ZINF(:) * ZRHOI(:) * XLVTT * (PQ_CAN(:) - ZQCOOL_TRGT(:))
    !
    ! No target humidity for heating
    !
    DMT%XLE_BLD_HEAT(:) = 0.0
    !
  ENDWHERE
  !
  ! Autosize calculations: Removed as these only applied for Toulouse case.
  !
  ! Calculation of system efficiency
  !
  DMT%XM_SYS  (:) = B%XM_SYS_RAT  (:)
  DMT%XCOP    (:) = B%XCOP_RAT    (:)
  DMT%XCAP_SYS(:) = B%XCAP_SYS_RAT(:)
  !
  ! Calculation of mixing conditions
  ! Robert: At the moment the mixing conditions are equal
  !         to the indoor conditions since there is no 
  !         cooling/heating when there is ventilation.
  !
  ZV_VENT(:) = 0. ! Not used
  ZXMIX (:) = ZV_VENT(:) * ZRHOI(:) / DMT%XM_SYS(:)
  ZT_MIX(:) = ZXMIX(:) * PT_CAN(:) + (1.-ZXMIX(:)) * B%XTI_BLD(:)
  ZQ_MIX(:) = ZXMIX(:) * PQ_CAN(:) + (1.-ZXMIX(:)) * B%XQI_BLD(:)
  !
  ! ###########################################
  ! Calculations related to cooling
  ! Only for grid points with cooling required
  ! ###########################################
  !
  ! Calculation of cooling system performance and waste heat
  !
  IF (BOP%CCOOL_COIL=='IDEAL') THEN
    DO JJ=1,SIZE(DMT%XH_BLD_COOL)
      IF (DMT%XH_BLD_COOL(JJ) >= 0.0) THEN
        !
        ! Ideal system
        !
        ! desactivation of LE_BLD_COOL impact on HVAC_COOL calculation
        ! following too much impact in VURCA simulation (23/01/2012)
        ! this would be the case for a vaporization system !
        ! PT_BLD_COOL(JJ)=PH_BLD_COOL(JJ)+PLE_BLD_COOL(JJ)
        ! PHVAC_COOL (JJ) = PT_BLD_COOL(JJ) / PCOP_RAT(JJ)
        !
        DMT%XT_BLD_COOL(JJ) = DMT%XH_BLD_COOL(JJ) + DMT%XLE_BLD_COOL(JJ)
        !desactivation of LE_BLD_COOL impact on HVAC_COOL calculation
        !following too much impact in VURCA simulation (23/01/2012)
        !this would be the case for a vaporization system !
        !DMT%XHVAC_COOL (JJ) = DMT%XT_BLD_COOL(JJ) / B%XCOP_RAT(JJ)
        DMT%XHVAC_COOL (JJ) = DMT%XH_BLD_COOL(JJ) / B%XCOP_RAT(JJ)
        IF (DMT%XHVAC_COOL(JJ) < 0.0) DMT%XHVAC_COOL(JJ) = 0.0
        !
        DMT%XT_SYS(JJ) = ZT_MIX(JJ) - DMT%XH_BLD_COOL (JJ)  /DMT%XM_SYS(JJ) / XCPD
        !DMT%XQ_SYS(JJ) = ZQ_MIX(JJ) - DMT%XLE_BLD_COOL(JJ) / DMT%XM_SYS(JJ)/ XLVTT
        !desactivation following too much impact in VURCA simulation
        !(23/01/2012)
        DMT%XQ_SYS(JJ) = ZQ_MIX(JJ)
        !
        DMT%XH_WASTE(JJ)  = DMT%XHVAC_COOL(JJ) * (1.+B%XCOP_RAT(JJ)) * (1. - B%XF_WATER_COND(JJ))
        DMT%XLE_WASTE(JJ) = DMT%XHVAC_COOL(JJ) * (1.+B%XCOP_RAT(JJ)) * B%XF_WATER_COND(JJ)
      ENDIF
    ENDDO
  ELSEIF (BOP%CCOOL_COIL=='DXCOIL') THEN
    DO JJ=1,SIZE(DMT%XH_BLD_COOL)
      IF (DMT%XH_BLD_COOL(JJ).GT.0.0) THEN
        !
        ! FIXME: PQ_CAN is speicifc humidity, model wants mixing ratio. Also output.
        !        Do as done in MinimalDX below.
        CALL DX_AIR_COOLING_COIL_CV(PT_CAN(JJ), PQ_CAN(JJ), PPS(JJ),  ZRHOI(JJ), ZT_MIX(JJ), &
                                    ZQ_MIX(JJ), B%XCOP_RAT(JJ), B%XCAP_SYS_RAT(JJ),          &
                                    B%XT_ADP(JJ), B%XF_WATER_COND(JJ), DMT%XM_SYS(JJ),      &
                                    DMT%XH_BLD_COOL(JJ), DMT%XH_WASTE(JJ), DMT%XLE_WASTE(JJ), &
                                    DMT%XCOP(JJ), DMT%XCAP_SYS(JJ), DMT%XT_SYS(JJ), & 
                                    DMT%XQ_SYS(JJ), DMT%XHVAC_COOL(JJ), DMT%XT_BLD_COOL(JJ) )
      ENDIF
    ENDDO
  ELSEIF (BOP%CCOOL_COIL=='MinimalDX') THEN
    DO JJ=1,SIZE(DMT%XH_BLD_COOL)
      IF (DMT%XH_BLD_COOL(JJ) >= 0.0) THEN
        ! The appratus dew point temperature (B%XT_ADP(JJ); PT_ADP in dx_air_cooling_coil_cv.F90)
        ! should not be a parameter but a variable computed by the coil model. Thefire this is no longer required.
        !
        ! The (water) fraction of evaporation for the condensers (B%XF_WATER_COND(JJ); PF_WATER_COND in 
        ! dx_air_cooling_coil_cv.F90) is set to zero as MinimalDX only supoorts dry evaporators.
        !
        ! XM_SYS is used in BEM to compute the evolution of internal temperature and humidity.
        ! As this is not an output argument in MinimalDX but is modified internally, we manually
        ! preprocess XM_SYS in the same way as done in MinimalDX to ensure consinstent results.
        IF (DMT%XM_SYS(JJ) / ZRHOI(JJ) / B%XCAP_SYS_RAT(JJ) < 0.00004027) THEN
          DMT%XM_SYS(JJ) = 0.00004027 * B%XCAP_SYS_RAT(JJ) * ZRHOI(JJ)
        ELSE IF (DMT%XM_SYS(JJ) / ZRHOI(JJ) / B%XCAP_SYS_RAT(JJ) > 0.00006041) THEN
          DMT%XM_SYS(JJ) = 0.00006041 * B%XCAP_SYS_RAT(JJ) * ZRHOI(JJ)
        END IF
        !
        ! In TEB the humidity is specified in terms of speicifc humidity, in MinimalDX the mixing ratio is used instead.
        OutdoorHumRatio = GetHumRatioFromSpecificHum(PQ_CAN(JJ))
        InletHumRatio = GetHumRatioFromSpecificHum(ZQ_MIX(JJ))
        !
        ! In TEB the temperature is specified in K, in MinimalDX in degree C.
        OutdoorTDryBulb = PT_CAN(JJ) - XTT
        InletTDryBulb = ZT_MIX(JJ)  - XTT
        !
        CALL SimMinimalDXCooling( OutdoorTDryBulb,      &   ! PT_CANYON   : OutdoorTDryBulb
                                  OutdoorHumRatio,      &   ! PQ_CANYON   : OutdoorSpecificHum
                                  PPS(JJ),              &   ! PPS         : OutdoorPressure
                                  InletTDryBulb,        &   ! PT_IN       : InletTDryBulb
                                  InletHumRatio,        &   ! PQ_IN       : InletHumRatio
                                  B%XCOP_RAT(JJ),       &   ! PCOP_RAT    : RatedCOP
                                  B%XCAP_SYS_RAT(JJ),   &   ! PCAP_SYS_RAT: RatedTotCap
                                  DMT%XH_BLD_COOL(JJ),  &   ! PH_BLD_COOL : SensibleCoolingLoad
                                  DMT%XM_SYS(JJ),       &   ! PM_SYS      : RatedAirMassFlowRate
                                  DMT%XCOP(JJ),         &   ! PCOP        : COP
                                  DMT%XCAP_SYS(JJ),     &   ! PCAP_SYS    : TotalCoolingCapacity
                                  OutletTemperature,    &   ! PT_OUT      : OutletTemperature
                                  OutletHumRatio,       &   ! PQ_OUT      : OutletHumRatio
                                  DMT%XHVAC_COOL(JJ),   &   ! PDX_POWER   : ElecCoolingPower
                                  DMT%XLE_WASTE(JJ),    &   ! LE_WASTE    : LatCoolingEnergyRate
                                  DMT%XT_BLD_COOL(JJ),  &   ! PT_BLD_COOL : TotalCoolingEnergyRate
                                  DMT%XH_WASTE(JJ) )        ! PH_WASTE    : TotalSensibleHeatOut
        !
        DMT%XT_SYS(JJ) = OutletTemperature + XTT
        DMT%XQ_SYS(JJ) = GetSpecificHumFromHumRatio(OutletHumRatio)
      ENDIF
    ENDDO
  ELSE
    PRINT *, BOP%CCOOL_COIL
    STOP "This type of cooling system is not implemented"
  ENDIF
  !
  ! Cooling system without atmospheric releases.
  ! e.g. releases in soil/water F_WATER_COND < 0
  ! NOTE: In this case the energy leaving the system would need to
  !       be considered.
  !
  DO JJ=1,SIZE(DMT%XH_BLD_COOL)
    IF (DMT%XH_BLD_COOL(JJ).GT.0.0) THEN
      IF(B%XF_WATER_COND(JJ).LT.0.0) THEN
        DMT%XH_WASTE(JJ) = 0.
        DMT%XLE_WASTE(JJ) = 0.
        ELSE
        ! From EP Engineering Reference (p. 647)
        ! PFAN_POWER(JJ)=PM_SYS(JJ)*PFAN_AP(JJ)*PFAN_EFF(JJ)*ZRHOI(JJ)
        DMT%XH_BLD_HEAT (JJ) = 0.0
        DMT%XLE_BLD_HEAT(JJ) = 0.0
        DMT%XHVAC_HEAT  (JJ) = 0.0
      ENDIF
    ENDIF
  ENDDO
  !
  ! ---------------------------------------------
  ! * HEATING system : Performance and Waste heat
  ! ---------------------------------------------
  !
  IF (MAXVAL(DMT%XH_BLD_HEAT).GT.0.0) THEN
    IF (BOP%CHEAT_COIL=='MinimalDX') THEN
      DO JJ=1,SIZE(DMT%XH_BLD_HEAT)
        IF (DMT%XH_BLD_HEAT(JJ) .GT. 0.0) THEN
          STOP 'Placeholder: feature not implemented in this version'
        END IF
      ENDDO
    ENDIF
    !
    ! If the heating system is not ideal, the delivered
    ! heating is limited to the capacity of the system
    !
    IF (BOP%CHEAT_COIL .NE. 'IDEAL' .AND. BOP%CHEAT_COIL .NE. 'MinimalDX') THEN
      DO JJ=1,SIZE(DMT%XH_BLD_HEAT)
        DMT%XH_BLD_HEAT(JJ)=MIN(DMT%XH_BLD_HEAT(JJ), B%XCAP_SYS_HEAT(JJ))
      ENDDO
    ENDIF
    !
    IF (BOP%CHEAT_COIL .NE. 'MinimalDX') THEN
      WHERE (DMT%XH_BLD_HEAT(:).GT.0.0)
        !
        DMT%XT_SYS(:) = ZT_MIX(:) + DMT%XH_BLD_HEAT(:) / DMT%XM_SYS(:) / XCPD
        DMT%XQ_SYS(:) = ZQ_MIX(:)
        !
        ! The heating energy demand is equal to the
        ! building energy demand divided by the efficiency
        !
        ! For CONVENTIONAL heating systems, the heating energy
        ! demand is HIGHER than the building energy demand
        ! since the efficiency is lower than 1.0.
        ! For HEAT PUMPS, the heating energy demand
        ! is LOWER than the building energy demand since
        ! the efficiency is larger than 1.0.
        !
        DMT%XHVAC_HEAT  (:) = DMT%XH_BLD_HEAT(:) / B%XEFF_HEAT(:)
        ! The waste heat is the difference between the
        ! heating energy demand and the building energy demand.
        ! For EFF_HEAT < 1.0 the waste heat is positive (e.g. waste heat of furnace)
        ! For EFF_HEAT > 1.0 the waste heat is negative since
        !                    the heat is taken from the outside air (air heat pumps)
        !
        ! CAUTION: for a soil heat pump, this equation needs to be
        ! changed since the heat source is no longer the outdoor
        ! air, but the deep soil or a water reservoir.
        !
        DMT%XH_WASTE(:) = DMT%XHVAC_HEAT(:) - DMT%XH_BLD_HEAT(:)
        !
        ! Here it is assumed that there is not latent waste.
        ! This might not be true when there is fuel combustion
        ! for heating.
        !
        !
        ! No cooling
        !
      END WHERE
    END IF
    WHERE (DMT%XH_BLD_HEAT(:).GT.0.0)
      DMT%XH_BLD_COOL(:)=0.0
      DMT%XLE_BLD_COOL(:)=0.0
      DMT%XT_BLD_COOL(:)=0.0
      DMT%XHVAC_COOL(:)=0.0
      !
      ! From EP Engineering Reference (p. 647)
      ! PFAN_POWER(:)=PM_SYS(:)*PFAN_AP(:)*(PFAN_EFF(:)*ZRHOI(:))
      !
    END WHERE
  ENDIF
  !
  ! ------------------------------
  ! * NEITHER COOLING NOR HEATING
  ! ------------------------------
  !
  WHERE ( (DMT%XH_BLD_HEAT(:).LE.0.0).AND.(DMT%XH_BLD_COOL(:).LE.0.0) )
  DMT%XH_BLD_COOL (:) = 0.0
  DMT%XH_BLD_HEAT (:) = 0.0 
  DMT%XLE_BLD_COOL(:) = 0.0
  DMT%XLE_BLD_HEAT(:) = 0.0
  DMT%XT_BLD_COOL (:) = 0.0 
  DMT%XHVAC_COOL  (:) = 0.0
  DMT%XHVAC_HEAT  (:) = 0.0
  !
  DMT%XT_SYS      (:) = ZT_MIX(:)
  DMT%XQ_SYS      (:) = ZQ_MIX(:)
  DMT%XM_SYS      (:) = 0.0
  !
  DMT%XH_WASTE    (:) = 0.0
  DMT%XLE_WASTE   (:) = 0.0
  DMT%XFAN_POWER  (:) = 0.0
  !
  END WHERE
  !
  !
  !---------------------------------------------------------------------------------
  ! ENERGY DEMAND COMPUTED
  !################################################################################
  !
  ! EVOLUTION OF THE INTERNAL TEMPERATURE AND HUMIDITY
  !
  !
  ! Update of indoor temperature
  !
  B%XTI_BLD(:) = ( B%XTI_BLD(:)                           + &
    PTSTEP/(ZRHOI(:)*XCPD*T%XBLD_HEIGHT(:))               * &
    ( T%XWALL_O_BLD(:) * PCONV_WL_BLD(:)                  + &
    B%XGR(:) * PCONV_WIN_BLD(:)                           + &
    B%XMASS_O_BLD(:) * ZCONV_MA_BLD(:)                    + &
    PCONV_RF_BLD(:)                                       + &
    ZCONV_FL_BLD(:)                                       + &
    ZQIN(:)*(1.0-B%XQIN_FRAD(:))*(1.0-B%XQIN_FLAT(:)) )   + &
    PTSTEP/T%XBLD_HEIGHT(:)                               * &
    ((ZINF(:) + ZNAT_VENT(:)) * PT_CAN(:)                 + &
    DMT%XM_SYS(:) / ZRHOI(:)     * (DMT%XT_SYS(:)) ) )    / &
    (1. + PTSTEP/T%XBLD_HEIGHT(:)                         * &
    (ZINF(:) + ZNAT_VENT(:) + DMT%XM_SYS(:) / ZRHOI(:)) )
  !
  ! Update of indoor humidity
  !
  B%XQI_BLD(:) = ( B%XQI_BLD(:) +  PTSTEP/T%XBLD_HEIGHT(:)  * &
    (  ZQIN(:)*B%XQIN_FLAT(:)/(ZRHOI(:)*XLVTT)              + &
    (ZINF(:)+ZNAT_VENT(:))*PQ_CAN(:)                        + &
    DMT%XM_SYS(:)/ZRHOI(:)   *(DMT%XQ_SYS   (:)) ))         / &
    ( 1. + PTSTEP/T%XBLD_HEIGHT(:)                          * &
    (ZINF(:) + ZNAT_VENT(:) + DMT%XM_SYS(:)/ZRHOI(:)) )
  !
  ! Waste heat due to infiltration/ventilation
  !
  ZWASTE   (:) = (ZINF(:)+ZV_VENT(:)+ZNAT_VENT(:)) * ZRHOI(:) 
  DMT%XH_WASTE (:) = DMT%XH_WASTE (:) + ZWASTE(:) * XCPD  * (B%XTI_BLD(:) - PT_CAN(:))
  DMT%XLE_WASTE(:) = DMT%XLE_WASTE(:) + ZWASTE(:) * XLVTT * (B%XQI_BLD(:) - PQ_CAN(:)) 
  !
  IF (LHOOK) CALL DR_HOOK('BEM',1,ZHOOK_HANDLE)
  !
  CONTAINS
    !
    SUBROUTINE GET_NAT_VENT(PPTI_BLD,PPT_CANYON,PPU_CANYON,PPGLAZ_O_BLD, &
      PSHADE,PSHGC_SH,PNAT_VENT,PFOPEN)
      !
      IMPLICIT NONE
      !
      REAL, DIMENSION(:) , INTENT(IN)  :: PPTI_BLD
      REAL, DIMENSION(:) , INTENT(IN)  :: PPT_CANYON
      REAL, DIMENSION(:) , INTENT(IN)  :: PPU_CANYON
      REAL, DIMENSION(:) , INTENT(IN)  :: PPGLAZ_O_BLD
      REAL, DIMENSION(:) , INTENT(IN)  :: PSHGC_SH
      REAL, DIMENSION(:) , INTENT(IN)  :: PFOPEN
      REAL, DIMENSION(:) , INTENT(IN)  :: PSHADE
      REAL, DIMENSION(:) , INTENT(OUT) :: PNAT_VENT
      REAL(KIND=JPRB)   :: ZHOOK_HANDLE
      !
      ! Local variables
      !
      REAL :: ZHWIN,ZCT,ZCW,ZCST,ZALPHA,ZFAC
      REAL, DIMENSION(SIZE(PPTI_BLD)) :: ZMODVENT
      !
      IF (LHOOK) CALL DR_HOOK('BEM:GET_NAT_VENT',0,ZHOOK_HANDLE)
      !
      ! Robert: Ventilation flow rate through single sided windows
      ! taken from the European directive prEN 15242.
      ! "Ventilation for buildings — Calculation methods for the
      !  determination of air flow rates in buildings including infiltration"
      !
      ! The exchanges due to turbulence, wind force and buouyancy are considered
      !
      ZCT    = 0.01   ! Coefficient for exchange due to turbulence
      ZCW    = 0.001  ! Coefficient for exchange due to wind force
      ZCST   = 0.0035 ! Coefficient for exchange due to buouyancy
      ZALPHA = 180.0  ! Window opening angle (0° for closed, 180° for fully open)
      ZHWIN  = 1.5    ! Height of the window (m)
      !
      ! Empirical factor for dependency of air-flow on window opening angle
      !
      ZFAC = 2.60E-7*ZALPHA**3-1.19E-4*ZALPHA**2+1.86E-2*ZALPHA
      !
      ! Natural ventilation flow rate [m^3/m^2(bld)/s] (single sided ventilation)
      !
      PNAT_VENT(:)=ZFAC*0.5*PPGLAZ_O_BLD(:)*SQRT(ZCT+ZCW*PPU_CANYON(:)**2+ZCST*ZHWIN*ABS(PPTI_BLD(:)-PPT_CANYON(:)))
      !
      ! In the case of shading, the natural ventilation is influenced by
      ! the shading elements. However, the degree of modification is highly
      ! dependent on the building and window geometry, the type of shades, 
      ! the meteorological conditions and so on. A general formula
      ! is not available. For this reason, an ad-hoc modification of the 
      ! natural ventilation is made here based on the transmissivity of the
      ! shading elements.
      !
      ! No modification of the natural ventilation flow rate is made
      ! for shades with relatively high transmissivity.
      ! Indeed, the shades can even enhance ventilation.
      ! For shades with very low transmissivity it is assumed
      ! that only one-tenth of the natural ventilation takes place
      !
      DO JJ=1,SIZE(PSHGC_SH)
        IF (PSHGC_SH(JJ).GT.0.5) THEN
          ZMODVENT(JJ)=1.0
        ELSE
          ZMODVENT(JJ)=2*PSHGC_SH(JJ)
          ZMODVENT(JJ)=MAX(ZMODVENT(JJ),0.1)
        ENDIF
      ENDDO
      !
      ! The modulation factor is applied according
      ! to the fraction of closed blinds
      !
      PNAT_VENT(:)=PSHADE(:)*(ZMODVENT(:)*PNAT_VENT(:))+(1.0-PSHADE(:))*PNAT_VENT(:)
      !
      ! Multiplication of the ventilation flow rate
      ! by the fraction of windows opened
      !
      PNAT_VENT(:)=PNAT_VENT(:)*PFOPEN(:)
      !
      IF (LHOOK) CALL DR_HOOK('BEM:GET_NAT_VENT',1,ZHOOK_HANDLE)
      !
    END SUBROUTINE GET_NAT_VENT
    !
    SUBROUTINE GET_INFILTRATION(PRHOA,PPN50,PPTI_BLD,PPT_CANYON,PPU_CANYON,PPBLD_HEIGHT,PINFCALC)
      !
      USE MODD_CSTS, ONLY : XG,XRV,XRD
      !
      ! Robert: Calculation of the air exchange rate due to infiltration [1/h]
      ! as a funtion of the simulated canyon wind speed and the simulated 
      ! outdoor and indoor temperatures.
      ! The equations are similar to the Alberta Air Infiltration Model 
      ! (e.g. Wang et al. 2009, Building and Environment).
      !
      IMPLICIT NONE
      !
      REAL, INTENT(IN), DIMENSION(:) :: PRHOA
      REAL, INTENT(IN), DIMENSION(:) :: PPN50
      REAL, INTENT(IN), DIMENSION(:) :: PPTI_BLD
      REAL, INTENT(IN), DIMENSION(:) :: PPT_CANYON
      REAL, INTENT(IN), DIMENSION(:) :: PPU_CANYON
      REAL, INTENT(IN), DIMENSION(:) :: PPBLD_HEIGHT
      !
      REAL, INTENT(OUT), DIMENSION(:) :: PINFCALC
      !
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      !
      ! Local variables
      !
      REAL :: ZFLOW
      REAL :: ZBETA
      REAL :: ZR
      REAL :: ZX
      REAL :: ZY
      REAL :: ZHF
      REAL :: ZZF
      REAL :: ZM
      REAL :: ZXC
      REAL :: ZF
      REAL :: ZJ
      !
      REAL, DIMENSION(SIZE(PPTI_BLD)) :: ZCAIRTIGHT
      REAL, DIMENSION(SIZE(PPTI_BLD)) :: ZRHOOUT
      REAL, DIMENSION(SIZE(PPTI_BLD)) :: ZFSTACK
      REAL, DIMENSION(SIZE(PPTI_BLD)) :: ZFWIND
      REAL, DIMENSION(SIZE(PPTI_BLD)) :: ZQSTACK
      REAL, DIMENSION(SIZE(PPTI_BLD)) :: ZQWIND
      REAL, DIMENSION(SIZE(PPTI_BLD)) :: ZQTOTAL
      !
      IF (LHOOK) CALL DR_HOOK('BEM:GET_INFILTRATION',0,ZHOOK_HANDLE)
      !
      ZFLOW = 2./3. ! Flow exponent (1.0 = lam.; 0.5 = turb.; 2/3 for typical infil. flow)
      ZBETA = -0.33 ! Empirical constant
      ZR    = 0.0   ! All leakage concentrated in the walls
      ZX    = 0.0   ! No sealing-floor difference since all leakage in walls
      ZY    = 0.2   ! Typical value for flue fraction (Walker and Wilson, 1998)
      ZZF   = 1.5   ! Typical value for normalised flue height (Walker and Wilson, 1998)
      !
      ! Stack factor based on Walker and Wilson (1998)
      !
      ZM=MIN(1.0,(ZX+(2*ZFLOW+1)*ZY)**2/(2.0-ZR))
      !
      ZXC=ZR+2.0*(1.0-ZR-ZY)/(ZFLOW+1.0)-2.0*ZY*(ZZF-1.0)**ZFLOW
      !
      ZF=ZFLOW*ZY*(ZZF-1)**(ZFLOW-1./3.)*(1.0-(3*(ZXC-ZX)**2*(ZR)**(1.0-ZFLOW))/(2*(ZZF+1)))
      !
      ZFSTACK=((1+ZFLOW*ZR)/(ZFLOW+1))*(0.5-0.5*(ZM)**(5./4.))**(ZFLOW+1)+ZF
      !
      ! Infiltration factor based on Walker and Wilson (1998), no crawl space
      !
      ZJ=0.5*(ZX+ZR+2*ZY)
      !
      ZFWIND=0.19*(2.0-ZFLOW)*(1.0-(0.5*(ZX+ZR))**(1.5-ZY))-0.25*ZY*(ZJ-2.0*ZY*(ZJ)**4)
      !
      ! Conversion of PPN50 [1/h] to CAIRTIGHT [1/h]
      !
      ZCAIRTIGHT(:)=PPN50(:)/(50)**ZFLOW
      !
      ! The infiltration rate due to the stack effect
      !
      ZQSTACK(:)=ZCAIRTIGHT(:)*ZFSTACK*(PRHOA(:)*XG*PPBLD_HEIGHT(:)*(ABS(PPTI_BLD(:)-PPT_CANYON(:))/PPTI_BLD(:)))**ZFLOW
      !
      ! The infiltration rate due to the wind effect
      ! A shelter factor of 1.0 is assumed since the canyon wind speed is used.
      !
      ZQWIND(:)=ZCAIRTIGHT(:)*ZFWIND*(PRHOA(:)*0.5*(PPU_CANYON(:))**2)**ZFLOW
      !
      ! Total infiltration rate in 1/h
      !
      PINFCALC(:)=((ZQSTACK(:))**(1.0/ZFLOW)+(ZQWIND(:))**(1.0/ZFLOW)+ZBETA*(ZQSTACK(:)*ZQWIND(:))**(1.0/(2.0*ZFLOW)))**ZFLOW
      !
      IF (LHOOK) CALL DR_HOOK('BEM:GET_INFILTRATION',1,ZHOOK_HANDLE)
    !
    END SUBROUTINE GET_INFILTRATION
!
END SUBROUTINE BEM
