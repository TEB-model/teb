!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
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
!     Awall   =  WALL_O_HOR * (1 - GR) / BLD [m2_wall/m2_bld]
!     Awin    =  WALL_O_HOR * GR / BLD       [m2_win/m2_bld]   
!     Amass   =  2 * N_FLOOR                  [m2_mass/m2_bld]  
!     N_FLOOR  =  BLD_HEIGHT / FLOOR_HEIGHT   [#]
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
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZFAN_AP      ! Fan design pressure increase [Pa]
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZFAN_EFF     ! Fan total efficiency
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
LOGICAL, DIMENSION(SIZE(B%XTI_BLD)):: GNAT_VENT  ! Is Natural ventilation active ? 
REAL, DIMENSION(SIZE(B%XTI_BLD)):: ZNAT_VENT     ! Nat.vent airflow rate [m3 s-1 m-2(bld)]
REAL,DIMENSION(SIZE(B%XTI_BLD)) :: ZTI_BLD       ! Indoor air temperature at time step t + dt [K]
REAL,DIMENSION(SIZE(B%XTI_BLD)) :: ZTI_BLD_OPEN  ! Indoor air temperature if windows opened
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
INTEGER :: JJ                                  ! Loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!!REAL :: ZEXPL = 0.5 !explicit coefficient for internal temperature evol.
!!REAL :: ZIMPL = 0.5 !implicit coef..
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('BEM',0,ZHOOK_HANDLE)
!
!*      1.   Initializations
!            ---------------
!
ZRHOI(:) = PPS(:) / (XRD * B%XTI_BLD(:) * ( 1.+((XRV/XRD)-1.)*B%XQI_BLD(:) ) )
! *Temperal definitions for nat.vent*
!ZTCOMF_MAX  = 26. + 273.16
!
! *Definitions
ZFAN_AP (:) = 600.0
ZFAN_EFF(:) = 0.7
!
! *Other calcs
IRF  = SIZE(T%XT_ROOF,2)
IWL  = SIZE(T%XT_WALL_A,2)
!
!
! initial condition of QI_BLD equivalent to 50% RH
IF (ANY(B%XQI_BLD(:) <= 1E-6)) B%XQI_BLD = 0.5 * QSAT(B%XTI_BLD, PPS)
!
! *Temperal definitions for shedule*
GSCHED(:) = .FALSE.
WHERE (GSCHED(:))
 ZF_NIGHT(:)  = 0.8
 ZF_DAY(:)    = 1.2
ELSE WHERE
 ZF_NIGHT(:)  = 1.
 ZF_DAY(:)    = 1.
END WHERE
!
! *Int.gains schedule
!
ZQIN = DMT%XQIN * B%XN_FLOOR
WHERE (PSUNTIME(:) > 0. .AND. PSUNTIME(:) < 25200.) ! night between 0000 and 0700
  ZQIN(:) = ZQIN(:) * ZF_NIGHT(:)
ELSEWHERE
  ZQIN(:) = ZQIN(:) * ZF_DAY(:)
END WHERE

! *Change of units AC/H -> [m3 s-1 m-2(bld)]
ZV_VENT(:) = B%XV_VENT(:) * T%XBLD_HEIGHT(:) / 3600.
ZINF   (:) = B%XINF   (:) * T%XBLD_HEIGHT(:) / 3600.  
!
!*      2.   heat balance for building floor and mass
!            ----------------------------------------
!
!*      2.1 total load on the internal mass or floor
ZLOAD_FL(:) = (ZQIN(:) * B%XQIN_FRAD(:) * (1.-B%XQIN_FLAT(:)) + DMT%XTR_SW_WIN(:)) / (B%XMASS_O_BLD(:)+1.)
WHERE (B%XN_FLOOR(:) > 1.)
   ZLOAD_MA(:) = ZLOAD_FL(:)
ELSEWHERE
   ZLOAD_MA(:) = 0.
ENDWHERE
!
!*      2.2 FLOOR HEAT BALANCE
!
 CALL FLOOR_LAYER_E_BUDGET(B, PTSTEP, PFLX_BLD_FL, ZDQS_FL, ZIMB_FL, PRADHT_IN, &
                           PRAD_WL_FL, PRAD_RF_FL, PRAD_WIN_FL, PLOAD_IN_FL,    &
                           ZRAD_FL_MA, ZCONV_FL_BLD)
!
!*      2.3 MASS HEAT BALANCE
!
 CALL MASS_LAYER_E_BUDGET(B, PTSTEP, PFLX_BLD_MA, ZDQS_MA, ZIMB_MA, PRADHT_IN, &
                          PRAD_WL_MA, PRAD_RF_MA, PRAD_WIN_MA, PLOAD_IN_MA,    &
                          ZRAD_FL_MA, ZCONV_MA_BLD  )
!
!
ZAC_IN_WL_A_COOL = CHTC_VERT_DOE(T%XT_WALL_A(:,IWL), DMT%XTCOOL_TARGET)
ZAC_IN_WL_B_COOL = CHTC_VERT_DOE(T%XT_WALL_B(:,IWL), DMT%XTCOOL_TARGET)
ZAC_IN_WIN_COOL  = CHTC_VERT_DOE(B%XT_WIN2         , DMT%XTCOOL_TARGET)
ZAC_IN_MA_COOL   = CHTC_VERT_DOE(B%XT_MASS(:,1)    , DMT%XTCOOL_TARGET)
ZAC_IN_RF_COOL   = CHTC_DOWN_DOE(T%XT_ROOF(:,IRF)  , DMT%XTCOOL_TARGET)
ZAC_IN_FL_COOL   = CHTC_UP_DOE  (B%XT_FLOOR(:,1)   , DMT%XTCOOL_TARGET)

ZAC_IN_WL_A_HEAT = CHTC_VERT_DOE(T%XT_WALL_A(:,IWL), DMT%XTHEAT_TARGET)
ZAC_IN_WL_B_HEAT = CHTC_VERT_DOE(T%XT_WALL_B(:,IWL), DMT%XTHEAT_TARGET)
ZAC_IN_WIN_HEAT  = CHTC_VERT_DOE(B%XT_WIN2         , DMT%XTHEAT_TARGET)
ZAC_IN_MA_HEAT   = CHTC_VERT_DOE(B%XT_MASS(:,1)    , DMT%XTHEAT_TARGET)
ZAC_IN_RF_HEAT   = CHTC_DOWN_DOE(T%XT_ROOF(:,IRF)  , DMT%XTHEAT_TARGET)
ZAC_IN_FL_HEAT   = CHTC_UP_DOE  (B%XT_FLOOR(:,1)   , DMT%XTHEAT_TARGET)

DO JJ=1,SIZE(ZAC_IN_WIN_COOL)
   ZAC_IN_WL_A_COOL(JJ) = MAX(1., ZAC_IN_WL_A_COOL(JJ))
   ZAC_IN_WL_B_COOL(JJ) = MAX(1., ZAC_IN_WL_B_COOL(JJ))
   ZAC_IN_WIN_COOL (JJ) = MAX(1., ZAC_IN_WIN_COOL(JJ))
   ZAC_IN_MA_COOL  (JJ) = MAX(1., ZAC_IN_MA_COOL (JJ))
   ZAC_IN_RF_COOL  (JJ) = MAX(1., ZAC_IN_RF_COOL (JJ))
   ZAC_IN_FL_COOL  (JJ) = MAX(1., ZAC_IN_FL_COOL (JJ))
   
   ZAC_IN_WL_A_HEAT(JJ) = MAX(1., ZAC_IN_WL_A_HEAT(JJ))
   ZAC_IN_WL_B_HEAT(JJ) = MAX(1., ZAC_IN_WL_B_HEAT(JJ))
   ZAC_IN_WIN_HEAT (JJ) = MAX(1., ZAC_IN_WIN_HEAT(JJ))
   ZAC_IN_MA_HEAT  (JJ) = MAX(1., ZAC_IN_MA_HEAT (JJ))
   ZAC_IN_RF_HEAT  (JJ) = MAX(1., ZAC_IN_RF_HEAT (JJ))
   ZAC_IN_FL_HEAT  (JJ) = MAX(1., ZAC_IN_FL_HEAT (JJ))
ENDDO

!*      4.   Indoor energy balance calculation
!            ---------------------------------
!
DO JJ=1,SIZE(PT_CAN)
  ! *first guess of indoor temperature

  ZTI_BLD(JJ) = B%XTI_BLD(JJ) + PTSTEP/(ZRHOI(JJ) * XCPD * T%XBLD_HEIGHT(JJ))    & 
          * (  T%XWALL_O_BLD(JJ) * PCONV_WL_BLD(JJ) + B%XGLAZ_O_BLD (JJ) * PCONV_WIN_BLD(JJ)  &
             + B%XMASS_O_BLD(JJ) * ZCONV_MA_BLD(JJ) + PCONV_RF_BLD(JJ) + ZCONV_FL_BLD(JJ)   &
             + ZQIN(JJ) * (1 - B%XQIN_FRAD(JJ))  * (1 - B%XQIN_FLAT(JJ)) )
  !
  !################################################################################
  ! *is natural surventilation active at the current time step ?
  !---------------------------------------------------------------------------------
  !
  !    *no surventilation possible

  IF (B%CNATVENT(JJ)=='NONE') THEN
    !
    GNAT_VENT(JJ) = .FALSE.
    !
  !    *automatic management of surventilation
  ELSEIF (B%CNATVENT(JJ)=='AUTO' .OR. B%CNATVENT(JJ)=='MECH') THEN
    !
    IF (MOD(PSUNTIME(JJ), 3600.) .LT. PTSTEP) THEN
      !
      IF ( B%XTI_BLD(JJ).GT. PT_CAN(JJ) + 1 ) THEN ! condition to enable the
        IF (B%CNATVENT(JJ)=='AUTO') THEN
        ! natural surventilation rate calculation (window opening)
          CALL GET_NAT_VENT(B%XTI_BLD(JJ), PT_CAN(JJ), PU_CAN(JJ), B%XGR(JJ), &
                            B%XFLOOR_HW_RATIO(JJ), T%XBLD_HEIGHT(JJ), ZNAT_VENT(JJ))
        ELSE IF (B%CNATVENT(JJ)=='MECH') THEN
        ! mechanical surventilation rate calculation : 5 volumes/hour
          ZNAT_VENT(JJ) =  5.0*T%XBLD_HEIGHT(JJ)/3600.
        END IF
        !
        ZTI_BLD_OPEN  (JJ) = ZTI_BLD(JJ) &
                + ZNAT_VENT(JJ)            * PTSTEP/T%XBLD_HEIGHT(JJ) * (PT_CAN(JJ) - B%XTI_BLD(JJ)) 
        ZTI_BLD_CLOSED(JJ) = ZTI_BLD(JJ) &
                + (ZINF(JJ) + ZV_VENT(JJ)) * PTSTEP/T%XBLD_HEIGHT(JJ) * (PT_CAN(JJ) - B%XTI_BLD(JJ)) 
        !
        GNAT_VENT(JJ) = (ZTI_BLD_OPEN(JJ) <= DMT%XTCOOL_TARGET (JJ) .AND. &            
                         ZTI_BLD_OPEN(JJ) <  ZTI_BLD_CLOSED(JJ) .AND. &
                         ZTI_BLD_OPEN(JJ) >  DMT%XTHEAT_TARGET (JJ) + 4.)
        !
      ELSE
        GNAT_VENT(JJ) = .FALSE.
      ENDIF
      B%LNATVENT_NIGHT(JJ) = GNAT_VENT(JJ)
    ELSE 
      GNAT_VENT(JJ) = B%LNATVENT_NIGHT(JJ)
    ENDIF
    !
  !    *manual management of surventilation
  ELSEIF (B%CNATVENT(JJ)=='MANU') THEN
    !
    B%LNATVENT_NIGHT(JJ) = B%LNATVENT_NIGHT(JJ) .AND. &
                         .NOT. ( PSUNTIME(JJ) > 5.*3600 .AND. PSUNTIME(JJ) < 18.*3600 )
    !
    GNAT_VENT(JJ) = ( PSUNTIME(JJ) > 18.*3600. .AND. PSUNTIME(JJ) < 21.*3600.  &
                      .AND. PT_CAN(JJ) < B%XTI_BLD(JJ)+2.       &
                      .AND. PT_CAN(JJ) > DMT%XTHEAT_TARGET(JJ)    & 
                      .AND. ( B%XTI_BLD(JJ) > DMT%XTHEAT_TARGET(JJ)+5. &
                       .OR. B%XTI_BLD(JJ) == DMT%XTCOOL_TARGET(JJ) )   ) 
    GNAT_VENT(JJ) = GNAT_VENT(JJ) .OR. B%LNATVENT_NIGHT(JJ)
    !
  ENDIF
  !
  ! Decicion about natural surventilation OK
  !################################################################################
  !
  !
  !################################################################################
  ! COMPUTE ENERGY DEMAND
  !---------------------------------------------------------------------------------

  ! *If natural surventilation ACTIVE
  IF (GNAT_VENT(JJ)) THEN
     ! 
     CALL GET_NAT_VENT(B%XTI_BLD(JJ), PT_CAN(JJ), PU_CAN(JJ), B%XGR(JJ), &
                       B%XFLOOR_HW_RATIO(JJ), T%XBLD_HEIGHT(JJ), ZNAT_VENT(JJ)     )
     !
     ZV_VENT     (JJ) = 0.
     ZINF        (JJ) = 0.
     !
     DMT%XH_BLD_COOL (JJ) = 0.0         ! No HVAC consumption
     DMT%XH_BLD_HEAT (JJ) = 0.0    
     DMT%XLE_BLD_COOL(JJ) = 0.0         ! No HVAC consumption
     DMT%XLE_BLD_HEAT(JJ) = 0.0   
     !    
     DMT%XT_BLD_COOL (JJ) = 0.0         ! No HVAC consumption
     DMT%XHVAC_COOL  (JJ) = 0.0    
     DMT%XT_SYS      (JJ) = B%XTI_BLD(JJ) ! No mechanical ventilation
     DMT%XQ_SYS      (JJ) = B%XQI_BLD(JJ) ! 
     DMT%XH_WASTE    (JJ) = 0.0
     DMT%XLE_WASTE   (JJ) = 0.0     
     DMT%XFAN_POWER  (JJ) = 0.0    
     DMT%XHVAC_HEAT  (JJ) = 0.0
     !
     DMT%XM_SYS  (JJ) = 0.0
     DMT%XCOP    (JJ) = 0.0
     DMT%XCAP_SYS(JJ) = 0.0
     !
  ! *If natural surventilation INACTIVE
  ELSE 
    !
    ZNAT_VENT(JJ) = 0.
    !
    ! ------------------------------------------------
    ! * Building energy demand for heating and cooling
    ! ------------------------------------------------
    !
    DMT%XH_BLD_COOL(JJ) = T%XWALL_O_BLD(JJ)/2. * (ZAC_IN_WL_A_COOL(JJ) * (T%XT_WALL_A(JJ,IWL) - DMT%XTCOOL_TARGET(JJ))  &
                                            + ZAC_IN_WL_B_COOL(JJ) * (T%XT_WALL_B(JJ,IWL) - DMT%XTCOOL_TARGET(JJ))) &
                         + B%XGLAZ_O_BLD(JJ) * ZAC_IN_WIN_COOL(JJ) * (B%XT_WIN2(JJ)       - DMT%XTCOOL_TARGET(JJ))  &    
                          + ZAC_IN_MA_COOL (JJ)* B%XMASS_O_BLD(JJ) * (B%XT_MASS(JJ,1)     - DMT%XTCOOL_TARGET(JJ))  &
                                              + ZAC_IN_RF_COOL(JJ) * (T%XT_ROOF(JJ,IRF)   - DMT%XTCOOL_TARGET(JJ))  &
                                              + ZAC_IN_FL_COOL(JJ) * (B%XT_FLOOR(JJ,1)    - DMT%XTCOOL_TARGET(JJ))  &
                                 + ZQIN(JJ) * (1 - B%XQIN_FRAD(JJ)) * (1 - B%XQIN_FLAT(JJ))                   &
                              + (ZINF(JJ) + ZV_VENT(JJ)) * ZRHOI(JJ) * XCPD * (PT_CAN(JJ) - DMT%XTCOOL_TARGET(JJ))
    !
    DMT%XH_BLD_HEAT(JJ) = - ( T%XWALL_O_BLD(JJ)/2. * (ZAC_IN_WL_A_HEAT(JJ) * (T%XT_WALL_A(JJ,IWL) - DMT%XTHEAT_TARGET(JJ))  &
                                                + ZAC_IN_WL_B_HEAT(JJ) * (T%XT_WALL_B(JJ,IWL) - DMT%XTHEAT_TARGET(JJ))) &
                             + B%XGLAZ_O_BLD(JJ) * ZAC_IN_WIN_HEAT(JJ) * (B%XT_WIN2(JJ)       - DMT%XTHEAT_TARGET(JJ))  &    
                              +  ZAC_IN_MA_HEAT(JJ)* B%XMASS_O_BLD(JJ) * (B%XT_MASS(JJ,1)     - DMT%XTHEAT_TARGET(JJ))  &
                                                  + ZAC_IN_RF_HEAT(JJ) * (T%XT_ROOF(JJ,IRF)   - DMT%XTHEAT_TARGET(JJ))  &
                                                  + ZAC_IN_FL_HEAT(JJ) * (B%XT_FLOOR(JJ,1)    - DMT%XTHEAT_TARGET(JJ))  &
                                   + ZQIN(JJ) * (1 - B%XQIN_FRAD(JJ))* (1 - B%XQIN_FLAT(JJ))                      &
                                  + (ZINF(JJ) + ZV_VENT(JJ)) * ZRHOI(JJ) * XCPD * (PT_CAN(JJ) - DMT%XTHEAT_TARGET(JJ)))
    !
    ZQCOOL_TRGT(JJ) = 0.62198 * B%XHR_TARGET(JJ) * PSAT(DMT%XTCOOL_TARGET(JJ)) / &
                      (PPS(JJ)- B%XHR_TARGET(JJ) * PSAT(DMT%XTCOOL_TARGET(JJ)))    
    !
    DMT%XLE_BLD_COOL(JJ) = ZQIN(JJ) * B%XQIN_FLAT(JJ)                                           &
               + (ZINF(JJ) + ZV_VENT(JJ)) * ZRHOI(JJ) * XLVTT * (PQ_CAN(JJ) - ZQCOOL_TRGT(JJ)) 
    !
    ZQHEAT_TRGT(JJ) = 0.62198 * B%XHR_TARGET(JJ) * PSAT(DMT%XTHEAT_TARGET(JJ)) / &
                      (PPS(JJ)- B%XHR_TARGET(JJ) * PSAT(DMT%XTHEAT_TARGET(JJ)))    
    !
    DMT%XLE_BLD_HEAT(JJ) = ZQIN(JJ) * B%XQIN_FLAT(JJ)                                           &
              + (ZINF(JJ) + ZV_VENT(JJ)) * ZRHOI(JJ) * XLVTT * (PQ_CAN(JJ) - ZQHEAT_TRGT(JJ))       
    !
    ! * Autosize calculations
    !
    IF (BOP%LAUTOSIZE .AND. KDAY==15) THEN
      !
      IF (DMT%XH_BLD_COOL(JJ) > B%XAUX_MAX(JJ))  THEN
        !  
        B%XAUX_MAX    (JJ) = DMT%XH_BLD_COOL(JJ)
        !
        ! Cooling coil sensible heat rate 
        ZSHR        (JJ) = MIN(XCPD * (DMT%XTCOOL_TARGET(JJ) - B%XT_ADP(JJ)) /           &
                                  (ENTH_FN_T_Q(DMT%XTCOOL_TARGET(JJ),ZQCOOL_TRGT(JJ)) -  &
                                   ENTH_FN_T_Q(B%XT_ADP(JJ),QSAT(B%XT_ADP(JJ),PPS(JJ)))), 1.)
        ! Cooling Coil Capacity [W m-2(bld)]
        B%XCAP_SYS_RAT(JJ) = DMT%XH_BLD_COOL(JJ) / ZSHR(JJ) 
        !
        ! Cooling rated air flow rate [kg s-1 m-2(bld)]
        ZM_SYS_RAT  (JJ) = DMT%XH_BLD_COOL(JJ) / XCPD / (DMT%XTCOOL_TARGET(JJ)-(14.0+273.16))
        IF (ZM_SYS_RAT(JJ) > B%XM_SYS_RAT(JJ)) B%XM_SYS_RAT(JJ) = ZM_SYS_RAT(JJ)
        !
        ! Impose condition 
        IF (B%XM_SYS_RAT(JJ)/ZRHOI(JJ)/B%XCAP_SYS_RAT(JJ) < 0.00004027) THEN
          B%XCAP_SYS_RAT(JJ) = B%XM_SYS_RAT(JJ)/ZRHOI(JJ)/0.00004027
        ELSE IF (B%XM_SYS_RAT(JJ)/ZRHOI(JJ)/B%XCAP_SYS_RAT(JJ) > 0.00006041) THEN
          B%XCAP_SYS_RAT(JJ) = B%XM_SYS_RAT(JJ)/ZRHOI(JJ)/0.00006041
        END IF
        !
      END IF
      !
    END IF
    !
    ! * END Autosize calculations
    !
    ! * system efficiency
    ! ...................
    !
    DMT%XM_SYS  (JJ) = B%XM_SYS_RAT  (JJ)
    DMT%XCOP    (JJ) = B%XCOP_RAT    (JJ)
    DMT%XCAP_SYS(JJ) = B%XCAP_SYS_RAT(JJ)
    !
    ! * Mixing conditions
    ! .................
    !
    ZXMIX (JJ) = ZV_VENT(JJ) * ZRHOI(JJ) / DMT%XM_SYS(JJ)
    ZT_MIX(JJ) = ZXMIX(JJ) * PT_CAN(JJ) + (1.-ZXMIX(JJ)) * B%XTI_BLD(JJ)
    ZQ_MIX(JJ) = ZXMIX(JJ) * PQ_CAN(JJ) + (1.-ZXMIX(JJ)) * B%XQI_BLD(JJ)
    ! 
    ! ---------------------------------------------
    ! * COOLING system : Performance and Waste heat
    ! ---------------------------------------------
    !
    IF (DMT%XH_BLD_COOL(JJ) >= 0.0) THEN
      !
      ! *ideal system
      IF (BOP%CCOOL_COIL=='IDEAL') THEN
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
          !
        ! *real system
        ELSEIF (BOP%CCOOL_COIL=='DXCOIL') THEN
          !
          CALL DX_AIR_COOLING_COIL_CV(PT_CAN(JJ), PQ_CAN(JJ), PPS(JJ),  ZRHOI(JJ), ZT_MIX(JJ), &
                                      ZQ_MIX(JJ), B%XCOP_RAT(JJ), B%XCAP_SYS_RAT(JJ),          &
                                      B%XT_ADP(JJ), B%XF_WATER_COND(JJ), DMT%XM_SYS(JJ),      &
                                      DMT%XH_BLD_COOL(JJ), DMT%XH_WASTE(JJ), DMT%XLE_WASTE(JJ), &
                                      DMT%XCOP(JJ), DMT%XCAP_SYS(JJ), DMT%XT_SYS(JJ), & 
                                      DMT%XQ_SYS(JJ), DMT%XHVAC_COOL(JJ), DMT%XT_BLD_COOL(JJ) )
          !
        ENDIF !end type of cooling system

        !!! case of system without atmospheric releases. I-e releases in soil/water F_WATER_COND < 0 
        IF (B%XF_WATER_COND(JJ) < 0) THEN
          DMT%XH_WASTE(JJ) = 0. 
          DMT%XLE_WASTE(JJ) = 0. 
        ENDIF
        !!!!
        !
        !         From EP Engineering Reference (p. 647)
        DMT%XFAN_POWER(JJ) = DMT%XM_SYS(JJ) * ZFAN_AP(JJ) * ZFAN_EFF(JJ) * ZRHOI(JJ)
        !
        DMT%XH_BLD_HEAT (JJ) = 0.0
        DMT%XLE_BLD_HEAT(JJ) = 0.0
        DMT%XHVAC_HEAT  (JJ) = 0.0
      !
      ! ---------------------------------------------
      ! * HEATING system : Performance and Waste heat
      ! ---------------------------------------------
      !
      ELSE IF (DMT%XH_BLD_HEAT(JJ) > 0.0) THEN
        !
        ! *specific computation for real heating system
        IF  (BOP%CHEAT_COIL .EQ. 'FINCAP') THEN
          IF (DMT%XH_BLD_HEAT(JJ) > B%XCAP_SYS_HEAT(JJ)) DMT%XH_BLD_HEAT(JJ) =  B%XCAP_SYS_HEAT(JJ)
        END IF
        !
        DMT%XT_SYS(JJ) = ZT_MIX(JJ) + DMT%XH_BLD_HEAT(JJ) / DMT%XM_SYS(JJ) / XCPD
        DMT%XQ_SYS(JJ) = ZQ_MIX(JJ)
        !
        DMT%XHVAC_HEAT  (JJ) = DMT%XH_BLD_HEAT(JJ) / B%XEFF_HEAT(JJ)
        DMT%XH_WASTE    (JJ) = DMT%XHVAC_HEAT(JJ) - DMT%XH_BLD_HEAT(JJ)  
        DMT%XLE_WASTE   (JJ) = 0.0
        DMT%XH_BLD_COOL (JJ) = 0.0
        DMT%XLE_BLD_COOL(JJ) = 0.0
        DMT%XT_BLD_COOL (JJ) = 0.0
        DMT%XHVAC_COOL  (JJ) = 0.0
!       From EP Engineering Reference (p. 647)
        DMT%XFAN_POWER(JJ) = DMT%XM_SYS(JJ)*ZFAN_AP(JJ)*(ZFAN_EFF(JJ)*ZRHOI(JJ))
     !
     ! ------------------------------
     ! * NEITHEIR COOLING NOR HEATING
     ! ------------------------------
     ! 
     ELSE
        !
        DMT%XH_BLD_COOL (JJ) = 0.0
        DMT%XH_BLD_HEAT (JJ) = 0.0 
        DMT%XLE_BLD_COOL(JJ) = 0.0
        DMT%XLE_BLD_HEAT(JJ) = 0.0
        !
        DMT%XT_BLD_COOL (JJ) = 0.0 
        DMT%XHVAC_COOL  (JJ) = 0.0
        DMT%XT_SYS      (JJ) = ZT_MIX(JJ)
        DMT%XQ_SYS      (JJ) = ZQ_MIX(JJ)
        DMT%XH_WASTE    (JJ) = 0.0
        DMT%XLE_WASTE   (JJ) = 0.0
        DMT%XFAN_POWER  (JJ) = 0.0     
        DMT%XHVAC_HEAT  (JJ) = 0.0
        !
     END IF !end for heating/cooling sytem
     !
  END IF
  !
  !---------------------------------------------------------------------------------
  ! ENERGY DEMAND COMPUTED
  !################################################################################
ENDDO
!
!---------------------------------------------------
! EVOLUTION OF THE INTERNAL TEMPERATURE AND HUMIDITY
!###################################################
!
ZTI_BLD(:) = ( ZTI_BLD(:) + PTSTEP/T%XBLD_HEIGHT(:) *                   & 
            ((ZINF(:) + ZNAT_VENT(:)) * PT_CAN(:) + DMT%XM_SYS(:) / ZRHOI(:) * (DMT%XT_SYS(:) ) )) &
          / (1. + PTSTEP/T%XBLD_HEIGHT(:)*(ZINF(:) + ZNAT_VENT(:) + DMT%XM_SYS(:) / ZRHOI(:))  )
ZQI_BLD(:) = ( B%XQI_BLD(:) +  PTSTEP/T%XBLD_HEIGHT(:) *                    & 
             (ZQIN(:) * B%XQIN_FLAT(:) / (ZRHOI(:) * XLVTT) + (ZINF(:) + ZNAT_VENT(:)) * (PQ_CAN(:))   &
              + DMT%XM_SYS(:) / ZRHOI(:)     * (DMT%XQ_SYS(:) ) ))&
          / (1. + PTSTEP/T%XBLD_HEIGHT(:)* (ZINF(:) + ZNAT_VENT(:) + DMT%XM_SYS(:) / ZRHOI(:))  )
!
! Update variables
B%XTI_BLD(:) = ZTI_BLD(:)
B%XQI_BLD(:) = ZQI_BLD(:)
!
! Waste heat due to infiltration/ventilation
ZWASTE   (:) = (ZINF(:)+ZV_VENT(:)+ZNAT_VENT(:)) * ZRHOI(:) 
DMT%XH_WASTE (:) = DMT%XH_WASTE (:) + ZWASTE(:) * XCPD  * (B%XTI_BLD(:) - PT_CAN(:))
DMT%XLE_WASTE(:) = DMT%XLE_WASTE(:) + ZWASTE(:) * XLVTT * (B%XQI_BLD(:) - PQ_CAN(:)) 
!
!
IF (LHOOK) CALL DR_HOOK('BEM',1,ZHOOK_HANDLE)
!
CONTAINS
!
SUBROUTINE GET_NAT_VENT(PTI_BLD, PPT_CAN, PPU_CAN, PGR, PF_AUX, PPBLD_HEIGHT, PNAT_VENT)
!
IMPLICIT NONE
!
REAL, INTENT(IN) :: PTI_BLD
REAL, INTENT(IN) :: PPT_CAN
REAL, INTENT(IN) :: PPU_CAN
REAL, INTENT(IN) :: PGR
REAL, INTENT(IN) :: PF_AUX
REAL, INTENT(IN) :: PPBLD_HEIGHT
REAL, INTENT(OUT) :: PNAT_VENT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('BEM:GET_NAT_VENT',0,ZHOOK_HANDLE)
!
PNAT_VENT = XG * (PTI_BLD - PPT_CAN)
IF (PNAT_VENT .LT. 0.) THEN ! exceptional case with MANU ventilation system
   PNAT_VENT= PPBLD_HEIGHT/3600. !minimum value
ELSE
   PNAT_VENT = 1./3. * (PNAT_VENT/PPT_CAN)**(1./2.)                          &
               * (1.5 + PTI_BLD/PNAT_VENT * 1./2. * PPU_CAN**2*0.1)**(3./2.) &
               * PGR * PF_AUX / 1.5 / 2.
   PNAT_VENT = MIN(PNAT_VENT, 5.0*PPBLD_HEIGHT/3600.)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('BEM:GET_NAT_VENT',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_NAT_VENT
!
END SUBROUTINE BEM
