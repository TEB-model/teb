!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE SNOW_COVER_1LAYER(PTSTEP, PANSMIN, PANSMAX, PTODRY, PRHOSMIN, PRHOSMAX,   &
                                 PRHOFOLD, OALL_MELT, PDRAIN_TIME, PWCRN, PZ0SN, PZ0HSN, &
                                 TPSNOW, PTG, PTG_COEFA, PTG_COEFB, PABS_SW, PLW1, PLW2, &
                                 PTA, PQA, PVMOD, PPS, PRHOA, PSR, PZREF, PUREF, PRNSNOW,&
                                 PHSNOW, PLESNOW, PGSNOW, PMELT, PDQS_SNOW, PABS_LW  )  
!   ##########################################################################
!
!!****  *SNOW_COVER_1LAYER*  
!!
!!    PURPOSE
!!    -------
!
!     One layer snow mantel scheme
!         
!     
!!**  METHOD
!     ------
!
!
! The temperature equation is written as:
!
!              b T+ = y
!
!
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
!!      Original    08/09/98 
!!      J. Escobar 24/10/2012 : BUF PGI10.X , rewrite some 1 line WHERE statement
!!      V. Masson  13/09/2013 : implicitation of coupling with roof below
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TYPE_SNOW, ONLY : SURF_SNOW
!
USE MODD_CSTS,       ONLY : XTT, XCI, XRHOLI, XRHOLW, XCPD, XLSTT, XLMTT, XDAY, XCONDI
USE MODD_SNOW_PAR,   ONLY : XEMISSN
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODE_THERMOS
!
USE MODI_SURFACE_RI
USE MODI_SURFACE_AERO_COND
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
REAL,                 INTENT(IN)    :: PTSTEP   ! time step
REAL,                 INTENT(IN)    :: PANSMIN  ! minimum snow albedo
REAL,                 INTENT(IN)    :: PANSMAX  ! maximum snow albedo
REAL,                 INTENT(IN)    :: PTODRY   ! snow albedo decreasing constant
REAL,                 INTENT(IN)    :: PRHOSMIN ! minimum snow density
REAL,                 INTENT(IN)    :: PRHOSMAX ! maximum snow density
REAL,                 INTENT(IN)    :: PRHOFOLD ! snow density increasing constant
LOGICAL,              INTENT(IN)    :: OALL_MELT! T --> all snow runs off if
                                                ! lower surf. temperature is
                                                ! positive
REAL,                 INTENT(IN)    :: PDRAIN_TIME ! drainage folding time (days)
REAL,                 INTENT(IN)    :: PWCRN    ! critical snow amount necessary
                                                ! to cover the considered surface
REAL,                 INTENT(IN)    :: PZ0SN    ! snow roughness length for momentum
REAL,                 INTENT(IN)    :: PZ0HSN   ! snow roughness length for heat
TYPE(SURF_SNOW), INTENT(INOUT) :: TPSNOW
REAL, DIMENSION(:), INTENT(IN)    :: PTG      ! underlying ground temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTG_COEFA! underlying ground temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTG_COEFB! implicit terms
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW  ! absorbed SW energy (Wm-2)
REAL, DIMENSION(:), INTENT(IN)    :: PLW1     ! LW coef independant of TSNOW
                                              ! (Wm-2)     usually equal to:
                                              !      emis_snow * LW_down
                                              !
REAL, DIMENSION(:), INTENT(IN)    :: PLW2     ! LW coef dependant   of TSNOW
                                              ! (Wm-2 K-4) usually equal to:
                                              ! -1 * emis_snow * stefan_constant
                                              !
REAL, DIMENSION(:), INTENT(IN)    :: PTA      ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA      ! specific humidity
                                                ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD    ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PPS      ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA    ! air density
                                                ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PSR      ! snow rate
REAL, DIMENSION(:), INTENT(IN)    :: PZREF    ! reference height of the first
                                              ! atmospheric level (temperature)
REAL, DIMENSION(:), INTENT(IN)    :: PUREF    ! reference height of the first
                                              ! atmospheric level (wind)
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSNOW  ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSNOW   ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESNOW  ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSNOW   ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT    ! snow melting rate (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_SNOW! heat storage inside snow
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW ! absorbed LW rad by snow (W/m2)
!
!
!*      0.2    declarations of local variables
!
REAL :: ZEXPL = 0.
REAL :: ZIMPL = 1.
!
REAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: ZEXNS, ZEXNA, ZDIRCOSZW
REAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: ZZ0      ! roughness length for momentum
REAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: ZZ0H     ! roughness length forheat
!
REAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: ZRI      ! Richardson number
REAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: ZAC      ! aerodynamical conductance
REAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: ZRA      ! aerodynamical resistance
REAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: ZCH      ! drag coefficient for heat
REAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: ZB, ZY   ! coefficients in Ts eq.
REAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: ZWSNOW   ! snow before evolution
REAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: ZSNOW_HC ! snow heat capacity
REAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: ZSNOW_TC ! snow thermal conductivity
REAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: ZSNOW_D  ! snow depth
REAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: ZMELT    ! snow melting rate (kg/m3/s)
REAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: ZTS_SNOW ! snow surface temperature
                                          ! at previous time-step
REAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: ZQSAT    ! specific humidity
!                                         ! for ice
REAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: ZDQSAT   ! d(specific humidity)/dT
!                                         ! for ice
!
REAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: ZSR1, ZSR2   ! norm. snow precip.
!
LOGICAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: GSNOWMASK ! where snow is
!                                             ! at previuos time-step
LOGICAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: GFLUXMASK ! where fluxes can
!                                             ! be computed at
!                                             ! new time-step
!                                             ! i.e. snow occurence
!                                             ! at previous time-step
!                                             ! OR snow fall
INTEGER, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: JSNOWMASK1, JSNOWMASK2, JSNOWMASK3 ! where snow is or not
!                                                                      ! at previuos time-step
INTEGER, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: JFLUXMASK ! where fluxes can
!                                             ! be computed at
!                                             ! new time-step
!                                             ! i.e. snow occurence
!                                             ! at previous time-step
!                                             ! OR snow fall
!
REAL :: ZWSNOW_MIN = 0.1 ! minimum value of snow content (kg/m2) for prognostic
!                        ! computations
!
REAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: ZEI_SNOW  ! internal energy of snow
REAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: ZPEI_SNOW ! internal energy of snow at t+
REAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: ZWORK1
REAL, DIMENSION(SIZE(TPSNOW%WSNOW,1)) :: ZDQSATI, ZQSATI
!
INTEGER                         :: JJ, JI, JCOMPT_SNOW1, JCOMPT_SNOW2, JCOMPT_SNOW3, JCOMPT_FLUX
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('SNOW_COVER_1LAYER',0,ZHOOK_HANDLE)
ZB(:)=0.
ZY(:)=0.
ZMELT  (:) = 0.
PMELT  (:) = 0.
PRNSNOW(:) = 0.
PHSNOW (:) = 0.
PLESNOW(:) = 0.
PGSNOW (:) = 0.
!RJ: workaround to prevent decomposition unstable xundef masks for Tx_LWA_SN_RD fields
!RJ: in TEB_DIAGNOSTICS.nc during TEB_GARDEN_GREENROOF_BEM_3L_IRRIG_* tests
!RJ: problem with decomposition handling somewhere else
#ifdef RJ_PFIX
PABS_LW(:) = 0.0
#else
PABS_LW(:) = XUNDEF
#endif
!
!* snow reservoir before evolution
!
ZWSNOW(:) = TPSNOW%WSNOW(:,1)
ZTS_SNOW(:) = MIN(XTT,PTG(:))
!
ZSNOW_D (:) = 0.
ZSNOW_TC(:) = 0.
ZSNOW_HC(:) = 0.
!
!-------------------------------------------------------------------------------
!
!*      1.1    most useful masks
!              -----------------
!
GSNOWMASK(:)=.FALSE.
GFLUXMASK(:)=.FALSE.
JSNOWMASK1(:)=0.
JSNOWMASK2(:)=0.
JSNOWMASK3(:)=0.
JFLUXMASK(:)=0.

  !*      1.2    drag
!              ----
!
!*      1.2.1  defaults
!              --------
!
!* variation of temperature with altitude is neglected
!
ZEXNS(:) = 1.
ZEXNA(:) = 1.
!
!* slope is neglected in drag computation
!
ZDIRCOSZW(:) = 1.
!
!* roughness length are imposed:
!
ZZ0   (:) = PZ0SN
ZZ0H  (:) = PZ0HSN

!
!
!*      1.1    most useful masks
!              -----------------
!* snow occurence at previous time-step
!
!* snow occurence during the time-step for fluxes computation
!
JCOMPT_SNOW1=0
JCOMPT_SNOW2=0
JCOMPT_SNOW3=0
JCOMPT_FLUX=0
DO JJ=1,SIZE(ZWSNOW)
  IF (ZWSNOW(JJ)>0.) THEN
    GSNOWMASK(JJ)=.TRUE.
    !* surface temperature
    ZTS_SNOW(JJ) = TPSNOW%TS(JJ)
    GFLUXMASK(JJ)=.TRUE.
    !gsnowmask=t
    JCOMPT_SNOW1=JCOMPT_SNOW1+1
    JSNOWMASK1(JCOMPT_SNOW1) = JJ
    !gfluxmask=t
    JCOMPT_FLUX=JCOMPT_FLUX+1
    JFLUXMASK(JCOMPT_FLUX) = JJ
    IF (ZWSNOW(JJ)>=ZWSNOW_MIN) THEN
      !second snow mask
      JCOMPT_SNOW3=JCOMPT_SNOW3+1
      JSNOWMASK3(JCOMPT_SNOW3)=JJ
    ELSE
      !lower limit of snow cover for prognostic computations
      !0.1 kg/m2 of snow water content
      TPSNOW%T(JJ,1) = MIN(PTG(JJ),XTT)
    ENDIF
  ELSE
    TPSNOW%T(JJ,1) = MIN(PTG(JJ),XTT)
    !gsnowmask=false
    JCOMPT_SNOW2=JCOMPT_SNOW2+1
    JSNOWMASK2(JCOMPT_SNOW2) = JJ
    IF (PSR(JJ)>0.) THEN
      GFLUXMASK(JJ)=.TRUE.
      JCOMPT_FLUX=JCOMPT_FLUX+1
      JFLUXMASK(JCOMPT_FLUX) = JJ
    ENDIF
  ENDIF
ENDDO
!
!-------------------------------------------------------------------------------
!
!*      1.2    drag
!              ----
!
!*      1.2.2   qsat (Tsnow)
!              ------------
!
ZQSAT(:) = QSATI(ZTS_SNOW(:), PPS(:) )
!
!*      1.2.3  Richardson number
!              -----------------
!
!* snow is present on all the considered surface.
!* computation occurs where snow is and/or falls.
!
 CALL SURFACE_RI(ZTS_SNOW, ZQSAT, ZEXNS, ZEXNA, PTA, PQA, &
                  PZREF, PUREF, ZDIRCOSZW, PVMOD, ZRI      )  
!
!*      1.2.4  Aerodynamical conductance
!              -------------------------
!
 CALL SURFACE_AERO_COND(ZRI, PZREF, PUREF, PVMOD, ZZ0, ZZ0H, ZAC, ZRA, ZCH)
!
!-------------------------------------------------------------------------------
!
!*      2.     snow thermal characteristics
!              ----------------------------
!cdir nodep
DO JJ=1,JCOMPT_SNOW1
  !
  JI = JSNOWMASK1(JJ)
  !
  !*      2.1    snow heat capacity
  ZSNOW_HC(JI) = TPSNOW%RHO(JI,1) * XCI * XRHOLI / XRHOLW
!*      2.2    snow depth
  ZSNOW_D(JI) = ZWSNOW(JI) / TPSNOW%RHO(JI,1)
!*      2.3    snow thermal conductivity
  ZSNOW_TC(JI) = XCONDI * (TPSNOW%RHO(JI,1)/XRHOLW)**1.885
!*      2.4    internal energy of snow
  ZEI_SNOW(JI) = ZSNOW_HC(JI)*ZSNOW_D(JI)*TPSNOW%T(JI,1)
  !
ENDDO
!
!cdir nodep
DO JJ=1,JCOMPT_SNOW2
  !
  JI = JSNOWMASK2(JJ)
  !
  !*      2.1    snow heat capacity
  ZSNOW_HC(JI) = PRHOSMIN * XCI * XRHOLI / XRHOLW
!*      2.2    snow depth
  ZSNOW_D(JI) = PTSTEP * PSR(JI) / PRHOSMIN
!*      2.3    snow thermal conductivity
  ZSNOW_TC(JI) = XCONDI * (PRHOSMIN /XRHOLW)**1.885
!*      2.4    internal energy of snow
  ZEI_SNOW(JI) = 0.
!
ENDDO
!
!-------------------------------------------------------------------------------
!
!*      3.     Snow temperature evolution
!              --------------------------
!
!*      3.5    dqsat/ dT (Tsnow)
!              -----------------
!
ZDQSATI = DQSATI(ZTS_SNOW(:),PPS(:),ZQSAT(:))
WHERE (GSNOWMASK(:) .AND. ZWSNOW(:)>=ZWSNOW_MIN)
  ZDQSAT(:) = ZDQSATI(:)
END WHERE
!
!*      3.1    coefficients from Temperature tendency
!              --------------------------------------
!
!cdir nodep
DO JJ=1,JCOMPT_SNOW3
!
  JI=JSNOWMASK3(JJ)

  ZWORK1(JI) = ZSNOW_D(JI) * ZSNOW_HC(JI) / PTSTEP
!
  ZB(JI) = ZB(JI) + ZWORK1(JI)
!
!*      3.2    coefficients from solar radiation
!          ---------------------------------
!    
  ZY(JI) = ZY(JI) + ZWORK1(JI) * TPSNOW%T(JI,1) + PABS_SW(JI)
!
!
!*      3.3    coefficients from infra-red radiation
!              -------------------------------------
!
  ZWORK1(JI) = PLW2(JI) * TPSNOW%T(JI,1)**3
!
  ZB(JI) = ZB(JI) - 4 * ZIMPL * ZWORK1(JI)
!
  ZY(JI) = ZY(JI) + PLW1(JI) + ZWORK1(JI) * (ZEXPL-3.*ZIMPL) * TPSNOW%T(JI,1)
!
!
!*      3.4    coefficients from sensible heat flux
!              ------------------------------------
!
 ZWORK1(JI) = XCPD * PRHOA(JI) * ZAC(JI)
! 
  ZB(JI) = ZB(JI) + ZWORK1(JI) *   ZIMPL
!
  ZY(JI) = ZY(JI) - ZWORK1(JI) * ( ZEXPL * TPSNOW%T(JI,1) - PTA(JI) )
!
!
!*      3.6    coefficients from latent heat flux
!              ----------------------------------
!
  ZWORK1(JI) =  XLSTT * PRHOA(JI) * ZAC(JI)
!
  ZB(JI) = ZB(JI) + ZWORK1(JI) *  ZIMPL * ZDQSAT(JI)
!
  ZY(JI) = ZY(JI) - ZWORK1(JI) * (  ZQSAT(JI) - PQA(JI) - ZIMPL * ZDQSAT(JI)*TPSNOW%T(JI,1) )
!
!*      3.7    coefficients from conduction flux at snow base
!              ----------------------------------------------
!
  ZWORK1(JI) = ZSNOW_TC(JI)/(0.5*ZSNOW_D(JI))
!
  ZB(JI) = ZB(JI) + ZWORK1(JI) *  ZIMPL / ( 1. + ZWORK1(JI)*PTG_COEFA(JI) )
!
  ZY(JI) = ZY(JI) - ZWORK1(JI) * (ZEXPL * TPSNOW%T(JI,1) - PTG_COEFB(JI)) &
                   / ( 1. + ZWORK1(JI)*PTG_COEFA(JI) )
!
!*      3.8    guess of snow temperature before accumulation and melting
!              ---------------------------------------------------------
!
  TPSNOW%T(JI,1) = ZY(JI) / ZB(JI)
!
ENDDO
!
!-------------------------------------------------------------------------------
!
!*      4.     Snow melt
!              ---------
!
!*      4.1    melting
!              -------
!
!cdir nodep
DO JJ=1,JCOMPT_SNOW1
!
  JI = JSNOWMASK1(JJ)
!
  ZMELT(JI)  = MAX( TPSNOW%T(JI,1) - XTT , 0. ) * ZSNOW_HC(JI) /  XLMTT / PTSTEP
!
  ZMELT(JI)  = MIN( ZMELT(JI) , ZWSNOW(JI) / ZSNOW_D(JI) / PTSTEP )
!
  TPSNOW%T(JI,1) = MIN( TPSNOW%T(JI,1) , XTT )
!
ENDDO
!
!*      4.2    run-off of all snow if lower surface temperature is positive
!              ------------------------------------------------------------
!
!* this option is used when snow is located on sloping roofs for example.
!
IF (OALL_MELT) THEN
  WHERE ( GSNOWMASK(:) .AND. PTG(:)>XTT .AND. ZWSNOW(:)>=ZWSNOW_MIN )
    PMELT(:) = PMELT(:) + ZWSNOW(:) / PTSTEP
  END WHERE
END IF
!
!*      4.3    output melting in kg/m2/s
!              -------------------------
!
PMELT(:) = ZMELT(:) * ZSNOW_D(:)
!
!-------------------------------------------------------------------------------
!
!*      5.     fluxes
!              ------
!
!*      5.3    qsat (Tsnow)
!              ------------
!
ZQSATI = QSATI(TPSNOW%T(:,1),PPS(:))
WHERE (GFLUXMASK(:)) 
   ZQSAT(:) = ZQSATI(:)
END WHERE
!
!*      5.1    net radiation (with Ts lin. extrapolation)
!              -------------
!
!cdir nodep
DO JJ = 1, JCOMPT_FLUX
!
  JI = JFLUXMASK(JJ)
!
  PABS_LW(JI) =  PLW1(JI) + PLW2(JI) * TPSNOW%T(JI,1)**4
!
  PRNSNOW(JI) = PABS_SW(JI) + PABS_LW(JI)
!
!
!*      5.2    sensible heat flux
!              ------------------
!
  PHSNOW(JI) = XCPD * PRHOA(JI) * ZAC(JI) * ( TPSNOW%T(JI,1) - PTA(JI) )
!
!
!*      5.4    latent heat flux
!              ----------------
!
  PLESNOW(JI) = XLSTT * PRHOA(JI) * ZAC(JI) * ( ZQSAT(JI) - PQA(JI) )
  !
!
!*      5.5    Conduction heat flux
!              --------------------
!
  !PGSNOW(JI) = ZSNOW_TC(JI)/(0.5*ZSNOW_D(JI)) * ( TPSNOW%T(JI,1) - PTG(JI) )
  PGSNOW(JI) = ZSNOW_TC(JI)/(0.5*ZSNOW_D(JI)) * ( TPSNOW%T(JI,1) - PTG_COEFB(JI) ) &
             / ( 1. + ZSNOW_TC(JI)/(0.5*ZSNOW_D(JI))*PTG_COEFA(JI) )
!
!
!*      5.6    If ground T>0 C, Melting is estimated from conduction heat flux
!              ---------------------------------------------------------------
!
  IF (PTG(JI)>XTT)  PMELT(JI) = MAX (PMELT(JI), -PGSNOW(JI)/XLMTT)
!
ENDDO
!
!-------------------------------------------------------------------------------
!
!*      6.     reservoir evolution
!              -------------------
!
!cdir nodep
DO JJ = 1, SIZE(TPSNOW%WSNOW,1)
!
!*      6.1    snow fall
!              ---------
!
  TPSNOW%WSNOW(JJ,1) = TPSNOW%WSNOW(JJ,1) + PTSTEP * PSR(JJ)
!
!
!*      6.2    sublimation
!              -----------
!
  PLESNOW(JJ) = MIN( PLESNOW(JJ), XLSTT*TPSNOW%WSNOW(JJ,1)/PTSTEP )
!
  TPSNOW%WSNOW(JJ,1)  = MAX( TPSNOW%WSNOW(JJ,1) - PTSTEP * PLESNOW(JJ)/XLSTT , 0.)
!
  IF ( TPSNOW%WSNOW(JJ,1)<1.E-8 * PTSTEP ) TPSNOW%WSNOW(JJ,1) = 0.
!
!*      6.3    melting
!              -------
!
  PMELT(JJ) = MIN( PMELT(JJ), TPSNOW%WSNOW(JJ,1)/PTSTEP )
!
  TPSNOW%WSNOW(JJ,1)= MAX( TPSNOW%WSNOW(JJ,1) - PTSTEP * PMELT(JJ) , 0.)
!
  IF ( TPSNOW%WSNOW(JJ,1)<1.E-8 * PTSTEP ) TPSNOW%WSNOW(JJ,1) = 0.
!
  IF (TPSNOW%WSNOW(JJ,1)==0.) PGSNOW(JJ) = MAX ( PGSNOW(JJ), - PMELT(JJ)*XLMTT )
!
ENDDO
!
!*      6.4    time dependent drainage
!              -----------------------
!
IF (PDRAIN_TIME>0.) THEN
  WHERE ( TPSNOW%WSNOW(:,1)>0.)
    TPSNOW%WSNOW(:,1) = TPSNOW%WSNOW(:,1) * EXP(-PTSTEP/PDRAIN_TIME/XDAY)
  END WHERE
END IF
!
!*      6.5    melting of last 1mm of snow depth
!              ---------------------------------
!
WHERE ( TPSNOW%WSNOW(:,1)<ZWSNOW_MIN .AND. PMELT(:)>0. .AND. PSR(:)==0. )
  PMELT(:) = PMELT(:) + TPSNOW%WSNOW(:,1) / PTSTEP
  TPSNOW%WSNOW(:,1)=0.
END WHERE
!
WHERE ( TPSNOW%WSNOW(:,1)<1.E-8 * PTSTEP ) 
   TPSNOW%WSNOW(:,1) = 0.
END WHERE
!
!-------------------------------------------------------------------------------
!
!*      7.     albedo evolution
!              ----------------
!
!*      7.1    If melting occurs or not
!              -----------------------
!
!
!cdir nodep
DO JJ=1,JCOMPT_SNOW1
!
  JI = JSNOWMASK1(JJ)
!
  IF (PMELT(JI) > 0. ) THEN
!
    TPSNOW%ALB(JI) = (TPSNOW%ALB(JI)-PANSMIN)*EXP(-PRHOFOLD*PTSTEP/XDAY) + PANSMIN   &
                     + PSR(JI)*PTSTEP/PWCRN*PANSMAX  
!
  ELSEIF (PMELT(JI)==0.) THEN 
    TPSNOW%ALB(JI) = TPSNOW%ALB(JI) - PTODRY*PTSTEP/XDAY  + PSR(JI)*PTSTEP/PWCRN*PANSMAX  
!
  ENDIF
!
ENDDO
!
!-------------------------------------------------------------------------------
!
!*      8.     density evolution
!              -----------------
!
!*      8.1    old snow
!              --------
!
!cdir nodep
DO JJ = 1, JCOMPT_SNOW1
! 
  JI = JSNOWMASK1(JJ)
!
  IF (TPSNOW%WSNOW(JI,1)>0. ) THEN
!
    ZSR1(JI) = MAX( TPSNOW%WSNOW(JI,1) , PSR(JI) * PTSTEP )
!
    TPSNOW%RHO(JI,1) = (TPSNOW%RHO(JI,1)-PRHOSMAX)*EXP(-PRHOFOLD*PTSTEP/XDAY) + PRHOSMAX
    TPSNOW%RHO(JI,1) = ( (ZSR1(JI)-PSR(JI)*PTSTEP) * TPSNOW%RHO(JI,1)    &
                  + (PSR(JI)*PTSTEP) * PRHOSMIN ) / ZSR1(JI)  
  ENDIF
!
ENDDO
!
!*      8.2    fresh snow
!              ----------
!
!cdir nodep
DO JJ=1,SIZE(TPSNOW%WSNOW,1)
  IF (  TPSNOW%WSNOW(JJ,1)>0. ) THEN
    TPSNOW%ALB(JJ) = MAX(TPSNOW%ALB(JJ),PANSMIN)
    TPSNOW%ALB(JJ) = MIN(TPSNOW%ALB(JJ),PANSMAX)
    IF (ZWSNOW(JJ)==0.) THEN
      TPSNOW%ALB  (JJ) = PANSMAX
      TPSNOW%EMIS (JJ) = XEMISSN
      TPSNOW%RHO(JJ,1) = PRHOSMIN
    ENDIF
  ENDIF
    ENDDO
!
!-------------------------------------------------------------------------------
!
!*      9.     fresh snow accumulation (if more than 1mm of snow depth)
!              -----------------------
!
!cdir nodep
DO JJ=1,JCOMPT_SNOW3
!
  JI = JSNOWMASK3(JJ)
!
  IF (PSR(JI)>0. .AND. TPSNOW%WSNOW(JI,1)>0.) THEN
!
    ZSR2(JI) = MIN( TPSNOW%WSNOW(JI,1) , PSR(JI) * PTSTEP )
!
    TPSNOW%T(JI,1) =( ( TPSNOW%WSNOW(JI,1) - ZSR2(JI) ) *  TPSNOW%T(JI,1)        &
              +   ZSR2(JI)  * MIN( PTA   (JI) ,XTT )) / ( TPSNOW%WSNOW(JI,1) )  
  ENDIF
!
ENDDO
!
!-------------------------------------------------------------------------------
!
!*     10.     Surface temperature
!              -------------------
!
!* note that if the relation between snow pack temperature and its
!  surface temperature is modified, think to modify it also in
!  subroutine init_snow_lw.f90
!
WHERE (GSNOWMASK(:) )
  TPSNOW%TS(:) = TPSNOW%T(:,1)
END WHERE
!
!-------------------------------------------------------------------------------
!
!*     11.     bogus values
!              ------------
!
!*     11.1    snow characteristics where snow IS present at current time-step
!              ---------------------------------------------------------------
!
WHERE (TPSNOW%WSNOW(:,1)==0.)
  TPSNOW%T(:,1)   = XUNDEF
  TPSNOW%RHO(:,1) = XUNDEF
  TPSNOW%ALB(:)   = XUNDEF
  TPSNOW%TS(:)    = XUNDEF
  TPSNOW%EMIS(:)  = XUNDEF
END WHERE
!
!
!-------------------------------------------------------------------------------
!
!*     12.     Heat storage inside snow pack
!
WHERE (GSNOWMASK(:))
  ZPEI_SNOW(:) = ZSNOW_HC(:)*ZSNOW_D(:)*TPSNOW%T(:,1)
ELSEWHERE
  ZPEI_SNOW(:) = 0.
END WHERE
PDQS_SNOW(:) = (ZPEI_SNOW(:)-ZEI_SNOW(:))/PTSTEP
!
IF (LHOOK) CALL DR_HOOK('SNOW_COVER_1LAYER',1,ZHOOK_HANDLE)

!------------------------------------------------------------------------------- !
END SUBROUTINE SNOW_COVER_1LAYER

