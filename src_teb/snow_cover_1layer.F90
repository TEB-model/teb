!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE SNOW_COVER_1LAYER(PTSTEP, PANSMIN, PANSMAX, PTODRY,         &
                                   PRHOSMIN, PRHOSMAX, PRHOFOLD, OALL_MELT,  &
                                   PDRAIN_TIME, PWCRN, PZ0SN, PZ0HSN,        &
                                   PTSNOW, PASNOW, PRSNOW, PWSNOW, PTS_SNOW, &
                                   PESNOW,                                   &
                                   PTG, PTG_COEFA, PTG_COEFB,                &
                                   PABS_SW, PLW1, PLW2,                      &
                                   PTA, PQA, PVMOD, PPS, PRHOA, PSR,         &
                                   PZREF, PUREF,                             &
                                   PRNSNOW, PHSNOW, PLESNOW, PGSNOW, PMELT,  &
                                   PDQS_SNOW, PABS_LW                        )  
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
REAL, DIMENSION(:), INTENT(INOUT) :: PWSNOW   ! snow reservoir (kg/m2)
REAL, DIMENSION(:), INTENT(INOUT) :: PTSNOW   ! snow temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PASNOW   ! snow albedo
REAL, DIMENSION(:), INTENT(INOUT) :: PRSNOW   ! snow density
REAL, DIMENSION(:), INTENT(INOUT) :: PTS_SNOW ! snow surface temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PESNOW   ! snow emissivity
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
REAL, DIMENSION(SIZE(PWSNOW)) :: ZEXNS, ZEXNA, ZDIRCOSZW
REAL, DIMENSION(SIZE(PWSNOW)) :: ZZ0      ! roughness length for momentum
REAL, DIMENSION(SIZE(PWSNOW)) :: ZZ0H     ! roughness length forheat
!
REAL, DIMENSION(SIZE(PWSNOW)) :: ZRI      ! Richardson number
REAL, DIMENSION(SIZE(PWSNOW)) :: ZAC      ! aerodynamical conductance
REAL, DIMENSION(SIZE(PWSNOW)) :: ZRA      ! aerodynamical resistance
REAL, DIMENSION(SIZE(PWSNOW)) :: ZCH      ! drag coefficient for heat
REAL, DIMENSION(SIZE(PWSNOW)) :: ZB, ZY   ! coefficients in Ts eq.
REAL, DIMENSION(SIZE(PWSNOW)) :: ZWSNOW   ! snow before evolution
REAL, DIMENSION(SIZE(PWSNOW)) :: ZSNOW_HC ! snow heat capacity
REAL, DIMENSION(SIZE(PWSNOW)) :: ZSNOW_TC ! snow thermal conductivity
REAL, DIMENSION(SIZE(PWSNOW)) :: ZSNOW_D  ! snow depth
REAL, DIMENSION(SIZE(PWSNOW)) :: ZMELT    ! snow melting rate (kg/m3/s)
REAL, DIMENSION(SIZE(PWSNOW)) :: ZTS_SNOW ! snow surface temperature
                                          ! at previous time-step
REAL, DIMENSION(SIZE(PWSNOW)) :: ZQSAT    ! specific humidity
!                                         ! for ice
REAL, DIMENSION(SIZE(PWSNOW)) :: ZDQSAT   ! d(specific humidity)/dT
!                                         ! for ice
!
REAL, DIMENSION(SIZE(PWSNOW)) :: ZSR1, ZSR2   ! norm. snow precip.
!
LOGICAL, DIMENSION(SIZE(PWSNOW)) :: GSNOWMASK ! where snow is
!                                             ! at previuos time-step
LOGICAL, DIMENSION(SIZE(PWSNOW)) :: GFLUXMASK ! where fluxes can
!                                             ! be computed at
!                                             ! new time-step
!                                             ! i.e. snow occurence
!                                             ! at previous time-step
!                                             ! OR snow fall
INTEGER, DIMENSION(SIZE(PWSNOW)) :: JSNOWMASK1, JSNOWMASK2, JSNOWMASK3 ! where snow is or not
!                                                                      ! at previuos time-step
INTEGER, DIMENSION(SIZE(PWSNOW)) :: JFLUXMASK ! where fluxes can
!                                             ! be computed at
!                                             ! new time-step
!                                             ! i.e. snow occurence
!                                             ! at previous time-step
!                                             ! OR snow fall
!
REAL :: ZWSNOW_MIN = 0.1 ! minimum value of snow content (kg/m2) for prognostic
!                        ! computations
!
REAL, DIMENSION(SIZE(PWSNOW)) :: ZEI_SNOW  ! internal energy of snow
REAL, DIMENSION(SIZE(PWSNOW)) :: ZPEI_SNOW ! internal energy of snow at t+
REAL, DIMENSION(SIZE(PWSNOW)) :: ZWORK1
REAL, DIMENSION(SIZE(PWSNOW)) :: ZDQSATI, ZQSATI
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
ZWSNOW(:) = PWSNOW(:)
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
    ZTS_SNOW(JJ)=PTS_SNOW(JJ)
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
      PTSNOW(JJ)=MIN(PTG(JJ),XTT)
    ENDIF
  ELSE
    PTSNOW(JJ)=MIN(PTG(JJ),XTT)
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
  ZSNOW_HC(JI) = PRSNOW(JI) * XCI * XRHOLI / XRHOLW
!*      2.2    snow depth
  ZSNOW_D(JI) = ZWSNOW(JI) / PRSNOW(JI)
!*      2.3    snow thermal conductivity
  ZSNOW_TC(JI) = XCONDI * (PRSNOW(JI)/XRHOLW)**1.885
!*      2.4    internal energy of snow
  ZEI_SNOW(JI) = ZSNOW_HC(JI)*ZSNOW_D(JI)*PTSNOW(JI)
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
  ZY(JI) = ZY(JI) + ZWORK1(JI) * PTSNOW(JI) + PABS_SW(JI)
!
!
!*      3.3    coefficients from infra-red radiation
!              -------------------------------------
!
  ZWORK1(JI) = PLW2(JI) * PTSNOW(JI)**3
!
  ZB(JI) = ZB(JI) - 4 * ZIMPL * ZWORK1(JI)
!
  ZY(JI) = ZY(JI) + PLW1(JI) + ZWORK1(JI) * (ZEXPL-3.*ZIMPL) * PTSNOW(JI)
!
!
!*      3.4    coefficients from sensible heat flux
!              ------------------------------------
!
 ZWORK1(JI) = XCPD * PRHOA(JI) * ZAC(JI)
! 
  ZB(JI) = ZB(JI) + ZWORK1(JI) *   ZIMPL
!
  ZY(JI) = ZY(JI) - ZWORK1(JI) * ( ZEXPL * PTSNOW(JI) - PTA(JI) )
!
!
!*      3.6    coefficients from latent heat flux
!              ----------------------------------
!
  ZWORK1(JI) =  XLSTT * PRHOA(JI) * ZAC(JI)
!
  ZB(JI) = ZB(JI) + ZWORK1(JI) *  ZIMPL * ZDQSAT(JI)
!
  ZY(JI) = ZY(JI) - ZWORK1(JI) * (  ZQSAT(JI) - PQA(JI) - ZIMPL * ZDQSAT(JI)*PTSNOW(JI) )
!
!*      3.7    coefficients from conduction flux at snow base
!              ----------------------------------------------
!
  ZWORK1(JI) = ZSNOW_TC(JI)/(0.5*ZSNOW_D(JI))
!
  ZB(JI) = ZB(JI) + ZWORK1(JI) *  ZIMPL / ( 1. + ZWORK1(JI)*PTG_COEFA(JI) )
!
  ZY(JI) = ZY(JI) - ZWORK1(JI) * (ZEXPL * PTSNOW(JI) - PTG_COEFB(JI)) &
                   / ( 1. + ZWORK1(JI)*PTG_COEFA(JI) )
!
!*      3.8    guess of snow temperature before accumulation and melting
!              ---------------------------------------------------------
!
  PTSNOW(JI) = ZY(JI) / ZB(JI)
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
  ZMELT(JI)  = MAX( PTSNOW(JI) - XTT , 0. ) * ZSNOW_HC(JI) /  XLMTT / PTSTEP
!
  ZMELT(JI)  = MIN( ZMELT(JI) , ZWSNOW(JI) / ZSNOW_D(JI) / PTSTEP )
!
  PTSNOW(JI) = MIN( PTSNOW(JI) , XTT )
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
ZQSATI = QSATI(PTSNOW(:),PPS(:))
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
  PABS_LW(JI) =  PLW1(JI) + PLW2(JI) * PTSNOW(JI)**4
!
  PRNSNOW(JI) = PABS_SW(JI) + PABS_LW(JI)
!
!
!*      5.2    sensible heat flux
!              ------------------
!
  PHSNOW(JI) = XCPD * PRHOA(JI) * ZAC(JI) * ( PTSNOW(JI) - PTA(JI) )
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
  !PGSNOW(JI) = ZSNOW_TC(JI)/(0.5*ZSNOW_D(JI)) * ( PTSNOW(JI) - PTG(JI) )
  PGSNOW(JI) = ZSNOW_TC(JI)/(0.5*ZSNOW_D(JI)) * ( PTSNOW(JI) - PTG_COEFB(JI) ) &
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
DO JJ = 1, SIZE(PWSNOW)
!
!*      6.1    snow fall
!              ---------
!
  PWSNOW(JJ) = PWSNOW(JJ) + PTSTEP * PSR(JJ)
!
!
!*      6.2    sublimation
!              -----------
!
  PLESNOW(JJ) = MIN( PLESNOW(JJ), XLSTT*PWSNOW(JJ)/PTSTEP )
!
  PWSNOW(JJ)  = MAX( PWSNOW(JJ) - PTSTEP * PLESNOW(JJ)/XLSTT , 0.)
!
  IF ( PWSNOW(JJ)<1.E-8 * PTSTEP ) PWSNOW(JJ) = 0.
!
!*      6.3    melting
!              -------
!
  PMELT(JJ) = MIN( PMELT(JJ), PWSNOW(JJ)/PTSTEP )
!
  PWSNOW(JJ)= MAX( PWSNOW(JJ) - PTSTEP * PMELT(JJ) , 0.)
!
  IF ( PWSNOW(JJ)<1.E-8 * PTSTEP ) PWSNOW(JJ) = 0.
!
  IF (PWSNOW(JJ)==0.) PGSNOW(JJ) = MAX ( PGSNOW(JJ), - PMELT(JJ)*XLMTT )
!
ENDDO
!
!*      6.4    time dependent drainage
!              -----------------------
!
IF (PDRAIN_TIME>0.) THEN
  WHERE ( PWSNOW(:)>0.)
    PWSNOW(:) = PWSNOW(:) * EXP(-PTSTEP/PDRAIN_TIME/XDAY)
  END WHERE
END IF
!
!*      6.5    melting of last 1mm of snow depth
!              ---------------------------------
!
WHERE ( PWSNOW(:)<ZWSNOW_MIN .AND. PMELT(:)>0. .AND. PSR(:)==0. )
  PMELT(:) = PMELT(:) + PWSNOW(:) / PTSTEP
  PWSNOW(:)=0.
END WHERE
!
WHERE ( PWSNOW(:)<1.E-8 * PTSTEP ) 
   PWSNOW(:) = 0.
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
    PASNOW(JI) = (PASNOW(JI)-PANSMIN)*EXP(-PRHOFOLD*PTSTEP/XDAY) + PANSMIN   &
                + PSR(JI)*PTSTEP/PWCRN*PANSMAX  
!
  ELSEIF (PMELT(JI)==0.) THEN 
    PASNOW(JI) = PASNOW(JI) - PTODRY*PTSTEP/XDAY                          &
                + PSR(JI)*PTSTEP/PWCRN*PANSMAX  
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
  IF (PWSNOW(JI)>0. ) THEN
!
    ZSR1(JI) = MAX( PWSNOW(JI) , PSR(JI) * PTSTEP )
!
    PRSNOW(JI) = (PRSNOW(JI)-PRHOSMAX)*EXP(-PRHOFOLD*PTSTEP/XDAY) + PRHOSMAX
    PRSNOW(JI) = ( (ZSR1(JI)-PSR(JI)*PTSTEP) * PRSNOW(JI)    &
                  + (PSR(JI)*PTSTEP) * PRHOSMIN ) / ZSR1(JI)  
  ENDIF
!
ENDDO
!
!*      8.2    fresh snow
!              ----------
!
!cdir nodep
DO JJ=1,SIZE(PWSNOW)
  IF (  PWSNOW(JJ)>0. ) THEN
    PASNOW(JJ) = MAX(PASNOW(JJ),PANSMIN)
    PASNOW(JJ) = MIN(PASNOW(JJ),PANSMAX)
    IF (ZWSNOW(JJ)==0.) THEN
      PASNOW(JJ) = PANSMAX
      PESNOW(JJ) = XEMISSN
      PRSNOW(JJ) = PRHOSMIN
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
  IF (PSR(JI)>0. .AND. PWSNOW(JI)>0.) THEN
!
    ZSR2(JI) = MIN( PWSNOW(JI) , PSR(JI) * PTSTEP )
!
    PTSNOW(JI) =( ( PWSNOW(JI) - ZSR2(JI) ) *      PTSNOW(JI)        &
              +                ZSR2(JI)   * MIN( PTA   (JI) ,XTT ))&
                /(   PWSNOW(JI) )  
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
  PTS_SNOW(:) = PTSNOW(:)
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
WHERE (PWSNOW(:)==0.)
  PTSNOW  (:) = XUNDEF
  PRSNOW  (:) = XUNDEF
  PASNOW  (:) = XUNDEF
  PTS_SNOW(:) = XUNDEF
  PESNOW  (:) = XUNDEF
END WHERE
!
!
!-------------------------------------------------------------------------------
!
!*     12.     Heat storage inside snow pack
!
WHERE (GSNOWMASK(:))
  ZPEI_SNOW(:) = ZSNOW_HC(:)*ZSNOW_D(:)*PTSNOW(:)
ELSEWHERE
  ZPEI_SNOW(:) = 0.
END WHERE
PDQS_SNOW(:) = (ZPEI_SNOW(:)-ZEI_SNOW(:))/PTSTEP
!
IF (LHOOK) CALL DR_HOOK('SNOW_COVER_1LAYER',1,ZHOOK_HANDLE)

!------------------------------------------------------------------------------- !
END SUBROUTINE SNOW_COVER_1LAYER

