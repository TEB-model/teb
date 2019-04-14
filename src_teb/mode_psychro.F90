!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!##################
MODULE MODE_PSYCHRO
!##################
!
!!****  *MODE_PSYCHRO* -
!!
!!    PURPOSE
!!    -------
!      
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!       NONE          
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    12/04/11
!!      J.Escobar   11/13 :  remove space in ELSEWHERE statement
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
interface PE_FROM_PQ
        module procedure PE_FROM_PQ_0D
        module procedure PE_FROM_PQ_1D
end interface
interface TD_FROM_TQ
        module procedure TD_FROM_TQ_0D
        module procedure TD_FROM_TQ_1D
end interface
interface RV_FROM_TPTWB
        module procedure RV_FROM_TPTWB_0D
        module procedure RV_FROM_TPTWB_1D
end interface
interface TWB_FROM_TPQ
        module procedure TWB_FROM_TPQ_0D
        module procedure TWB_FROM_TPQ_1D
end interface
INTERFACE ENTH_FN_T_Q
  MODULE PROCEDURE ENTH_FN_T_Q
END INTERFACE
INTERFACE Q_FN_T_ENTH
  MODULE PROCEDURE Q_FN_T_ENTH
END INTERFACE

contains
!PE_FROM_PQ
!----------
function PE_FROM_PQ_0D(PP, PQ) RESULT(PE) 
!arguments and result
REAL, INTENT(IN) :: PP !atmos. pressure (Pa)
REAL, INTENT(IN) :: PQ !specific humidity (kg/kg)
REAL :: PE !water vapour pressure (Pa)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PSYCHRO:PE_FROM_PQ_0D',0,ZHOOK_HANDLE)
PE = PQ * PP /(0.622 + 0.378 * PQ)
IF (LHOOK) CALL DR_HOOK('MODE_PSYCHRO:PE_FROM_PQ_0D',1,ZHOOK_HANDLE)
end function PE_FROM_PQ_0D

function PE_FROM_PQ_1D(PP, PQ) RESULT(PE) 
!arguments and result
REAL, DIMENSION(:), INTENT(IN) :: PP !atmos. pressure (Pa)
REAL, DIMENSION(:), INTENT(IN) :: PQ !specific humidity (kg/kg)
REAL, DIMENSION(SIZE(PQ))      :: PE !water vapour pressure (Pa)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PSYCHRO:PE_FROM_PQ_1D',0,ZHOOK_HANDLE)
PE(:) = PQ(:) * PP(:) /(0.622 + 0.378 * PQ(:))
IF (LHOOK) CALL DR_HOOK('MODE_PSYCHRO:PE_FROM_PQ_1D',1,ZHOOK_HANDLE)
end function PE_FROM_PQ_1D
!-------------------------

!TD_FROM_TQ
function TD_FROM_TQ_0D(PT, PQ) RESULT(PTD)
USE MODD_CSTS
USE MODD_SURF_PAR, ONLY: XUNDEF
!arguments and result
REAL, INTENT(IN) :: PT !Air Temp. (K)
REAL, INTENT(IN) :: PQ !Specific humidity (kg/kg)
REAL :: PTD !Dew Point Air Temp. (K)
!local variables
REAL :: ALPHA
REAL :: ZPE !water vapour pressure
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('MODE_PSYCHRO:TD_FROM_TQ_0D',0,ZHOOK_HANDLE)
ZPE = PE_FROM_PQ(PT, PQ)
ALPHA = LOG(ZPE/1000.)
IF (PT .GE. XTT .AND. PT .GE. 93.+XTT) THEN
        PTD = XTT+6.54+14.526*ALPHA+0.7389*ALPHA*ALPHA+0.09486*ALPHA**3 &
              +0.4569*(ZPE/1000.)**0.1984
ELSE IF (PT .LT. XTT) THEN
        PTD = XTT+6.09+12.608*ALPHA+0.4959*ALPHA*ALPHA
ELSE
        PTD = XUNDEF
ENDIF
PTD = MIN(PTD, PT)
IF (LHOOK) CALL DR_HOOK('MODE_PSYCHRO:TD_FROM_TQ_0D',1,ZHOOK_HANDLE)
end function TD_FROM_TQ_0D

function TD_FROM_TQ_1D(PT, PQ) RESULT(PTD)
USE MODD_CSTS
USE MODD_SURF_PAR, ONLY: XUNDEF
!arguments and result
REAL, DIMENSION(:), INTENT(IN) :: PT !Air Temp. (K)
REAL, DIMENSION(:), INTENT(IN) :: PQ !Specific humidity (kg/kg)
REAL, DIMENSION(SIZE(PQ))      :: PTD !Dew Point Air Temp. (K)
!local variables
REAL, DIMENSION(SIZE(PQ)) :: ALPHA
REAL, DIMENSION(SIZE(PQ)) :: ZPE !water vapour pressure
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('MODE_PSYCHRO:TD_FROM_TQ_1D',0,ZHOOK_HANDLE)
ZPE = PE_FROM_PQ(PT, PQ)
ALPHA(:) = LOG(ZPE(:)/1000.)
WHERE (PT .GE. XTT .AND. PT .GE. 93.+XTT)
        PTD = XTT+6.54+14.526*ALPHA+0.7389*ALPHA*ALPHA+0.09486*ALPHA**3 &
              +0.4569*(ZPE/1000.)**0.1984
      ELSEWHERE (PT .LT. XTT)
        PTD = XTT+6.09+12.608*ALPHA+0.4959*ALPHA*ALPHA
ELSEWHERE
        PTD = XUNDEF
END WHERE
PTD(:) = MIN(PTD(:), PT(:))
IF (LHOOK) CALL DR_HOOK('MODE_PSYCHRO:TD_FROM_TQ_1D',1,ZHOOK_HANDLE)
end function TD_FROM_TQ_1D
!-------------------------

!RV_FROM_TPTWB
function RV_FROM_TPTWB_0D(PT, PP, PTWB) RESULT(PRV)
USE MODE_THERMOS
USE MODD_CSTS
!arguments and result
REAL, INTENT(IN) :: PT !Air temperature (K)
REAL, INTENT(IN) :: PP !Atmos. Pressure (Pa)
REAL, INTENT(IN) :: PTWB !Wet Bulb Temp. (K)
REAL :: PRV !water vapor mixing ratio (kg/kg)
REAL :: ZRVSAT !saturation water vapor mixing ratio
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('MODE_PSYCHRO:RV_FROM_TPTWB_0D',0,ZHOOK_HANDLE)
ZRVSAT = QSAT(PT, PP) / (1 - QSAT(PT, PP))
PRV = ((2501. - 2.326*(PTWB-XTT))*ZRVSAT - 1.006*(PT - PTWB)) &
       / (2501. + 1.86*(PT - XTT) -4.186*(PTWB - XTT))
IF (LHOOK) CALL DR_HOOK('MODE_PSYCHRO:RV_FROM_TPTWB_0D',1,ZHOOK_HANDLE)
end function RV_FROM_TPTWB_0D

function RV_FROM_TPTWB_1D(PT, PP, PTWB) RESULT(PRV)
USE MODE_THERMOS
USE MODD_CSTS
!arguments and result
REAL, DIMENSION(:), INTENT(IN) :: PT !Air temperature (K)
REAL, DIMENSION(:),INTENT(IN) :: PP !Atmos. Pressure (Pa)
REAL, DIMENSION(:),INTENT(IN) :: PTWB !Wet Bulb Temp. (K)
REAL, DIMENSION(SIZE(PT)) :: PRV !water vapor mixing ratio (kg/kg)
REAL, DIMENSION(SIZE(PT)) :: ZRVSAT !saturation water vapor mixing ratio
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('MODE_PSYCHRO:RV_FROM_TPTWB_1D',0,ZHOOK_HANDLE)
ZRVSAT = QSAT(PT, PP) / (1 - QSAT(PT, PP))
PRV(:) = ((2501. - 2.326*(PTWB(:)-XTT))*ZRVSAT(:) - 1.006*(PT(:) - PTWB(:))) &
       / (2501. + 1.86*(PT(:) - XTT) -4.186*(PTWB(:) - XTT))
IF (LHOOK) CALL DR_HOOK('MODE_PSYCHRO:RV_FROM_TPTWB_1D',1,ZHOOK_HANDLE)       
end function RV_FROM_TPTWB_1D
!----------------------------

!TWB_FROM_TPQ
!------------
function TWB_FROM_TPQ_0D(PT, PP, PQ) RESULT(PTWB)
!arguments and results
REAL, INTENT(IN) :: PT !air temperature (K)
REAL, INTENT(IN) :: PQ !mixing ratio (kg/kg)
REAL, INTENT(IN) :: PP !atmos. pressure (Pa)
REAL :: PTWB !Wet Bulb Temp. (K)
!local variable
REAL :: ZTD !Dew Point Temp. (K)
REAL :: ZTWBINF, ZTWBSUP, ZRV
INTEGER :: JITER
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PSYCHRO:TWB_FROM_TPQ_0D',0,ZHOOK_HANDLE)
JITER = 1
ZTD = TD_FROM_TQ(PT, PQ) 
!initial guess
ZTWBSUP = PT
ZTWBINF = ZTD
PTWB = 0.5 * (ZTWBSUP + ZTWBINF)
DO WHILE (ZTWBSUP - ZTWBINF > 0.001 .OR. JITER .LE. 50)
   ZRV = RV_FROM_TPTWB(PT, PP, PTWB)
   IF (ZRV .GT. PQ/(1 - PQ)) THEN
           ZTWBSUP = PTWB
   ELSE
           ZTWBINF = PTWB
   ENDIF
   PTWB = 0.5 * (ZTWBINF + ZTWBSUP)
   JITER = JITER + 1
END DO
IF (LHOOK) CALL DR_HOOK('MODE_PSYCHRO:TWB_FROM_TPQ_0D',1,ZHOOK_HANDLE)
end function TWB_FROM_TPQ_0D

function TWB_FROM_TPQ_1D(PT, PP, PQ) RESULT(PTWB)
!arguments and results
REAL, DIMENSION(:), INTENT(IN) :: PT !air temperature (K)
REAL, DIMENSION(:), INTENT(IN) :: PQ !humidity content (kg/kg)
REAL, DIMENSION(:), INTENT(IN) :: PP !atmos. pressure (Pa)
REAL, DIMENSION(SIZE(PT)) :: PTWB !Wet Bulb Temp. (K)
!local variable
REAL, DIMENSION(SIZE(PT)) :: ZTD !Dew Point Temp. (K)
REAL, DIMENSION(SIZE(PT)) :: ZTWBINF, ZTWBSUP, ZRV
INTEGER :: JITER, JI
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PSYCHRO:TWB_FROM_TPQ_1D',0,ZHOOK_HANDLE)
ZTD = TD_FROM_TQ(PT, PQ) 
!initial guess
ZTWBSUP = PT
ZTWBINF = ZTD
PTWB = 0.5 * (ZTWBSUP + ZTWBINF)
DO JI=1,SIZE(PT)
   JITER = 1
   DO WHILE (ZTWBSUP(JI) - ZTWBINF(JI) > 0.001 .OR. JITER .LE. 50)
      ZRV(JI) = RV_FROM_TPTWB(PT(JI), PP(JI), PTWB(JI))
      IF (ZRV(JI) .GT. PQ(JI)/(1 - PQ(JI))) THEN
              ZTWBSUP(JI) = PTWB(JI)
      ELSE
              ZTWBINF(JI) = PTWB(JI)
      ENDIF
      PTWB(JI) = 0.5 * (ZTWBINF(JI) + ZTWBSUP(JI))
   END DO
END DO
IF (LHOOK) CALL DR_HOOK('MODE_PSYCHRO:TWB_FROM_TPQ_1D',1,ZHOOK_HANDLE)
end function TWB_FROM_TPQ_1D
!-------------------------------------------------------------------------------
!
!     ######################################
      FUNCTION ENTH_FN_T_Q(PT,PQ) RESULT(PENTH)
!     ######################################
!
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the enthalpy function
!       of temperature and humidity content
!      
!
!!**  METHOD
!!    ------
!!  
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------  
!!      
!!    REFERENCE
!!    ---------
!!      
!!
!!
!!    AUTHOR
!!    ------
!!      
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    12/04/11 
!!      A. Alias    01/2013   compi. on Bull : must be 1.0E-5 instead of 1.0D-5
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, INTENT(IN)                :: PT     ! Temperature (K)
REAL, INTENT(IN)                :: PQ     ! Humidity content (kg/kg)
REAL                           :: PENTH  ! Enthalpy (J/kg)
!
!*       0.2   Declarations of local variables
!
REAL        :: ZT                          ! Temperature (C)
REAL        :: ZRV                         ! Mixing ratio (kg/kg_da)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
      IF (LHOOK) CALL DR_HOOK('MODE_PSYCHRO:ENTH_FN_T_Q',0,ZHOOK_HANDLE)
! calculate enthalpy
      ZT = PT - 273.15
      ZRV=MAX(PQ/(1-PQ),1.0E-5)
      PENTH=1.00484d3*ZT+ZRV*(2.50094d6+1.85895d3*ZT)
!
      IF (LHOOK) CALL DR_HOOK('MODE_PSYCHRO:ENTH_FN_T_Q',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END FUNCTION ENTH_FN_T_Q
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!     ######################################
      FUNCTION Q_FN_T_ENTH(PT,PENTH) RESULT(PQ)
!     ######################################
!
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the humidity content
!       as a function of temperature and enthalpy
!      
!
!!**  METHOD
!!    ------
!!  
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------  
!!      
!!    REFERENCE
!!    ---------
!!      
!!
!!
!!    AUTHOR
!!    ------
!!      
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    12/04/11 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, INTENT(IN)                :: PT     ! Temperature (K)
REAL, INTENT(IN)                :: PENTH  ! Enthalpy (J/kg)
REAL                           :: PQ     ! Humidity content (kg/kg)
!
!*       0.2   Declarations of local variables
!
REAL        :: ZT                          ! Temperature (C)
REAL        :: ZRV                         ! Mixing ratio (kg/kg_da)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
      IF (LHOOK) CALL DR_HOOK('MODE_PSYCHRO:Q_FN_T_ENTH',0,ZHOOK_HANDLE)
!
      ZT = PT - 273.15
!
!    calculate mixing ratio
      ZRV=(PENTH-1.00484d3*ZT)/(2.50094d6+1.85895d3*ZT)
!
!    validity test
      IF (ZRV < 0.0d0) THEN
        ZRV=1.d-5
      ENDIF
!
!    calculate humidity content
      PQ = ZRV/(1+ZRV)
!
     IF (LHOOK) CALL DR_HOOK('MODE_PSYCHRO:Q_FN_T_ENTH',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END FUNCTION Q_FN_T_ENTH

END MODULE MODE_PSYCHRO
