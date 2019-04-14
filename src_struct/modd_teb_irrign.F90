!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ################
      MODULE MODD_TEB_IRRIG_n
!     ################
!
!!****  *MODD_TEB_IRRIG_n - declaration of surface parameters for urban canopy
!!
!!    PURPOSE
!!    -------
!     Declaration of surface parameters
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       07/2006
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE TEB_IRRIG_t
!
  LOGICAL                        :: LTEB_IRRIG        ! flag to use irrigation for gardens or greenroofs
  LOGICAL                        :: LPAR_GD_IRRIG     ! flag to use prescribed irrigation for gardens
  LOGICAL                        :: LPAR_GR_IRRIG     ! flag to use prescribed irrigation for greenroofs
  LOGICAL                        :: LPAR_RD_IRRIG     ! flag to use prescribed irrigation for roads
  REAL,    POINTER, DIMENSION(:) :: XGD_START_MONTH   ! gardens : start month for irrigation (included)
  REAL,    POINTER, DIMENSION(:) :: XGD_END_MONTH     ! gardens : end   month for irrigation (included)
  REAL,    POINTER, DIMENSION(:) :: XGD_START_HOUR    ! gardens : start solar hour for irrigation (included, hour)
  REAL,    POINTER, DIMENSION(:) :: XGD_END_HOUR      ! gardens : end   solar hour for irrigation (excluded, hour)
  REAL,    POINTER, DIMENSION(:) :: XGD_24H_IRRIG     ! gardens : total irrigation over 24 hours (kg/m2)
  REAL,    POINTER, DIMENSION(:) :: XGR_START_MONTH   ! greenroofs : start month for irrigation (included)
  REAL,    POINTER, DIMENSION(:) :: XGR_END_MONTH     ! greenroofs : end   month for irrigation (included)
  REAL,    POINTER, DIMENSION(:) :: XGR_START_HOUR    ! greenroofs : start solar hour for irrigation (included, hour)
  REAL,    POINTER, DIMENSION(:) :: XGR_END_HOUR      ! greenroofs : end   solar hour for irrigation (excluded, hour)
  REAL,    POINTER, DIMENSION(:) :: XGR_24H_IRRIG     ! greenroofs : total irrigation over 24 hours (kg/m2)
  REAL,    POINTER, DIMENSION(:) :: XRD_START_MONTH   ! roads : start month for irrigation (included)
  REAL,    POINTER, DIMENSION(:) :: XRD_END_MONTH     ! roads : end   month for irrigation (included)
  REAL,    POINTER, DIMENSION(:) :: XRD_START_HOUR    ! roads : start solar hour for irrigation (included, hour)
  REAL,    POINTER, DIMENSION(:) :: XRD_END_HOUR      ! roads : end   solar hour for irrigation (excluded, hour)
  REAL,    POINTER, DIMENSION(:) :: XRD_24H_IRRIG     ! roads : total irrigation over 24 hours (kg/m2)
!
END TYPE TEB_IRRIG_t



CONTAINS

!




SUBROUTINE TEB_IRRIG_INIT(YTEB_IRRIG)
TYPE(TEB_IRRIG_t), INTENT(INOUT) :: YTEB_IRRIG
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_IRRIG_N:TEB_IRRIG_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YTEB_IRRIG%XGD_START_MONTH)
  NULLIFY(YTEB_IRRIG%XGD_END_MONTH)
  NULLIFY(YTEB_IRRIG%XGD_START_HOUR)
  NULLIFY(YTEB_IRRIG%XGD_END_HOUR)
  NULLIFY(YTEB_IRRIG%XGD_24H_IRRIG)
  NULLIFY(YTEB_IRRIG%XGR_START_MONTH)
  NULLIFY(YTEB_IRRIG%XGR_END_MONTH)
  NULLIFY(YTEB_IRRIG%XGR_START_HOUR)
  NULLIFY(YTEB_IRRIG%XGR_END_HOUR)
  NULLIFY(YTEB_IRRIG%XGR_24H_IRRIG)
  NULLIFY(YTEB_IRRIG%XRD_START_MONTH)
  NULLIFY(YTEB_IRRIG%XRD_END_MONTH)
  NULLIFY(YTEB_IRRIG%XRD_START_HOUR)
  NULLIFY(YTEB_IRRIG%XRD_END_HOUR)
  NULLIFY(YTEB_IRRIG%XRD_24H_IRRIG)
YTEB_IRRIG%LTEB_IRRIG = .FALSE.
YTEB_IRRIG%LPAR_GD_IRRIG = .FALSE.
YTEB_IRRIG%LPAR_GR_IRRIG = .FALSE.
YTEB_IRRIG%LPAR_RD_IRRIG = .FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_TEB_IRRIG_N:TEB_IRRIG_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_IRRIG_INIT


END MODULE MODD_TEB_IRRIG_n
