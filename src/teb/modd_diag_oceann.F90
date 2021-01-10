!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ######################
      MODULE MODD_DIAG_OCEAN_n
!     ######################
!
!!****  *MODD_DIAG_OCEAN - declaration of diagnostics for ocean
!!
!!    PURPOSE
!!    -------
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
!!      C. Lebeaupin Brossier   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       02/2008
!
!*       0.   DECLARATIONS
!             ------------
!
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE DIAG_OCEAN_t
!------------------------------------------------------------------------------
!
  LOGICAL :: LDIAG_OCEAN                    ! flag for ocean variables
!
!* averaged variables
!
  REAL, POINTER, DIMENSION(:)   :: XTOCMOY  ! mean temperature in CMO (°C)
  REAL, POINTER, DIMENSION(:)   :: XSOCMOY  ! mean salinity in CMO (psu)
  REAL, POINTER, DIMENSION(:)   :: XUOCMOY  ! mean zonal current in CMO (m/s)
  REAL, POINTER, DIMENSION(:)   :: XVOCMOY  ! mean meridian current in CMO (m/s)
  REAL, POINTER, DIMENSION(:)   :: XDOCMOY  ! mean density in CMO (kg/m3)
!------------------------------------------------------------------------------
!

END TYPE DIAG_OCEAN_t



CONTAINS

!





SUBROUTINE DIAG_OCEAN_INIT(YDIAG_OCEAN)
TYPE(DIAG_OCEAN_t), INTENT(INOUT) :: YDIAG_OCEAN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_OCEAN_N:DIAG_OCEAN_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YDIAG_OCEAN%XTOCMOY)
  NULLIFY(YDIAG_OCEAN%XSOCMOY)
  NULLIFY(YDIAG_OCEAN%XUOCMOY)
  NULLIFY(YDIAG_OCEAN%XVOCMOY)
  NULLIFY(YDIAG_OCEAN%XDOCMOY)
YDIAG_OCEAN%LDIAG_OCEAN=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_OCEAN_N:DIAG_OCEAN_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_OCEAN_INIT


END MODULE MODD_DIAG_OCEAN_n
