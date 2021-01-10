!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ###########################
      MODULE MODD_DIAG_UTCI_TEB_n
!     ############################
!
!!****  *MODD_DIAG_UTCI_TEB - declaration of confort from TEB scheme
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
!!      P. Le Moigne   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       07/10/04
!
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

TYPE DIAG_UTCI_TEB_t
!------------------------------------------------------------------------------
!
  LOGICAL :: LUTCI                               ! flag to compute UTCI quantities
  REAL, POINTER, DIMENSION(:)   :: XUTCI_IN      ! UTCI for person indoor
  REAL, POINTER, DIMENSION(:)   :: XUTCI_OUTSUN  ! UTCI for person outdoor at sun
  REAL, POINTER, DIMENSION(:)   :: XUTCI_OUTSHADE! UTCI for person outdoor at shade
  REAL, POINTER, DIMENSION(:)   :: XTRAD_SUN! Mean radiant temperature seen by person at sun (K)
  REAL, POINTER, DIMENSION(:)   :: XTRAD_SHADE! Mean radiant temperature seen by person in shade (K)
  REAL, POINTER, DIMENSION(:,:) :: XUTCIC_IN      ! Cumulated UTCI stress for person indoor
  REAL, POINTER, DIMENSION(:,:) :: XUTCIC_OUTSUN  ! Cumulated UTCI stress for person outdoor at sun
  REAL, POINTER, DIMENSION(:,:) :: XUTCIC_OUTSHADE! Cumulated UTCI stress for person outdoor at shade
!
END TYPE DIAG_UTCI_TEB_t
!



CONTAINS

!




SUBROUTINE DIAG_UTCI_TEB_INIT(YDIAG_UTCI_TEB)
TYPE(DIAG_UTCI_TEB_t), INTENT(INOUT) :: YDIAG_UTCI_TEB
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_UTCI_TEB_N:DIAG_UTCI_TEB_INIT",0,ZHOOK_HANDLE)
YDIAG_UTCI_TEB%LUTCI=.FALSE.
  NULLIFY(YDIAG_UTCI_TEB%XUTCI_IN)
  NULLIFY(YDIAG_UTCI_TEB%XUTCI_OUTSUN)
  NULLIFY(YDIAG_UTCI_TEB%XUTCI_OUTSHADE)
  NULLIFY(YDIAG_UTCI_TEB%XTRAD_SUN)
  NULLIFY(YDIAG_UTCI_TEB%XTRAD_SHADE)
  NULLIFY(YDIAG_UTCI_TEB%XUTCIC_IN)
  NULLIFY(YDIAG_UTCI_TEB%XUTCIC_OUTSUN)
  NULLIFY(YDIAG_UTCI_TEB%XUTCIC_OUTSHADE)
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_UTCI_TEB_N:DIAG_UTCI_TEB_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_UTCI_TEB_INIT



END MODULE MODD_DIAG_UTCI_TEB_n
