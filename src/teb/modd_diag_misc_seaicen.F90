!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ######################
      MODULE MODD_DIAG_MISC_SEAICE_n
!     ######################
!
!!****  *MODD_DIAG_SEAICE - declaration of diagnostics for sea ice model 
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
!!      S.Sénési   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       08/2013
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

TYPE DIAG_MISC_SEAICE_t
!------------------------------------------------------------------------------
!
  LOGICAL :: LDIAG_MISC_SEAICE                    ! flag for SeaIce variables
!
!* averaged variables
!
  REAL, POINTER, DIMENSION(:)   :: XSIT  ! Sea ice thickness 
  REAL, POINTER, DIMENSION(:)   :: XSND  ! Sea ice snow depth 
  REAL, POINTER, DIMENSION(:)   :: XMLT  ! Sea mixed layer temp. seen by Gelato 
!------------------------------------------------------------------------------
!

END TYPE DIAG_MISC_SEAICE_t

CONTAINS

!





SUBROUTINE DIAG_MISC_SEAICE_INIT(YDIAG_MISC_SEAICE)
TYPE(DIAG_MISC_SEAICE_t), INTENT(INOUT) :: YDIAG_MISC_SEAICE
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_SEAICE_N:DIAG_MISC_SEAICE_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YDIAG_MISC_SEAICE%XSIT)
  NULLIFY(YDIAG_MISC_SEAICE%XSND)
  NULLIFY(YDIAG_MISC_SEAICE%XMLT)
YDIAG_MISC_SEAICE%LDIAG_MISC_SEAICE=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_SEAICE_N:DIAG_MISC_SEAICE_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_MISC_SEAICE_INIT


END MODULE MODD_DIAG_MISC_SEAICE_n
