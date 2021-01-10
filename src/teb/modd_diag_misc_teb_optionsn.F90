!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ############################
      MODULE MODD_DIAG_MISC_TEB_OPTIONS_n
!     ############################
!
!!****  *MODD_DIAG_MISC_TEB - declaration of packed surface parameters for TEB scheme
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
!!      C de Munck        02/13  adding runoff contributions for teb garden  
!!      V. Masson      06/2013 splits module in two
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

TYPE DIAG_MISC_TEB_OPTIONS_t
!------------------------------------------------------------------------------
!
  LOGICAL :: LSURF_MISC_BUDGET   ! flag for miscellaneous terms of teb scheme
  LOGICAL :: LSURF_EVAP_BUDGET   ! flag for all terms of evaporation
  LOGICAL :: LSURF_DIAG_ALBEDO   ! flag to write out diagnostic albedo
!
END TYPE DIAG_MISC_TEB_OPTIONS_t
!
CONTAINS
!
SUBROUTINE DIAG_MISC_TEB_OPTIONS_INIT(YDIAG_MISC_TEB_OPTIONS)
TYPE(DIAG_MISC_TEB_OPTIONS_t), INTENT(INOUT) :: YDIAG_MISC_TEB_OPTIONS
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_TEB_N:DIAG_MISC_TEB_OPTIONS_INIT",0,ZHOOK_HANDLE)
YDIAG_MISC_TEB_OPTIONS%LSURF_MISC_BUDGET=.FALSE.
YDIAG_MISC_TEB_OPTIONS%LSURF_EVAP_BUDGET=.FALSE.
YDIAG_MISC_TEB_OPTIONS%LSURF_DIAG_ALBEDO=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_TEB_N:DIAG_MISC_TEB_OPTIONS_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_MISC_TEB_OPTIONS_INIT


END MODULE MODD_DIAG_MISC_TEB_OPTIONS_n
