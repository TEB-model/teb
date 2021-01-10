!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #####################
      MODULE MODD_CH_SEAFLUX_n
!     ######################
!
!!
!!    PURPOSE
!!    -------
!     
!   
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!
!!    AUTHOR
!!    ------
!!  P. Tulet   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!  16/07/03 (P. Tulet)  restructured for externalization
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_SV_n, ONLY : SV_t, SV_INIT
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE CH_SEAFLUX_t
!
  CHARACTER(LEN=6)              :: CCH_DRY_DEP              !  deposition scheme
  REAL, DIMENSION(:,:), POINTER :: XDEP                     ! final dry deposition  
                                                            ! velocity  for sea
  TYPE(SV_t) :: SVS  

  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CCH_NAMES      ! NAME OF CHEMICAL SPECIES
                                                            ! (FOR DIAG ONLY)
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CDSTNAMES
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CSLTNAMES
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CAER_NAMES                                                            
!
END TYPE CH_SEAFLUX_t



CONTAINS

!




SUBROUTINE CH_SEAFLUX_INIT(YCH_SEAFLUX)
TYPE(CH_SEAFLUX_t), INTENT(INOUT) :: YCH_SEAFLUX
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_SEAFLUX_N:CH_SEAFLUX_INIT",0,ZHOOK_HANDLE)
NULLIFY(YCH_SEAFLUX%XDEP)
NULLIFY(YCH_SEAFLUX%CCH_NAMES)
NULLIFY(YCH_SEAFLUX%CAER_NAMES)
NULLIFY(YCH_SEAFLUX%CDSTNAMES)
NULLIFY(YCH_SEAFLUX%CSLTNAMES)
YCH_SEAFLUX%CCH_DRY_DEP=' '
CALL SV_INIT(YCH_SEAFLUX%SVS)
IF (LHOOK) CALL DR_HOOK("MODD_CH_SEAFLUX_N:CH_SEAFLUX_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE CH_SEAFLUX_INIT


END MODULE MODD_CH_SEAFLUX_n
