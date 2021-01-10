!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #####################
      MODULE MODD_CH_WATFLUX_n
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

TYPE CH_WATFLUX_t
!
  CHARACTER(LEN=6)              :: CCH_DRY_DEP              ! deposition scheme
  REAL, DIMENSION(:,:), POINTER :: XDEP                     ! final dry deposition  
                                                            ! velocity  for lakes
  TYPE(SV_t) :: SVW

  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CCH_NAMES      ! NAME OF CHEMICAL SPECIES
                                                            ! (FOR DIAG ONLY)
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CDSTNAMES
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CSLTNAMES
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CAER_NAMES
!
END TYPE CH_WATFLUX_t
!
CONTAINS
!
SUBROUTINE CH_WATFLUX_INIT(YCH_WATFLUX)
TYPE(CH_WATFLUX_t), INTENT(INOUT) :: YCH_WATFLUX
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_WATFLUX_N:CH_WATFLUX_INIT",0,ZHOOK_HANDLE)
NULLIFY(YCH_WATFLUX%XDEP)
NULLIFY(YCH_WATFLUX%CCH_NAMES)
NULLIFY(YCH_WATFLUX%CAER_NAMES)
NULLIFY(YCH_WATFLUX%CDSTNAMES)
NULLIFY(YCH_WATFLUX%CSLTNAMES)
YCH_WATFLUX%CCH_DRY_DEP=' '
CALL SV_INIT(YCH_WATFLUX%SVW)
IF (LHOOK) CALL DR_HOOK("MODD_CH_WATFLUX_N:CH_WATFLUX_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE CH_WATFLUX_INIT


END MODULE MODD_CH_WATFLUX_n
