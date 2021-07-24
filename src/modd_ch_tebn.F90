!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL version 2.1
!SFX_LIC version 1. See LICENSE, Licence_CeCILL_V2.1-en.txt and Licence_CeCILL_V2.1-fr.txt  
!SFX_LIC for details. version 1.
!     #####################
      MODULE MODD_CH_TEB_n
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

TYPE CH_TEB_t
!
  CHARACTER(LEN=28)  :: CCHEM_SURF_FILE  ! name of general (chemical) purpose ASCII input file
  CHARACTER(LEN=6)                :: CCH_DRY_DEP            !  deposition scheme
  REAL, DIMENSION(:,:),   POINTER :: XDEP                   ! final dry deposition  
                                                            ! velocity  for nature
  REAL, DIMENSION(:,:),   POINTER :: XSOILRC_SO2            ! for SO2
  REAL, DIMENSION(:,:),   POINTER :: XSOILRC_O3             ! for O3                                                            
  LOGICAL                         :: LCH_BIO_FLUX           ! flag for the calculation of
                                                            ! biogenic fluxes
  LOGICAL                         :: LCH_NO_FLUX            ! flag for the calculation of
                                                            ! biogenic NO fluxes
  TYPE(SV_t) :: SVT
  
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CCH_NAMES      ! NAME OF CHEMICAL SPECIES
                                                            ! (FOR DIAG ONLY)
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CAER_NAMES     ! NAME OF CHEMICAL SPECIES
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CDSTNAMES      ! NAME OF CHEMICAL SPECIES
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CSLTNAMES      ! NAME OF CHEMICAL SPECIES
!
END TYPE CH_TEB_t



CONTAINS

!





SUBROUTINE CH_TEB_INIT(YCH_TEB)
TYPE(CH_TEB_t), INTENT(INOUT) :: YCH_TEB
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_TEB_N:CH_TEB_INIT",0,ZHOOK_HANDLE)
NULLIFY(YCH_TEB%XDEP)
NULLIFY(YCH_TEB%XSOILRC_SO2)
NULLIFY(YCH_TEB%XSOILRC_O3)
NULLIFY(YCH_TEB%CCH_NAMES)
NULLIFY(YCH_TEB%CAER_NAMES)
NULLIFY(YCH_TEB%CDSTNAMES)
NULLIFY(YCH_TEB%CSLTNAMES)
YCH_TEB%CCHEM_SURF_FILE=' '
YCH_TEB%CCH_DRY_DEP=' '
YCH_TEB%LCH_BIO_FLUX=.FALSE.
YCH_TEB%LCH_NO_FLUX=.FALSE.
CALL SV_INIT(YCH_TEB%SVT)
IF (LHOOK) CALL DR_HOOK("MODD_CH_TEB_N:CH_TEB_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE CH_TEB_INIT


END MODULE MODD_CH_TEB_n
