!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL version 2.1
!SFX_LIC version 1. See LICENSE, Licence_CeCILL_V2.1-en.txt and Licence_CeCILL_V2.1-fr.txt  
!SFX_LIC for details. version 1.
!     #####################
      MODULE MODD_CH_ISBA_n
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

TYPE CH_ISBA_t
!
  CHARACTER(LEN=28)  :: CCHEM_SURF_FILE  ! name of general (chemical) purpose ASCII input file
  CHARACTER(LEN=6)                :: CCH_DRY_DEP            !  deposition scheme
  REAL, DIMENSION(:,:), POINTER :: XDEP                   ! final dry deposition  
                                                            ! velocity  for nature
  REAL, DIMENSION(:),   POINTER :: XSOILRC_SO2            ! for SO2
  REAL, DIMENSION(:),   POINTER :: XSOILRC_O3             ! for O3                                                            
  LOGICAL                         :: LCH_BIO_FLUX           ! flag for the calculation of
                                                            ! biogenic fluxes
  LOGICAL                         :: LCH_NO_FLUX            ! flag for the calculation of
                                                            ! biogenic NO fluxes
  TYPE(SV_t) :: SVI

  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CCH_NAMES      ! NAME OF CHEMICAL SPECIES
                                                            ! (FOR DIAG ONLY)
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CAER_NAMES     ! NAME OF CHEMICAL SPECIES
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CDSTNAMES      ! NAME OF CHEMICAL SPECIES
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CSLTNAMES      ! NAME OF CHEMICAL SPECIES                                                            
!
END TYPE CH_ISBA_t
!
TYPE CH_ISBA_NP_t
!
TYPE(CH_ISBA_t), DIMENSION(:), POINTER :: AL=>NULL()
!
END TYPE CH_ISBA_NP_t
!
CONTAINS
!
SUBROUTINE CH_ISBA_INIT(YCH_ISBA)
TYPE(CH_ISBA_t), INTENT(INOUT) :: YCH_ISBA
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_ISBA_N:CH_ISBA_INIT",0,ZHOOK_HANDLE)
NULLIFY(YCH_ISBA%XDEP)
NULLIFY(YCH_ISBA%XSOILRC_SO2)
NULLIFY(YCH_ISBA%XSOILRC_O3)
NULLIFY(YCH_ISBA%CCH_NAMES)
NULLIFY(YCH_ISBA%CAER_NAMES)
NULLIFY(YCH_ISBA%CDSTNAMES)
NULLIFY(YCH_ISBA%CSLTNAMES)
YCH_ISBA%CCHEM_SURF_FILE=' '
YCH_ISBA%CCH_DRY_DEP=' '
YCH_ISBA%LCH_BIO_FLUX=.FALSE.
YCH_ISBA%LCH_NO_FLUX=.FALSE.
CALL SV_INIT(YCH_ISBA%SVI)
IF (LHOOK) CALL DR_HOOK("MODD_CH_ISBA_N:CH_ISBA_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE CH_ISBA_INIT
!
SUBROUTINE CH_ISBA_NP_INIT(YCH_ISBA_NP,KPATCH)
TYPE(CH_ISBA_NP_t), INTENT(INOUT) :: YCH_ISBA_NP
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_ISBA_N:CH_ISBA_NP_INIT",0,ZHOOK_HANDLE)
!
IF (ASSOCIATED(YCH_ISBA_NP%AL)) THEN
  DO JP = 1,KPATCH
    CALL CH_ISBA_INIT(YCH_ISBA_NP%AL(JP))
  ENDDO        
  DEALLOCATE(YCH_ISBA_NP%AL)
ELSE
  ALLOCATE(YCH_ISBA_NP%AL(KPATCH))
  DO JP = 1,KPATCH
    CALL CH_ISBA_INIT(YCH_ISBA_NP%AL(JP))
  ENDDO
ENDIF
!
IF (LHOOK) CALL DR_HOOK("MODD_CH_ISBA_N:CH_ISBA_NP_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE CH_ISBA_NP_INIT
!
END MODULE MODD_CH_ISBA_n
