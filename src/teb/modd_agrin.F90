!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ##################
      MODULE MODD_AGRI_n
!     ##################
!
!!****  *MODD_AGRI_n - declaration of SEEDING date for summer crops 
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
!!      P. LE MOIGNE   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       06/2006
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
TYPE AGRI_t
!                                          
INTEGER, POINTER, DIMENSION (:)   :: NIRRINUM       
                                        ! Stage for Irrigation (4 stages)
!
LOGICAL, POINTER,DIMENSION(:)     :: LIRRIGATE 
                                        ! True if irrigation performed
!
LOGICAL, POINTER,DIMENSION(:)     :: LIRRIDAY 
                                        ! True if irrigation occurs during present day
!                                          
REAL, POINTER, DIMENSION(:)       :: XTHRESHOLDSPT 
                                        ! Spatialized threshold

END TYPE AGRI_t
!-------------------------------------------------------------------------------
!
TYPE AGRI_NP_t
!
TYPE(AGRI_t), DIMENSION(:), POINTER :: AL=>NULL()
!
END TYPE AGRI_NP_t
!
CONTAINS
!
!
SUBROUTINE AGRI_INIT(AG)
TYPE(AGRI_t), INTENT(INOUT) :: AG
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_AGRI_N:AGRI_INIT",0,ZHOOK_HANDLE)
  NULLIFY(AG%NIRRINUM)
  NULLIFY(AG%LIRRIGATE)
  NULLIFY(AG%LIRRIDAY)
  NULLIFY(AG%XTHRESHOLDSPT)
IF (LHOOK) CALL DR_HOOK("MODD_AGRI_N:AGRI_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE AGRI_INIT
!
SUBROUTINE AGRI_NP_INIT(YNAG,KPATCH)
TYPE(AGRI_NP_t), INTENT(INOUT) :: YNAG
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_AGRI_N:AGRI_NP_INIT",0,ZHOOK_HANDLE)
! 
IF (ASSOCIATED(YNAG%AL)) THEN
  DO JP = 1,KPATCH
    CALL AGRI_INIT(YNAG%AL(JP))
  ENDDO
  DEALLOCATE(YNAG%AL)
ELSE
  ALLOCATE(YNAG%AL(KPATCH))
  DO JP = 1,KPATCH
    CALL AGRI_INIT(YNAG%AL(JP))
  ENDDO
ENDIF
!
IF (LHOOK) CALL DR_HOOK("MODD_AGRI_N:AGRI_NP_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE AGRI_NP_INIT
!
!-------------------------------------------------------------------------------
!
END MODULE MODD_AGRI_n
