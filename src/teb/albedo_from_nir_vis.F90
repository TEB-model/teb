!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE ALBEDO_FROM_NIR_VIS(PSW_BANDS,PALBNIR,PALBVIS,PALBUV,PDIR_ALB,PSCA_ALB)
!     ###########################################################################
!
!!****  *ALBEDO_FROM_NIR_VIS* - routine to initialize albedo for
!!                              any wavelength from near-infra-red,
!!                              visible and UV albedo
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/2003 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ISBA_PAR,   ONLY : XRED_EDGE, XUV_EDGE
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
REAL, DIMENSION(:),   INTENT(IN) :: PSW_BANDS ! mean wavelength of each shortwave band (m)
REAL, DIMENSION(:),   INTENT(IN) :: PALBNIR   ! near infra-red albedo
REAL, DIMENSION(:),   INTENT(IN) :: PALBVIS   ! visible albedo
REAL, DIMENSION(:),   INTENT(IN) :: PALBUV    ! UV albedo
REAL, DIMENSION(:,:), INTENT(OUT):: PDIR_ALB  ! direct  albedo for each wavelength
REAL, DIMENSION(:,:), INTENT(OUT):: PSCA_ALB  ! diffuse albedo for each wavelength
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: ISWB ! number of SW spectral bands
INTEGER :: JSWB ! loop counter on number of SW spectral bands
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ALBEDO_FROM_NIR_VIS',0,ZHOOK_HANDLE)
ISWB = SIZE(PSW_BANDS)
!
PDIR_ALB(:,:) = XUNDEF
PSCA_ALB(:,:) = XUNDEF
!
IF (ISWB==1) THEN
  WHERE(PALBNIR(:)/= XUNDEF) PDIR_ALB(:,1) = 0.5*(PALBNIR(:)+PALBVIS(:))
ELSE 
  DO JSWB=1,ISWB
    IF (PSW_BANDS(JSWB)>XRED_EDGE ) THEN      ! XRED_EDGE=0.7 micro-m
      PDIR_ALB(:,JSWB) = PALBNIR(:)
    ELSE IF (PSW_BANDS(JSWB)<XUV_EDGE ) THEN  ! XUV_EDGE=0.25 micro-m
      PDIR_ALB(:,JSWB) = PALBUV (:)
    ELSE
      PDIR_ALB(:,JSWB) = PALBVIS(:)
    END IF
  END DO
END IF
!
PSCA_ALB(:,:) = PDIR_ALB(:,:)
IF (LHOOK) CALL DR_HOOK('ALBEDO_FROM_NIR_VIS',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE ALBEDO_FROM_NIR_VIS
