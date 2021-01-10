!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE ISBA_SNOW_FRAC(HSNOW, PWSNOW, PRSNOW, PASNOW,       &
                                PVEG, PLAI, PZ0, PPSN, PPSNV_A, PPSNG, PPSNV )  
!     ##########################################################################
!
!!****  *ISBA_SNOW_FRAC*  
!!
!!    PURPOSE
!!    -------
!
!     Calculates grid-averaged albedo and emissivity (according to snow scheme)
!         
!!    EXTERNAL
!!    --------
!!
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      
!!    AUTHOR
!!    ------
!!
!!      S. Belair           * Meteo-France *
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SNOW_PAR   , ONLY : XEMISSN, XEMCRIN, XSNOWDMIN, &
                             XRHOSMAX_ES, XRHOSMIN_ES, &
                             XWCRN_EXPL, XDCRN_EXPL
!
USE MODD_PREP_SNOW, ONLY : LSNOW_FRAC_TOT
!
USE MODE_SURF_SNOW_FRAC
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
 CHARACTER(LEN=*)    , INTENT(IN)   :: HSNOW      ! ISBA snow scheme
!
REAL, DIMENSION(:,:), INTENT(IN)   :: PWSNOW     ! snow reservoir
REAL, DIMENSION(:,:), INTENT(IN)   :: PRSNOW     ! snow density
REAL, DIMENSION(:)  , INTENT(IN)   :: PASNOW     ! snow albedo
REAL, DIMENSION(:)  , INTENT(IN)   :: PVEG
REAL, DIMENSION(:)  , INTENT(IN)   :: PLAI
REAL, DIMENSION(:)  , INTENT(IN)   :: PZ0
REAL, DIMENSION(:)  , INTENT(INOUT):: PPSN
REAL, DIMENSION(:)  , INTENT(INOUT):: PPSNV_A
REAL, DIMENSION(:)  , INTENT(INOUT):: PPSNG
REAL, DIMENSION(:)  , INTENT(INOUT):: PPSNV
!                                     PVEG = fraction of vegetation
!                                     PLAI = leaf area index
!                                     PZ0  = roughness length for momentum
!                                     PPSN = grid fraction covered by snow
!                                     PPSNG = fraction of the ground covered by snow 
!                                     PPSNV = fraction of the veg covered by snow 
!
!-------------------------------------------------------------------------------
!
!*      0.     Local variables
!              ---------------
!
INTEGER                          :: JLAYER
REAL, DIMENSION(SIZE(PVEG))      :: ZSNOWSWE
REAL, DIMENSION(SIZE(PVEG))      :: ZSNOWD
REAL, DIMENSION(SIZE(PVEG))      :: ZSNOWRHO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_SNOW_FRAC',0,ZHOOK_HANDLE)
!
!*      1.     Compute Total SWE (kg m-2) and snowpack average density (kg m-3)
!              ----------------------------------------------------------------
!
ZSNOWSWE(:) = 0.
!
DO JLAYER=1,SIZE(PWSNOW,2)
  ZSNOWSWE(:) = ZSNOWSWE(:) + PWSNOW(:,JLAYER)
END DO
!
IF (HSNOW == '3-L' .OR. HSNOW == 'CRO') THEN                  
   ZSNOWD  (:) = 0.
   DO JLAYER=1,SIZE(PWSNOW,2)
     ZSNOWD  (:) = ZSNOWD(:) + PWSNOW(:,JLAYER) / PRSNOW(:,JLAYER)
   END DO
   ZSNOWRHO(:) = ZSNOWSWE(:)/MAX(XSNOWDMIN,ZSNOWD(:))
   ZSNOWRHO(:) = MAX(XRHOSMIN_ES,MIN(XRHOSMAX_ES,ZSNOWRHO(:)))
ELSE
   ZSNOWRHO(:) = PRSNOW(:,1)
END IF
!
!*      2.     Snow fraction over ground
!              -------------------------
!
IF (HSNOW == 'CRO' .OR. HSNOW == '3-L') THEN
   PPSNG(:) = MIN(1.0, ZSNOWD(:)/XDCRN_EXPL)
ELSE
   PPSNG(:) = SNOW_FRAC_GROUND(ZSNOWSWE)
ENDIF
!
!*      3.     Snow fraction over vegetation
!              -----------------------------
!
IF (HSNOW == 'EBA' ) THEN
   PPSNV_A(:) = SNOW_FRAC_VEG_A (PPSNG,PLAI,PASNOW)
   PPSNV  (:) = PPSNV_A(:)
ELSE
   PPSNV  (:) = SNOW_FRAC_VEG   (PPSNG,ZSNOWSWE,PZ0,ZSNOWRHO)
ENDIF
!
!*      4.     Total snow fraction
!              -------------------
!
PPSN(:)       = SNOW_FRAC_NAT(ZSNOWSWE,PPSNG,PPSNV,PVEG)
!
IF (LSNOW_FRAC_TOT) THEN
  PPSN (:) = MIN(1.0, ZSNOWSWE(:)/XWCRN_EXPL)      
  PPSNG(:) = PPSN(:)
  PPSNV(:) = PPSN(:)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('ISBA_SNOW_FRAC',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ISBA_SNOW_FRAC
