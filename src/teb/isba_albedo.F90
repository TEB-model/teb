!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE ISBA_ALBEDO(PEK, OTR_ML, OMEB, PDIR_SW, PSCA_SW, PSW_BANDS, KSW, &
                             PFALB, PFFV, PFFG, PGLOBAL_SW,           &
                             PMEB_SCA_SW, PALBNIR_TVEG, PALBVIS_TVEG,               &
                             PALBNIR_TSOIL, PALBVIS_TSOIL               )
!     ##########################################################################
!
!!****  *ISBA_ALBEDO*  
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
!!	S. Belair           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    
!!      P. Samuelsson  02/2012  MEB
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_n, ONLY : ISBA_PE_t
!
USE MODD_SURF_PAR,     ONLY : XUNDEF
!
USE MODI_ALBEDO_FROM_NIR_VIS
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
LOGICAL,              INTENT(IN)   :: OTR_ML
LOGICAL,              INTENT(IN)   :: OMEB        ! True = patch with multi-energy balance 
!                                                 ! False = patch with classical ISBA
!
REAL, DIMENSION(:,:), INTENT(IN)   :: PDIR_SW            ! direct incoming solar radiation
REAL, DIMENSION(:,:), INTENT(IN)   :: PSCA_SW            ! diffus incoming solar radiation
REAL, DIMENSION(:)  , INTENT(IN)   :: PSW_BANDS          ! mean wavelength of each shortwave band (m)
INTEGER,              INTENT(IN)   :: KSW                ! number of short-wave spectral bands
!
TYPE(ISBA_PE_t), INTENT(INOUT) :: PEK
!
REAL, DIMENSION(:)  , INTENT(IN)   :: PFALB              ! Floodplain albedo
REAL, DIMENSION(:)  , INTENT(IN)   :: PFFV               ! Floodplain fraction over vegetation
REAL, DIMENSION(:)  , INTENT(IN)   :: PFFG               ! Floodplain fraction over the ground
!
REAL, DIMENSION(:)  , INTENT(OUT)  :: PGLOBAL_SW         ! global incoming SW rad.
REAL, DIMENSION(:)  , INTENT(OUT)  :: PMEB_SCA_SW        ! diffuse incoming SW rad.
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBNIR_TVEG       ! nearIR  veg tot albedo
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBVIS_TVEG       ! visible veg tot albedo
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBNIR_TSOIL      ! nearIR  soil tot albedo
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBVIS_TSOIL      ! visible soil tot albedo
!
!-------------------------------------------------------------------------------
!
!*      0.     Local variables
!              ---------------
!
INTEGER                          :: JLAYER
INTEGER                          :: JSWB
REAL, DIMENSION(SIZE(PEK%XALBNIR))      :: ZSW_UP
REAL, DIMENSION(SIZE(PEK%XALBNIR),KSW)  :: ZDIR_ALB_WITHOUT_SNOW
REAL, DIMENSION(SIZE(PEK%XALBNIR),KSW)  :: ZSCA_ALB_WITHOUT_SNOW
REAL, DIMENSION(SIZE(PEK%XALBNIR),KSW)  :: ZDIR_ALB_VEG_WITHOUT_SNOW
REAL, DIMENSION(SIZE(PEK%XALBNIR),KSW)  :: ZSCA_ALB_VEG_WITHOUT_SNOW
REAL, DIMENSION(SIZE(PEK%XALBNIR),KSW)  :: ZDIR_ALB_SOIL_WITHOUT_SNOW
REAL, DIMENSION(SIZE(PEK%XALBNIR),KSW)  :: ZSCA_ALB_SOIL_WITHOUT_SNOW
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*      2.     Compute snow-free albedo
!              ------------------------
!
!* Snow-free surface albedo for each wavelength
!
IF (LHOOK) CALL DR_HOOK('ISBA_ALBEDO',0,ZHOOK_HANDLE)
!
IF (OTR_ML )THEN
  IF (OMEB) THEN
    PALBNIR_TVEG (:) =               PEK%XALBNIR_VEG(:)
    PALBNIR_TSOIL(:) = ( 1.-PFFG(:))*PEK%XALBNIR_SOIL(:) + PFFG(:)*PFALB(:)   
    PALBVIS_TVEG (:) =               PEK%XALBVIS_VEG(:)
    PALBVIS_TSOIL(:) = ( 1.-PFFG(:))*PEK%XALBVIS_SOIL(:) + PFFG(:)*PFALB(:)
  ELSE
    PALBNIR_TVEG (:) = PEK%XALBNIR_VEG(:)
    PALBNIR_TSOIL(:) = PEK%XALBNIR_SOIL(:) 
    PALBVIS_TVEG (:) = PEK%XALBVIS_VEG(:)
    PALBVIS_TSOIL(:) = PEK%XALBVIS_SOIL(:) 
  ENDIF
ELSE
  PALBNIR_TVEG (:) = XUNDEF
  PALBNIR_TSOIL(:) = XUNDEF
  PALBVIS_TVEG (:) = XUNDEF
  PALBVIS_TSOIL(:) = XUNDEF
ENDIF
!
 CALL ALBEDO_FROM_NIR_VIS(PSW_BANDS, PEK%XALBNIR(:), PEK%XALBVIS(:), PEK%XALBUV(:),  &
                           ZDIR_ALB_WITHOUT_SNOW, ZSCA_ALB_WITHOUT_SNOW )  
!
!* total shortwave incoming radiation
!
PGLOBAL_SW (:) = 0.
PMEB_SCA_SW(:) = 0.
DO JSWB=1,KSW
  PGLOBAL_SW (:) = PGLOBAL_SW(:) + (PDIR_SW(:,JSWB) + PSCA_SW(:,JSWB))
  PMEB_SCA_SW(:) = PMEB_SCA_SW(:) + (PSCA_SW(:,JSWB))
END DO
!
!* snow-free global albedo (needed by ISBA)
!
ZSW_UP(:) = 0. 
DO JSWB=1,KSW
  ZSW_UP(:) =  ZSW_UP(:)                                       &
               + ZDIR_ALB_WITHOUT_SNOW(:,JSWB) * PDIR_SW(:,JSWB) &
               + ZSCA_ALB_WITHOUT_SNOW(:,JSWB) * PSCA_SW(:,JSWB)  
END DO
PEK%XSNOWFREE_ALB(:) = XUNDEF
WHERE(PGLOBAL_SW(:)>0.)  
  PEK%XSNOWFREE_ALB(:) = ZSW_UP(:) / PGLOBAL_SW(:)
ELSEWHERE
  PEK%XSNOWFREE_ALB(:) = ZDIR_ALB_WITHOUT_SNOW(:,1)
END WHERE
!
IF(PEK%TSNOW%SCHEME == 'EBA') THEN
  CALL ALBEDO_FROM_NIR_VIS(PSW_BANDS,            &
            PEK%XALBNIR_VEG(:), PEK%XALBVIS_VEG(:), PEK%XALBUV_VEG(:), &
            ZDIR_ALB_VEG_WITHOUT_SNOW, ZSCA_ALB_VEG_WITHOUT_SNOW )  
  ZSW_UP(:) = 0.
  DO JSWB=1,KSW
     ZSW_UP(:) =  ZSW_UP(:)                                           &
                  + ZDIR_ALB_VEG_WITHOUT_SNOW(:,JSWB) * PDIR_SW(:,JSWB) &
                  + ZSCA_ALB_VEG_WITHOUT_SNOW(:,JSWB) * PSCA_SW(:,JSWB)  
  END DO
  PEK%XSNOWFREE_ALB_VEG(:) = XUNDEF
  WHERE(PGLOBAL_SW(:)>0.)  PEK%XSNOWFREE_ALB_VEG(:) = ZSW_UP(:) / PGLOBAL_SW(:)
!
  CALL ALBEDO_FROM_NIR_VIS(PSW_BANDS,               &
            PEK%XALBNIR_SOIL(:), PEK%XALBVIS_SOIL(:), PEK%XALBUV_SOIL(:), &
            ZDIR_ALB_SOIL_WITHOUT_SNOW, ZSCA_ALB_SOIL_WITHOUT_SNOW    )  
  ZSW_UP(:) = 0.
  DO JSWB=1,KSW
    ZSW_UP(:) =  ZSW_UP(:)                                            &
               + ZDIR_ALB_SOIL_WITHOUT_SNOW(:,JSWB) * PDIR_SW(:,JSWB) &
               + ZSCA_ALB_SOIL_WITHOUT_SNOW(:,JSWB) * PSCA_SW(:,JSWB)  
  END DO
  PEK%XSNOWFREE_ALB_SOIL(:) = XUNDEF
  WHERE(PGLOBAL_SW(:)>0.)  PEK%XSNOWFREE_ALB_SOIL(:) = ZSW_UP(:) / PGLOBAL_SW(:)             
ENDIF
!
IF (LHOOK) CALL DR_HOOK('ISBA_ALBEDO',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ISBA_ALBEDO
