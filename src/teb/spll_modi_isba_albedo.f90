!auto_modi:spll_isba_albedo.D
MODULE MODI_ISBA_ALBEDO
INTERFACE
      SUBROUTINE ISBA_ALBEDO(PEK, OTR_ML, OMEB, PDIR_SW, PSCA_SW, PSW_BANDS, KSW, &
                             PFALB, PFFV, PFFG, PGLOBAL_SW,           &
                             PMEB_SCA_SW, PALBNIR_TVEG, PALBVIS_TVEG,               &
                             PALBNIR_TSOIL, PALBVIS_TSOIL               )
USE MODD_ISBA_n, ONLY : ISBA_PE_t
IMPLICIT NONE
LOGICAL,              INTENT(IN)   :: OTR_ML
LOGICAL,              INTENT(IN)   :: OMEB        ! True = patch with multi-energy balance 
REAL, DIMENSION(:,:), INTENT(IN)   :: PDIR_SW            ! direct incoming solar radiation
REAL, DIMENSION(:,:), INTENT(IN)   :: PSCA_SW            ! diffus incoming solar radiation
REAL, DIMENSION(:)  , INTENT(IN)   :: PSW_BANDS          ! mean wavelength of each shortwave band (m)
INTEGER,              INTENT(IN)   :: KSW                ! number of short-wave spectral bands
TYPE(ISBA_PE_t), INTENT(INOUT) :: PEK
REAL, DIMENSION(:)  , INTENT(IN)   :: PFALB              ! Floodplain albedo
REAL, DIMENSION(:)  , INTENT(IN)   :: PFFV               ! Floodplain fraction over vegetation
REAL, DIMENSION(:)  , INTENT(IN)   :: PFFG               ! Floodplain fraction over the ground
REAL, DIMENSION(:)  , INTENT(OUT)  :: PGLOBAL_SW         ! global incoming SW rad.
REAL, DIMENSION(:)  , INTENT(OUT)  :: PMEB_SCA_SW        ! diffuse incoming SW rad.
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBNIR_TVEG       ! nearIR  veg tot albedo
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBVIS_TVEG       ! visible veg tot albedo
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBNIR_TSOIL      ! nearIR  soil tot albedo
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBVIS_TSOIL      ! visible soil tot albedo
END SUBROUTINE ISBA_ALBEDO
END INTERFACE
END MODULE MODI_ISBA_ALBEDO
