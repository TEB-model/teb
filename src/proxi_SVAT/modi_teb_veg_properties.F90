!auto_modi:spll_teb_veg_properties.D
MODULE MODI_TEB_VEG_PROPERTIES
INTERFACE
      SUBROUTINE TEB_VEG_PROPERTIES (PMASK, IO, PEK, PDIR_SW, PSCA_SW, PSW_BANDS, KSW, &
                                     PTS, PEMIS, PALB, PTA, PALBNIR_TVEG, PALBVIS_TVEG,&
                                     PALBNIR_TSOIL, PALBVIS_TSOIL      )  
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_PE_t
IMPLICIT NONE
REAL, DIMENSION(:), INTENT(IN) :: PMASK
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_PE_t), INTENT(INOUT) :: PEK
REAL, DIMENSION(:,:), INTENT(IN)   :: PDIR_SW            ! direct incoming solar radiation
REAL, DIMENSION(:,:), INTENT(IN)   :: PSCA_SW            ! diffus incoming solar radiation
REAL, DIMENSION(:)  , INTENT(IN)   :: PSW_BANDS          ! mean wavelength of each shortwave band (m)
INTEGER,              INTENT(IN)   :: KSW                ! number of short-wave spectral bands
REAL, DIMENSION(:)  , INTENT(OUT)  :: PTS                ! radiative surface temperature
REAL, DIMENSION(:)  , INTENT(OUT)  :: PEMIS              ! green areas emissivity
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALB               ! green areas albedo
REAL, DIMENSION(:)  , INTENT(IN), OPTIONAL :: PTA        ! Air temperature (K)
REAL, DIMENSION(:)  , INTENT(OUT), OPTIONAL  :: PALBNIR_TVEG       ! nearIR  veg tot albedo
REAL, DIMENSION(:)  , INTENT(OUT), OPTIONAL  :: PALBVIS_TVEG       ! visible veg tot albedo
REAL, DIMENSION(:)  , INTENT(OUT), OPTIONAL  :: PALBNIR_TSOIL      ! nearIR  soil tot albedo
REAL, DIMENSION(:)  , INTENT(OUT), OPTIONAL  :: PALBVIS_TSOIL      ! visible soil tot albedo
END SUBROUTINE TEB_VEG_PROPERTIES
END INTERFACE
END MODULE MODI_TEB_VEG_PROPERTIES
