!auto_modi:spll_isba_properties.D
MODULE MODI_ISBA_PROPERTIES
INTERFACE
      SUBROUTINE ISBA_PROPERTIES(IO, PEK, PDIR_SW, PSCA_SW, PSW_BANDS, KSW,  &
                                 PASNOW, PANOSNOW, PESNOW, PENOSNOW,         &
                                 PTSSNOW, PTSNOSNOW,                         &
                                 PALBNIR_TVEG, PALBVIS_TVEG, PALBNIR_TSOIL, PALBVIS_TSOIL )  
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_PE_t
USE MODD_TYPE_SNOW
IMPLICIT NONE
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_PE_t), INTENT(INOUT) :: PEK
REAL, DIMENSION(:,:), INTENT(IN)   :: PDIR_SW            ! direct incoming solar radiation
REAL, DIMENSION(:,:), INTENT(IN)   :: PSCA_SW            ! diffus incoming solar radiation
REAL, DIMENSION(:)  , INTENT(IN)   :: PSW_BANDS          ! mean wavelength of each shortwave band (m)
INTEGER,              INTENT(IN)   :: KSW                ! number of short-wave spectral bands            
REAL, DIMENSION(:)  , INTENT(OUT)  :: PASNOW    ! = snow albedo
REAL, DIMENSION(:)  , INTENT(OUT)  :: PANOSNOW  ! = snow free albedo 
REAL, DIMENSION(:)  , INTENT(OUT)  :: PESNOW    ! = snow emissivity
REAL, DIMENSION(:)  , INTENT(OUT)  :: PENOSNOW  ! = snow free emissivity
REAL, DIMENSION(:)  , INTENT(OUT)  :: PTSSNOW   ! = snow radiative temperature
REAL, DIMENSION(:)  , INTENT(OUT)  :: PTSNOSNOW ! = snow free radiative temperature
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBNIR_TVEG       ! nearIR  veg tot albedo
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBVIS_TVEG       ! visible veg tot albedo
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBNIR_TSOIL      ! nearIR  soil tot albedo
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBVIS_TSOIL      ! visible soil tot albedo
END SUBROUTINE ISBA_PROPERTIES
END INTERFACE
END MODULE MODI_ISBA_PROPERTIES
