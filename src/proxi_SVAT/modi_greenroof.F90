!auto_modi:spll_greenroof.D
MODULE MODI_GREENROOF
INTERFACE
    SUBROUTINE GREENROOF (DTCO, G, T, TOP, TIR, DTV, GB, DK, DEK, DMK, GRO, S, K, P, PEK,    &
                          HIMPLICIT_WIND, TPTIME, PTSUN, PPEW_A_COEF, PPEW_B_COEF,  &
                          PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,       &
                          PTSTEP, PZREF, PUREF, PALB_GD, PTA, PQA, PEXNS, PEXNA, PRHOA,      &
                          PCO2, PPS, PRR, PSR, PZENITH, PSW, PLW, PVMOD,            &
                          PALBNIR_TVEG, PALBVIS_TVEG, PALBNIR_TSOIL, PALBVIS_TSOIL, &
                          PRN, PH, PLE, PGFLUX, PSFCO2, PEVAP, PUW, PRUNOFF, PDRAIN,&
                          PAC, PQSAT, PTSRAD, PAC_AGG, PHU_AGG, PDEEP_FLUX, PIRRIG )  
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_SFX_GRID_n, ONLY : GRID_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_IRRIG_n, ONLY : TEB_IRRIG_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_t
USE MODD_DIAG_n, ONLY : DIAG_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_S_t, ISBA_K_t, ISBA_P_t, ISBA_PE_t
USE MODD_TYPE_DATE_SURF,       ONLY: DATE_TIME
IMPLICIT NONE
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(GRID_t), INTENT(INOUT) :: G
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_IRRIG_t), INTENT(INOUT) :: TIR
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTV
TYPE(GR_BIOG_t), INTENT(INOUT) :: GB
TYPE(DIAG_t), INTENT(INOUT) :: DK
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DEK
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DMK
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: GRO
TYPE(ISBA_S_t), INTENT(INOUT) :: S
TYPE(ISBA_K_t), INTENT(INOUT) :: K
TYPE(ISBA_P_t), INTENT(INOUT) :: P
TYPE(ISBA_PE_t), INTENT(INOUT) :: PEK
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
TYPE(DATE_TIME)     , INTENT(IN)    :: TPTIME             ! current date and time from teb
REAL, DIMENSION(:)  , INTENT(IN)    :: PTSUN              ! solar time      (s from midnight)
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEW_A_COEF        ! implicit coefficients
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEW_B_COEF        ! for wind coupling
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEQ_A_COEF        ! implicit coefficients
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEQ_B_COEF        ! for humidity
REAL, DIMENSION(:)  , INTENT(IN)    :: PPET_A_COEF        ! implicit coefficients
REAL, DIMENSION(:)  , INTENT(IN)    :: PPET_B_COEF        ! for temperature
REAL                , INTENT(IN)    :: PTSTEP             ! time step
REAL, DIMENSION(:)  , INTENT(IN)    :: PZREF              ! height of the first atmospheric level
REAL, DIMENSION(:)  , INTENT(IN)    :: PUREF              ! reference height for the wind
REAL, DIMENSION(:)  , INTENT(IN)    :: PALB_GD            ! green areas albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PTA                ! temperature at first atm. level 
REAL, DIMENSION(:)  , INTENT(IN)    :: PQA                ! specific humidity at first atm. level
REAL, DIMENSION(:)  , INTENT(IN)    :: PPS                ! pressure at the surface
REAL, DIMENSION(:)  , INTENT(IN)    :: PEXNA              ! Exner function at first atm. level
REAL, DIMENSION(:)  , INTENT(IN)    :: PEXNS              ! surface Exner function
REAL, DIMENSION(:)  , INTENT(IN)    :: PRHOA              ! air density at the lowest level
REAL, DIMENSION(:)  , INTENT(IN)    :: PCO2               ! CO2 concentration in the air    (kg/m3)
REAL, DIMENSION(:)  , INTENT(IN)    :: PRR                ! rain rate
REAL, DIMENSION(:)  , INTENT(IN)    :: PSR                ! snow rate
REAL, DIMENSION(:)  , INTENT(IN)    :: PZENITH            ! solar zenithal angle
REAL, DIMENSION(:)  , INTENT(IN)    :: PSW                ! incoming total solar rad on an horizontal surface
REAL, DIMENSION(:)  , INTENT(IN)    :: PLW                ! atmospheric infrared radiation
REAL, DIMENSION(:)  , INTENT(IN)    :: PVMOD              ! module of horizontal wind near first atm. level
REAL, DIMENSION(:)  , INTENT(IN)    :: PALBNIR_TVEG       ! nearIR  veg tot albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PALBVIS_TVEG       ! visible veg tot albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PALBNIR_TSOIL      ! nearIR  soil tot albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PALBVIS_TSOIL      ! visible soil tot albedo
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN         ! net radiation over greenroofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH          ! sensible heat flux over greenroofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE         ! latent heat flux over greenroofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX      ! flux through the greenroofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PSFCO2      ! flux of greenroof CO2       (m/s*kg_CO2/kg_air)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PEVAP       ! total evaporation over greenroofs (kg/m2/s)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PUW         ! friction flux (m2/s2)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRUNOFF     ! greenroof surface runoff
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDRAIN      ! greenroof surface drainage
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC         ! greenroof aerodynamical conductance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PQSAT       ! saturation humidity
REAL, DIMENSION(:)  , INTENT(OUT)   :: PTSRAD      ! greenroof radiative surface temp. (snow free)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_AGG     ! aggreg. aeodynamic resistance for greenroofs for latent heat flux
REAL, DIMENSION(:)  , INTENT(OUT)   :: PHU_AGG     ! aggreg. relative humidity for greenroofs for latent heat flux
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDEEP_FLUX  ! Heat Flux at the bottom layer of the greenroof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PIRRIG      ! greenroof summer irrigation rate
END SUBROUTINE GREENROOF
END INTERFACE
END MODULE MODI_GREENROOF
