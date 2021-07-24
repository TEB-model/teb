!auto_modi:spll_isba.D
MODULE MODI_ISBA
INTERFACE
      SUBROUTINE ISBA(IO, KK, PK, PEK, G, AG, DK, DEK, DMK, TPTIME, PPOI, PABC, PIACAN, &
                      OMEB, PTSTEP, HIMPLICIT_WIND, PZREF, PUREF, PDIRCOSZW,     &
                      PTA, PQA, PEXNA, PRHOA, PPS, PEXNS, PRR, PSR, PZENITH,     &
                      PSCA_SW, PSW_RAD, PLW_RAD, PVMOD, PPEW_A_COEF, PPEW_B_COEF,&
                      PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,        &
                      PALBNIR_TVEG, PALBVIS_TVEG, PALBNIR_TSOIL, PALBVIS_TSOIL,  &
                      PPALPHAN, PZ0G_WITHOUT_SNOW, PZ0_MEBV, PZ0H_MEBV,          &
                      PZ0EFF_MEBV, PZ0_MEBN, PZ0H_MEBN, PZ0EFF_MEBN, PTDEEP_A,   &
                      PCSP, PFFG_NOSNOW, PFFV_NOSNOW, PEMIST, PUSTAR, PAC_AGG,   &
                      PHU_AGG, PRESP_BIOMASS_INST, PDEEP_FLUX, PIRRIG_GR     )
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_K_t, ISBA_P_t, ISBA_PE_t
USE MODD_SFX_GRID_n, ONLY : GRID_t
USE MODD_AGRI_n, ONLY : AGRI_t
USE MODD_DIAG_n, ONLY : DIAG_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
USE MODD_TYPE_DATE_SURF, ONLY : DATE_TIME
IMPLICIT NONE
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_K_t), INTENT(INOUT) :: KK
TYPE(ISBA_P_t), INTENT(INOUT) :: PK
TYPE(ISBA_PE_t), INTENT(INOUT) :: PEK
TYPE(GRID_t), INTENT(INOUT) :: G
TYPE(AGRI_t), INTENT(INOUT) :: AG
TYPE(DIAG_t), INTENT(INOUT) :: DK
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DEK
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DMK
TYPE(DATE_TIME), INTENT(IN)       :: TPTIME     ! current date and time
REAL, DIMENSION(:),    INTENT(IN) :: PPOI       ! Gaussian weights (as above)
REAL, DIMENSION(:), INTENT(INOUT) :: PABC       ! abscissa needed for integration
REAL, DIMENSION(:,:),   INTENT(OUT) :: PIACAN   ! PAR in the canopy at different gauss level
LOGICAL, INTENT(IN)                 :: OMEB     ! True = patch with multi-energy balance 
REAL,                 INTENT(IN)  :: PTSTEP     ! timestep of the integration
CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
REAL, DIMENSION(:),   INTENT(IN) :: PZREF       ! normal distance of the first
REAL, DIMENSION(:),   INTENT(IN) :: PUREF       ! reference height of the wind
REAL, DIMENSION(:),   INTENT(IN) ::  PDIRCOSZW  ! Director Cosinus along z
REAL, DIMENSION(:), INTENT(IN)  :: PTA        ! Temperature
REAL, DIMENSION(:), INTENT(IN)  :: PQA        ! specific humidity
REAL, DIMENSION(:), INTENT(IN)  :: PEXNA      ! Exner function
REAL, DIMENSION(:), INTENT(IN)  :: PRHOA      ! air density
REAL, DIMENSION(:), INTENT(IN)  :: PPS        ! Pressure
REAL, DIMENSION(:), INTENT(IN)  :: PEXNS      ! Exner function
REAL, DIMENSION(:), INTENT(IN)  :: PRR        ! Rain rate (in kg/m2/s)
REAL, DIMENSION(:), INTENT(IN)  :: PSR        ! Snow rate (in kg/m2/s)
REAL, DIMENSION(:), INTENT(IN)  :: PZENITH    ! solar zenith angle
REAL, DIMENSION(:), INTENT(IN)  :: PSW_RAD    ! solar   incoming radiation
REAL, DIMENSION(:), INTENT(IN)  :: PSCA_SW    ! solar diffuse incoming radiation
REAL, DIMENSION(:), INTENT(IN)  :: PLW_RAD    ! thermal incoming radiation
REAL, DIMENSION(:), INTENT(IN)  :: PVMOD      ! modulus of the wind
REAL, DIMENSION(:), INTENT(IN)  :: PPEW_A_COEF, PPEW_B_COEF, &
                                   PPET_A_COEF, PPEQ_A_COEF, &
                                   PPET_B_COEF, PPEQ_B_COEF  
REAL, DIMENSION(:), INTENT(IN)  :: PALBNIR_TVEG  ! tot albedo of vegetation in NIR (needed for LM_TR)
REAL, DIMENSION(:), INTENT(IN)  :: PALBVIS_TVEG  ! tot albedo of vegetation in VIS
REAL, DIMENSION(:), INTENT(IN)  :: PALBNIR_TSOIL ! tot albedo of bare soil in NIR
REAL, DIMENSION(:), INTENT(IN)  :: PALBVIS_TSOIL ! tot albedo of bare soil in VIS
REAL, DIMENSION(:), INTENT(IN)    :: PPALPHAN           ! snow/canopy transition coefficient
REAL, DIMENSION(:), INTENT(IN)    :: PZ0G_WITHOUT_SNOW  ! roughness length for momentum at snow-free canopy floor
REAL, DIMENSION(:), INTENT(IN)    :: PZ0_MEBV           ! roughness length for momentum over MEB vegetation part of patch
REAL, DIMENSION(:), INTENT(IN)    :: PZ0H_MEBV          ! roughness length for heat over MEB vegetation part of path
REAL, DIMENSION(:), INTENT(IN)    :: PZ0EFF_MEBV        ! roughness length for momentum over MEB vegetation part of patch
REAL, DIMENSION(:), INTENT(IN)    :: PZ0_MEBN           ! roughness length for momentum over MEB snow part of patch
REAL, DIMENSION(:), INTENT(IN)    :: PZ0H_MEBN          ! roughness length for heat over MEB snow part of path
REAL, DIMENSION(:), INTENT(IN)    :: PZ0EFF_MEBN        ! roughness length for momentum over MEB snow part of patch
REAL, DIMENSION(:), INTENT(IN)  :: PTDEEP_A       ! Deep soil temperature (prescribed)
REAL, DIMENSION(:),    INTENT(IN) :: PCSP       ! atmospheric CO2 concentration
REAL, DIMENSION(:), INTENT(IN)   :: PFFG_NOSNOW ! Without snow (ES)
REAL, DIMENSION(:), INTENT(IN)   :: PFFV_NOSNOW ! Without snow (ES)
REAL, DIMENSION(:), INTENT(OUT) :: PEMIST     ! grid-area surface emissivity
REAL, DIMENSION(:), INTENT(OUT) :: PUSTAR     ! friction velocity
REAL, DIMENSION(:),  INTENT(OUT) :: PAC_AGG  ! aggregated aerodynamic conductance
REAL, DIMENSION(:),  INTENT(OUT) :: PHU_AGG  ! aggregated relative humidity
REAL, DIMENSION(:,:),   INTENT(OUT) :: PRESP_BIOMASS_INST  ! instantaneous biomass respiration (kgCO2/kgair m/s)
REAL, DIMENSION(:),     INTENT(OUT) :: PDEEP_FLUX ! Heat flux at bottom of ISBA (W/m2)
REAL   ,DIMENSION(:),INTENT(IN)    :: PIRRIG_GR ! ground irrigation rate (kg/m2/s)
END SUBROUTINE ISBA
END INTERFACE
END MODULE MODI_ISBA
