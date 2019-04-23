!auto_modi:spll_teb_garden.D
MODULE MODI_TEB_GARDEN
INTERFACE
    SUBROUTINE TEB_GARDEN (TOP, T, BOP, B, TPN, TIR, DMT,                                         &
                           HIMPLICIT_WIND, PTSUN, PT_CAN, PQ_CAN, PU_CAN, PT_LOWCAN, PQ_LOWCAN,   &
                           PU_LOWCAN, PZ_LOWCAN, PPEW_A_COEF, PPEW_B_COEF, PPEW_A_COEF_LOWCAN,    &
                           PPEW_B_COEF_LOWCAN, PPS, PPA, PEXNS, PEXNA, PTA, PQA, PRHOA, PCO2,     &
                           PLW_RAD, PDIR_SW, PSCA_SW, PSW_BANDS, KSW, PZENITH, PAZIM, PRR, PSR,   &
                           PZREF, PUREF, PVMOD, PH_TRAFFIC, PLE_TRAFFIC, PTSTEP, PLEW_RF, PLEW_RD,&
                           PLE_WL_A, PLE_WL_B, PRNSN_RF, PHSN_RF, PLESN_RF, PGSN_RF, PMELT_RF,    &
                           PRNSN_RD, PHSN_RD, PLESN_RD, PGSN_RD, PMELT_RD, PRN_GRND, PH_GRND,     &
                           PLE_GRND, PGFLX_GRND, PRN_TWN, PH_TWN, PLE_TWN, PGFLX_TWN, PEVAP_TWN,  &
                           PSFCO2, PUW_GRND, PUW_RF, PDUWDU_GRND, PDUWDU_RF,                      &
                           PUSTAR_TWN, PCD, PCDN, PCH_TWN, PRI_TWN, PTS_TWN, PEMIS_TWN,           &
                           PDIR_ALB_TWN, PSCA_ALB_TWN, PRESA_TWN, PAC_RD, PAC_GD, PAC_GR,         &
                           PAC_RD_WAT, PAC_GD_WAT, PAC_GR_WAT, KDAY, PEMIT_LW_FAC, PEMIT_LW_GRND, &
                           PT_RAD_IND, PREF_SW_GRND, PREF_SW_FAC, PHU_BLD, PTIME, PPROD_BLD  )
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_TEB_PANEL_n, ONLY : TEB_PANEL_t
USE MODD_TEB_IRRIG_n, ONLY : TEB_IRRIG_t
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_t
IMPLICIT NONE
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(BEM_OPTIONS_t), INTENT(INOUT) :: BOP
TYPE(BEM_t), INTENT(INOUT) :: B
TYPE(TEB_PANEL_t), INTENT(INOUT) :: TPN
TYPE(TEB_IRRIG_t), INTENT(INOUT) :: TIR
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: DMT
!
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND      ! wind implicitation option
!                                                         ! 'OLD' = direct
!                                                         ! 'NEW' = Taylor serie, order 1
REAL, DIMENSION(:),   INTENT(IN)    :: PTSUN              ! solar time   (s from midnight)
!
REAL, DIMENSION(:)  , INTENT(INOUT) :: PT_CAN             ! canyon air temperature
REAL, DIMENSION(:)  , INTENT(INOUT) :: PQ_CAN             ! canyon air specific humidity
REAL, DIMENSION(:)  , INTENT(IN)    :: PU_CAN             ! canyon hor. wind
REAL, DIMENSION(:)  , INTENT(IN)    :: PU_LOWCAN          ! wind near the road
REAL, DIMENSION(:)  , INTENT(IN)    :: PT_LOWCAN          ! temp. near the road
REAL, DIMENSION(:)  , INTENT(IN)    :: PQ_LOWCAN          ! hum. near the road
REAL, DIMENSION(:)  , INTENT(IN)    :: PZ_LOWCAN          ! height of atm. var. near the road
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEW_A_COEF        ! implicit coefficients
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEW_B_COEF        ! for wind coupling
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEW_A_COEF_LOWCAN ! implicit coefficients for wind coupling
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEW_B_COEF_LOWCAN ! between low canyon wind and road
REAL, DIMENSION(:)  , INTENT(IN)    :: PPS                ! pressure at the surface
REAL, DIMENSION(:)  , INTENT(IN)    :: PPA                ! pressure at the first atmospheric level
REAL, DIMENSION(:)  , INTENT(IN)    :: PEXNS              ! surface exner function
REAL, DIMENSION(:)  , INTENT(IN)    :: PTA                ! temperature at the lowest level
REAL, DIMENSION(:)  , INTENT(IN)    :: PQA                ! specific humidity at the lowest level
REAL, DIMENSION(:)  , INTENT(IN)    :: PVMOD              ! module of the horizontal wind
REAL, DIMENSION(:)  , INTENT(IN)    :: PEXNA              ! exner function at the lowest level
REAL, DIMENSION(:)  , INTENT(IN)    :: PRHOA              ! air density at the lowest level
REAL, DIMENSION(:)  , INTENT(IN)    :: PCO2               ! CO2 concentration in the air    (kg/m3)
REAL, DIMENSION(:)  , INTENT(IN)    :: PLW_RAD            ! atmospheric infrared radiation
REAL, DIMENSION(:,:), INTENT(IN)    :: PDIR_SW            ! incoming direct solar rad on an horizontal surface
REAL, DIMENSION(:,:), INTENT(IN)    :: PSCA_SW            ! scattered incoming solar rad.
REAL, DIMENSION(:)  , INTENT(IN)    :: PSW_BANDS          ! mean wavelength of each shortwave band (m)
INTEGER,              INTENT(IN)    :: KSW                ! number of short-wave spectral bands
REAL, DIMENSION(:)  , INTENT(IN)    :: PZENITH            ! solar zenithal angle
REAL, DIMENSION(:)  , INTENT(IN)    :: PAZIM              ! solar azimuthal angle
REAL, DIMENSION(:)  , INTENT(IN)    :: PRR                ! rain rate
REAL, DIMENSION(:)  , INTENT(IN)    :: PSR                ! snow rate
REAL, DIMENSION(:)  , INTENT(IN)    :: PH_TRAFFIC         ! anthropogenic sensible heat fluxes due to traffic
REAL, DIMENSION(:)  , INTENT(IN)    :: PLE_TRAFFIC        ! anthropogenic latent heat fluxes due to traffic
REAL, DIMENSION(:)  , INTENT(IN)    :: PZREF              ! reference height of the first atm level (temperature)
REAL, DIMENSION(:)  , INTENT(IN)    :: PUREF              ! reference height of the first atm level (wind)
REAL                , INTENT(IN)    :: PTSTEP             ! time step
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLEW_RF          ! latent heat flux over roof (snow)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLEW_RD          ! latent heat flux over road (snow)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_WL_A         ! latent heat flux over wall
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_WL_B         ! latent heat flux over wall
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRNSN_RF       ! net radiation over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PHSN_RF        ! sensible heat flux over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLESN_RF       ! latent heat flux over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGSN_RF        ! flux under the snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PMELT_RF         ! snow melt
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRNSN_RD       ! net radiation over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PHSN_RD        ! sensible heat flux over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLESN_RD       ! latent heat flux over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGSN_RD        ! flux under the snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PMELT_RD       ! snow melt
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_GRND           ! net radiation over ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_GRND            ! sensible heat flux over ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_GRND           ! latent heat flux over ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLX_GRND        ! flux through the ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_TWN           ! net radiation over town
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_TWN            ! sensible heat flux over town
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_TWN           ! latent heat flux over town
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLX_TWN         ! flux through the ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PEVAP_TWN         ! evaporation flux (kg/m2/s)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PSFCO2            ! flux of CO2       (m/s*kg_CO2/kg_air)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PUW_GRND          ! momentum flux for ground built surf
REAL, DIMENSION(:)  , INTENT(OUT)   :: PUW_RF            ! momentum flux for roofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDUWDU_GRND       !
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDUWDU_RF         !
REAL, DIMENSION(:)  , INTENT(OUT)   :: PUSTAR_TWN        ! friciton velocity over town
REAL, DIMENSION(:)  , INTENT(OUT)   :: PCD               ! town averaged drag coefficient
REAL, DIMENSION(:)  , INTENT(OUT)   :: PCDN              ! town averaged neutral drag coefficient
REAL, DIMENSION(:)  , INTENT(OUT)   :: PCH_TWN           ! town averaged heat transfer coefficient
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRI_TWN           ! town averaged Richardson number
REAL, DIMENSION(:)  , INTENT(OUT)   :: PTS_TWN           ! town surface temperature
REAL, DIMENSION(:)  , INTENT(OUT)   :: PEMIS_TWN         ! town equivalent emissivity
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDIR_ALB_TWN      ! town equivalent direct albedo
REAL, DIMENSION(:)  , INTENT(OUT)   :: PSCA_ALB_TWN      ! town equivalent diffuse albedo
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRESA_TWN         ! town aerodynamical resistance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_RD            ! road conductance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_GD            ! green area conductance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_GR            ! green roof conductance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_RD_WAT        ! road conductance for latent heat
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_GD_WAT        ! green area conductance for latent heat
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_GR_WAT        ! green roof conductance for latent heat
INTEGER             , INTENT(IN)    :: KDAY               ! Simulation day
REAL, DIMENSION(:)  , INTENT(OUT)    :: PEMIT_LW_GRND     ! LW flux emitted by the ground (W/m2 ground)
REAL, DIMENSION(:)  , INTENT(OUT)    :: PEMIT_LW_FAC      ! LW flux emitted by the facade (W/m2 ground)
REAL, DIMENSION(:)  , INTENT(OUT)    :: PT_RAD_IND        ! Indoor mean radiant temperature [K]
REAL, DIMENSION(:)  , INTENT(OUT)    :: PREF_SW_GRND      ! total solar rad reflected from ground
REAL, DIMENSION(:)  , INTENT(OUT)    :: PREF_SW_FAC       ! total solar rad reflected from facade
REAL, DIMENSION(:)  , INTENT(OUT)    :: PHU_BLD           ! Indoor relative humidity 0 < (-) < 1
REAL                , INTENT(IN)     :: PTIME             ! current time since midnight (UTC, s)
REAL, DIMENSION(:),   INTENT(OUT)     :: PPROD_BLD        ! Averaged     Energy production of solar panel on roofs (W/m2 bld  )
END SUBROUTINE TEB_GARDEN
END INTERFACE
END MODULE MODI_TEB_GARDEN
