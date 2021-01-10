!auto_modi:spll_teb.D
MODULE MODI_TEB
INTERFACE
    SUBROUTINE TEB  (TOP, T, BOP, B, TIR, DMT, HIMPLICIT_WIND, PTSUN,                   &
                     PT_CANYON, PQ_CANYON, PU_CANYON, PT_LOWCAN, PQ_LOWCAN, PU_LOWCAN,  &
                     PZ_LOWCAN, PPEW_A_COEF, PPEW_B_COEF, PPEW_A_COEF_LOWCAN,           &
                     PPEW_B_COEF_LOWCAN, PPS, PPA, PEXNS, PEXNA, PTA, PQA, PRHOA,       &
                     PLW_RAD, PRR, PSR, PZREF, PUREF, PVMOD, PH_TRAFFIC, PLE_TRAFFIC,   &
                     PTSTEP, PDF_RF, PDN_RF, PDF_RD, PDN_RD, PQSAT_RF, PQSAT_RD,        &
                     PDELT_RF, PDELT_RD, PTS_GARDEN, PLEW_RF, PUW_GR, PLEW_RD, PLE_WL_A,&
                     PLE_WL_B, PRNSN_RF, PHSN_RF, PLESN_RF, PGSN_RF, PMELT_RF, PRN_GR, &
                     PH_GR, PLE_GR, PGFLUX_GR, PDRAIN_GR, PRUNOFF_GR, PRNSN_RD,    &
                     PHSN_RD, PLESN_RD, PGSN_RD, PMELT_RD, PUW_RD, PUW_RF, PDUWDU_RD,   &
                     PDUWDU_RF, PUSTAR_TWN, PCD, PCDN, PCH_TWN, PRI_TWN, PRESA_TWN,     &
                     PAC_RF, PAC_RD, PAC_WL, PAC_TOP, PAC_GARDEN, PAC_RF_WAT,           &
                     PAC_RD_WAT, PLW_WA_TO_WB, PLW_WA_TO_R, PLW_WB_TO_R, PLW_WA_TO_NR,  &
                     PLW_WB_TO_NR, PLW_R_TO_WA, PLW_R_TO_WB, PLW_G_TO_WA, PLW_G_TO_WB,  &
                     PLW_S_TO_WA, PLW_S_TO_WB, PLW_S_TO_R, PLW_S_TO_NR, PLW_NR_TO_WA,   &
                     PLW_NR_TO_WB, PLW_NR_TO_WIN, PLW_WA_TO_WIN, PLW_WB_TO_WIN,         &
                     PLW_G_TO_WIN, PLW_R_TO_WIN, PLW_S_TO_WIN, PLW_WIN_TO_WA,           &
                     PLW_WIN_TO_WB, PLW_WIN_TO_R, PLW_WIN_TO_NR, KDAY, PEMIT_LW_FAC,    &
                     PEMIT_LW_RD, PT_RAD_IND, PHU_BLD, PTIME, PE_SHADING )
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_TEB_IRRIG_n, ONLY : TEB_IRRIG_t
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_t
IMPLICIT NONE
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(BEM_OPTIONS_t), INTENT(INOUT) :: BOP
TYPE(BEM_t), INTENT(INOUT) :: B
TYPE(TEB_IRRIG_t), INTENT(INOUT) :: TIR
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: DMT
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
REAL, DIMENSION(:),   INTENT(IN)    :: PTSUN              ! solar time   (s from midnight)
REAL, DIMENSION(:), INTENT(INOUT) :: PT_CANYON     ! canyon air temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PQ_CANYON     ! canyon air specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PU_CANYON     ! canyon hor. wind
REAL, DIMENSION(:), INTENT(IN)    :: PU_LOWCAN     ! wind near the road
REAL, DIMENSION(:), INTENT(IN)    :: PT_LOWCAN     ! temp. near the road
REAL, DIMENSION(:), INTENT(IN)    :: PQ_LOWCAN     ! hum. near the road
REAL, DIMENSION(:), INTENT(IN)    :: PZ_LOWCAN     ! height of atm. var. near the road
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_A_COEF   ! implicit coefficients
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_B_COEF   ! for wind coupling
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_A_COEF_LOWCAN ! implicit coefficients for wind coupling
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_B_COEF_LOWCAN ! between low canyon wind and road
REAL, DIMENSION(:), INTENT(IN)    :: PPS           ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PPA           ! pressure at the first atmospheric level
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS         ! surface exner function
REAL, DIMENSION(:), INTENT(IN)    :: PTA           ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA           ! specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD         ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PH_TRAFFIC    ! anthropogenic sensible
REAL, DIMENSION(:), INTENT(IN)    :: PLE_TRAFFIC   ! anthropogenic latent
REAL, DIMENSION(:), INTENT(IN)    :: PEXNA         ! exner function
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA         ! air density
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD       ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PRR           ! rain rate
REAL, DIMENSION(:), INTENT(IN)    :: PSR           ! snow rate
REAL, DIMENSION(:), INTENT(IN)    :: PZREF         ! reference height of the first
REAL, DIMENSION(:), INTENT(IN)    :: PUREF         ! reference height of the first
REAL,               INTENT(IN)    :: PTSTEP        ! time step
REAL, DIMENSION(:), INTENT(INOUT) :: PDF_RF      ! snow-free    fraction on roofs
REAL, DIMENSION(:), INTENT(INOUT) :: PDN_RF      ! snow-covered fraction on roofs
REAL, DIMENSION(:), INTENT(INOUT) :: PDF_RD      ! snow-free    fraction on roads
REAL, DIMENSION(:), INTENT(INOUT) :: PDN_RD      ! snow-covered fraction on roads
REAL, DIMENSION(:), INTENT(OUT)   :: PQSAT_RF    ! hum at saturation over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PQSAT_RD    ! hum at saturation over road
REAL, DIMENSION(:), INTENT(OUT)   :: PDELT_RF    ! water fraction on roof
REAL, DIMENSION(:), INTENT(OUT)   :: PDELT_RD    ! water fraction on road
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GARDEN    ! GARDEN area surf temp.
REAL, DIMENSION(:), INTENT(OUT)   :: PLEW_RF    ! latent heat flux over roof (snow)
REAL, DIMENSION(:), INTENT(OUT)   :: PLEW_RD    ! latent heat flux over road (snow)
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WL_A   ! latent heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WL_B   ! latent heat flux over wall
REAL, DIMENSION(:), INTENT(IN)    :: PUW_GR     ! Momentum flux for greenroofs
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSN_RF ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSN_RF  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESN_RF ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSN_RF  ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_RF   ! snow melt
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSN_RD ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSN_RD  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESN_RD ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSN_RD  ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_RD   ! snow melt
REAL, DIMENSION(:), INTENT(OUT)   :: PUW_RD     ! Momentum flux for roads
REAL, DIMENSION(:), INTENT(OUT)   :: PUW_RF     ! Momentum flux for roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PDUWDU_RD  !
REAL, DIMENSION(:), INTENT(OUT)   :: PDUWDU_RF  !
REAL, DIMENSION(:), INTENT(OUT)   :: PUSTAR_TWN ! friciton velocity over town
REAL, DIMENSION(:), INTENT(IN)    :: PRN_GR     ! net radiation over greenroof
REAL, DIMENSION(:), INTENT(IN)    :: PH_GR      ! sensible heat flux over greenroof
REAL, DIMENSION(:), INTENT(IN)    :: PLE_GR     ! latent heat flux over greenroof
REAL, DIMENSION(:), INTENT(IN)    :: PGFLUX_GR  ! flux through the greenroof
REAL, DIMENSION(:), INTENT(IN)    :: PRUNOFF_GR ! runoff over green roofs
REAL, DIMENSION(:), INTENT(IN)    :: PDRAIN_GR  ! outlet drainage at base of green roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PCD          ! town averaged drag coefficient
REAL, DIMENSION(:), INTENT(OUT)   :: PCDN         ! town averaged neutral drag coefficient
REAL, DIMENSION(:), INTENT(OUT)   :: PCH_TWN     ! town averaged heat transfer
REAL, DIMENSION(:), INTENT(OUT)   :: PRI_TWN      ! town averaged Richardson number
REAL, DIMENSION(:), INTENT(OUT)   :: PRESA_TWN    ! town aerodynamical resistance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_RF      ! roof conductance
REAL, DIMENSION(:), INTENT(INOUT) :: PAC_RD      ! road conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_WL      ! wall conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_TOP       ! top conductance
REAL, DIMENSION(:), INTENT(IN)    :: PAC_GARDEN    ! garden conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_RF_WAT  ! roof water conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_RD_WAT  ! roof water conductance
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_WB      ! LW contrib. wall A (orB) -> wall B (or A)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_R         ! LW contrib. wall       -> road
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WB_TO_R         ! LW contrib. wall       -> road
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_NR        ! LW contrib. wall       -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WB_TO_NR        ! LW contrib. wall       -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_WA         ! LW contrib. road       -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_WB         ! LW contrib. road       -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_WA         ! LW contrib. GARDEN     -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_WB         ! LW contrib. GARDEN     -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_WA        ! LW contrib. road(snow) -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_WB        ! LW contrib. road(snow) -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_WA         ! LW contrib. sky        -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_WB         ! LW contrib. sky        -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_R          ! LW contrib. sky        -> road
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_NR         ! LW contrib. sky        -> road(snow)
INTEGER,            INTENT(IN)     :: KDAY         ! Simulation day
REAL, DIMENSION(:), INTENT(IN)   :: PLW_WA_TO_WIN ! Radiative heat trasfer coeff wall-window 
REAL, DIMENSION(:), INTENT(IN)   :: PLW_WB_TO_WIN ! Radiative heat trasfer coeff wall-window 
REAL, DIMENSION(:), INTENT(IN)   :: PLW_G_TO_WIN  ! Radiative heat trasfer coeff garden-window 
REAL, DIMENSION(:), INTENT(IN)   :: PLW_R_TO_WIN  ! Radiative heat trasfer coeff road-window 
REAL, DIMENSION(:), INTENT(IN)   :: PLW_S_TO_WIN ! Radiative heat trasfer coeff window-sky 
REAL, DIMENSION(:), INTENT(IN)   :: PLW_WIN_TO_WA! Radiative heat trasfer coeff window-wall
REAL, DIMENSION(:), INTENT(IN)   :: PLW_WIN_TO_WB! Radiative heat trasfer coeff window-wall
REAL, DIMENSION(:), INTENT(IN)   :: PLW_WIN_TO_R ! Radiative heat trasfer coeff window-road 
REAL, DIMENSION(:), INTENT(IN)   :: PLW_NR_TO_WIN! Radiative heat trasfer coeff road(snow)-win 
REAL, DIMENSION(:), INTENT(IN)   :: PLW_WIN_TO_NR! Radiative heat trasfer coeff win-road(snow) 
REAL, DIMENSION(:), INTENT(OUT) :: PEMIT_LW_RD ! LW fluxes emitted by road (W/m2 surf road)
REAL, DIMENSION(:), INTENT(OUT) :: PEMIT_LW_FAC  ! LW fluxes emitted by wall (W/m2 surf wall)
REAL, DIMENSION(:), INTENT(OUT) :: PT_RAD_IND    ! Indoor mean radiant temperature [K]
REAL, DIMENSION(:), INTENT(OUT) :: PHU_BLD       ! Indoor relative humidity 0 < (-) < 1
REAL,                INTENT(IN)  :: PTIME        ! current time since midnight (UTC, s)
REAL, DIMENSION(:), INTENT(IN)  :: PE_SHADING    !energy not ref., nor absorbed, nor
END SUBROUTINE TEB
END INTERFACE
END MODULE MODI_TEB
