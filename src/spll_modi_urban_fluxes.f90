!auto_modi:spll_urban_fluxes.D
MODULE MODI_URBAN_FLUXES
INTERFACE
    SUBROUTINE URBAN_FLUXES(TOP, T, B, DMT, HIMPLICIT_WIND, PT_CANYON, PPEW_A_COEF, PPEW_B_COEF,      &
                            PEXNS, PRHOA, PVMOD, PH_TRAFFIC, PLE_TRAFFIC, PAC_WL, PCD, PDF_RF,        &
                            PDN_RF, PDF_RD, PDN_RD, PRNSN_RF, PHSN_RF, PLESN_RF, PGSN_RF,             &
                            PRNSN_RD, PHSN_RD, PLESN_RD, PGSN_RD, PMELT_RF, PDQS_RF, PMELT_RD,        &
                            PDQS_RD, PDQS_WL_A, PDQS_WL_B, PFLX_BLD_RF, PFLX_BLD_WL_A,                &
                            PFLX_BLD_WL_B, PFLX_BLD_FL, PFLX_BLD_MA, PE_SHADING, PLEW_RF,             &
                            PRN_GR, PH_GR, PLE_GR, PGFLUX_GR,                                         &
                            PLEW_RD, PLE_WL_A, PLE_WL_B, PMELT_BLT, PUSTAR_TWN                        )
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_t
IMPLICIT NONE
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(BEM_t), INTENT(INOUT) :: B
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: DMT
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
REAL, DIMENSION(:), INTENT(IN)    :: PT_CANYON    ! canyon air temperature
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_A_COEF  ! implicit coefficients
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_B_COEF  ! for wind coupling
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS        ! surface exner function
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA        ! air density
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD        ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PH_TRAFFIC    ! anthropogenic sensible
REAL, DIMENSION(:), INTENT(IN)    :: PLE_TRAFFIC   ! anthropogenic latent
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WL     ! surface conductance
REAL, DIMENSION(:), INTENT(IN)    :: PCD          ! drag coefficient
REAL, DIMENSION(:), INTENT(IN)    :: PDF_RF     ! snow-free    roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDN_RF     ! snow-covered roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDF_RD     ! snow-free    road fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDN_RD     ! snow-covered road fraction
REAL, DIMENSION(:), INTENT(IN)    :: PRNSN_RF ! net radiation over snow
REAL, DIMENSION(:), INTENT(IN)    :: PHSN_RF  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PLESN_RF ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PGSN_RF  ! flux under the snow
REAL, DIMENSION(:), INTENT(IN)    :: PRNSN_RD ! net radiation over snow
REAL, DIMENSION(:), INTENT(IN)    :: PHSN_RD  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PLESN_RD ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PGSN_RD  ! flux under the snow
REAL, DIMENSION(:), INTENT(IN)    :: PMELT_RF   ! snow melting on roof
REAL, DIMENSION(:), INTENT(IN)    :: PDQS_RF    ! storage inside roofs
REAL, DIMENSION(:), INTENT(IN)    :: PMELT_RD   ! snow melting on road
REAL, DIMENSION(:), INTENT(IN)    :: PDQS_RD    ! storage inside roads
REAL, DIMENSION(:), INTENT(IN)    :: PDQS_WL_A  ! storage inside walls
REAL, DIMENSION(:), INTENT(IN)    :: PDQS_WL_B  ! storage inside walls
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_RF! heat flx from bld to roof
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_WL_A! heat flx from bld to wall A
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_WL_B! heat flx from bld to wall B
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_FL! heat flx from bld to floor
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_MA! heat flx from bld to mass
REAL, DIMENSION(:), INTENT(IN)    :: PE_SHADING   ! energy not ref., nor absorbed, nor
REAL, DIMENSION(:), INTENT(IN)   :: PLEW_RF     ! latent heat flux over snow-free roof
REAL, DIMENSION(:), INTENT(IN)   :: PLEW_RD     ! latent heat flux of snow-free road
REAL, DIMENSION(:), INTENT(IN)    :: PRN_GR     ! net radiation over greenroof
REAL, DIMENSION(:), INTENT(IN)    :: PH_GR      ! sensible heat flux over greenroof
REAL, DIMENSION(:), INTENT(IN)    :: PLE_GR     ! latent heat flux over greenroof
REAL, DIMENSION(:), INTENT(IN)    :: PGFLUX_GR  ! flux through the greenroof
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WL_A   ! latent heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WL_B   ! latent heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_BLT    ! snow melting for town
REAL, DIMENSION(:), INTENT(OUT)   :: PUSTAR_TWN  ! friction velocity over town
END SUBROUTINE URBAN_FLUXES
END INTERFACE
END MODULE MODI_URBAN_FLUXES
