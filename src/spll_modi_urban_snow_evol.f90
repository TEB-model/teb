!auto_modi:spll_urban_snow_evol.D
MODULE MODI_URBAN_SNOW_EVOL
INTERFACE
    SUBROUTINE URBAN_SNOW_EVOL(T, B, PT_LWCN, PQ_LWCN, PU_LWCN, PTS_RF, PTS_RD, PTS_WL_A,   &
                               PTS_WL_B, PPS, PTA, PQA, PRHOA, PLW_RAD, PSR, PZREF, PUREF,  &
                               PVMOD, PTSTEP, PZ_LWCN, PDN_RF, PABS_SW_SN_RF, PABS_LW_SN_RF,&
                               PDN_RD, PABS_SW_SN_RD, PABS_LW_SN_RD, PRNSN_RF, PHSN_RF,     &
                               PLESN_RF, PGSN_RF, PMELT_RF, PRNSN_RD, PHSN_RD, PLESN_RD,    &
                               PGSN_RD, PMELT_RD, PLW_WA_TO_NR , PLW_WB_TO_NR, PLW_S_TO_NR, &
                               PLW_WIN_TO_NR, PDQS_SN_RF, PDQS_SN_RD   )  
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_BEM_n, ONLY : BEM_t
IMPLICIT NONE
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(BEM_t), INTENT(INOUT) :: B
REAL, DIMENSION(:),   INTENT(IN)    :: PT_LWCN  ! LWCN air temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PQ_LWCN  ! LWCN air specific humidity
REAL, DIMENSION(:),   INTENT(IN)    :: PU_LWCN  ! LWCN hor. wind
REAL, DIMENSION(:),   INTENT(IN)    :: PTS_RF   ! roof surface temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PTS_RD   ! road surface temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PTS_WL_A ! wall surface temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PTS_WL_B ! wall surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PPS      ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PTA      ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA      ! specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD    ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA    ! air density at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD  ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PSR      ! snow rate
REAL, DIMENSION(:), INTENT(IN)    :: PZREF    ! reference height of the first
REAL, DIMENSION(:), INTENT(IN)    :: PUREF    ! reference height of the first
REAL,               INTENT(IN)    :: PTSTEP   ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PZ_LWCN  ! height of forcing
REAL, DIMENSION(:), INTENT(IN)    :: PDN_RF          ! snow-covered roof frac.
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_SN_RF ! SW absorbed by roof snow
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_SN_RF ! absorbed IR rad by snow on roof
REAL, DIMENSION(:), INTENT(INOUT) :: PDN_RD          ! snow-covered road frac.
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_SN_RD ! SW absorbed by road snow
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_SN_RD ! absorbed IR rad by snow on road
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
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_NR        ! LW contrib. wall       -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WB_TO_NR        ! LW contrib. wall       -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_NR         ! LW contrib. sky        -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WIN_TO_NR       ! LW contrib. win       -> road(snow)
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_SN_RF ! Heat storage in snowpack on roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_SN_RD ! Heat storage in snowpack on roads
END SUBROUTINE URBAN_SNOW_EVOL
END INTERFACE
END MODULE MODI_URBAN_SNOW_EVOL
