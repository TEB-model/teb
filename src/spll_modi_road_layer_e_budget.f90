!auto_modi:spll_road_layer_e_budget.D
MODULE MODI_ROAD_LAYER_E_BUDGET
INTERFACE
    SUBROUTINE ROAD_LAYER_E_BUDGET(T, B, PTSTEP, PDN_ROAD, PRHOA, PAC_ROAD, PAC_ROAD_WAT, &
                                   PLW_RAD, PPS, PQSAT_ROAD, PDELT_ROAD, PEXNS,           &
                                   PABS_SW_ROAD, PGSNOW_ROAD, PQ_LOWCAN, PT_LOWCAN,       &
                                   PTS_WALL_A, PTS_WALL_B, PTSNOW_ROAD, PTS_GARDEN,       &
                                   PLW_WA_TO_R, PLW_WB_TO_R, PLW_S_TO_R, PLW_WIN_TO_R,    &
                                   PEMIT_LW_ROAD, PDQS_ROAD, PABS_LW_ROAD, PHFREE_ROAD,   &
                                   PLEFREE_ROAD, PIMB_ROAD, PRR )
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_BEM_n, ONLY : BEM_t
IMPLICIT NONE
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(BEM_t), INTENT(INOUT) :: B
REAL,               INTENT(IN)    :: PTSTEP       ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROAD     ! road snow fraction
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA        ! rho
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROAD     ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROAD_WAT ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD      ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PPS          ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PQSAT_ROAD   ! q_sat(Ts)
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_ROAD   ! fraction of water
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS        ! surface Exner function
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROAD ! absorbed solar radiation
REAL, DIMENSION(:), INTENT(IN)    :: PGSNOW_ROAD  ! road snow conduction
REAL, DIMENSION(:), INTENT(IN)    :: PQ_LOWCAN    ! and specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PT_LOWCAN    ! low canyon air temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_WALL_A   ! wall surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_WALL_B   ! wall surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTSNOW_ROAD  ! road snow temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GARDEN   ! green area surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_R   ! LW interactions wall  -> road 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WB_TO_R   ! LW interactions wall  -> road 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_R    ! LW interactions sky   -> road 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WIN_TO_R ! LW interactions window -> road 
REAL, DIMENSION(:), INTENT(OUT)   :: PEMIT_LW_ROAD! LW flux emitted by the road (W/m2 of road)
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_ROAD    !heat storage inside the road
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_ROAD ! absorbed infrared rad.
REAL, DIMENSION(:), INTENT(OUT)   :: PHFREE_ROAD  ! sensible heat flux on the
REAL, DIMENSION(:), INTENT(OUT)   :: PLEFREE_ROAD ! latent heat flux on the
REAL, DIMENSION(:), INTENT(OUT)   :: PIMB_ROAD    ! road residual energy imbalance 
REAL, DIMENSION(:), INTENT(IN)    :: PRR          ! rain rate [kg m-2 s-1]
END SUBROUTINE ROAD_LAYER_E_BUDGET
END INTERFACE
END MODULE MODI_ROAD_LAYER_E_BUDGET
