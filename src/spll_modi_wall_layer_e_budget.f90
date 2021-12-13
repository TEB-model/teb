!auto_modi:spll_wall_layer_e_budget.D
MODULE MODI_WALL_LAYER_E_BUDGET
INTERFACE
    SUBROUTINE WALL_LAYER_E_BUDGET(TOP, T, B, PT_WL, PTS_WL_B, PTI_WL_B, PTSTEP, PDN_RD,     &
                                   PRHOA, PAC_WL, PAC_BLD, PLW_RAD, PPS, PEXNS, PABS_SW_WL,  &
                                   PT_CANYON, PTS_RD, PTSNOW_RD, PTS_GD, PTS_FL,             &
                                   PLW_WA_TO_WB, PLW_R_TO_W, PLW_G_TO_W, PLW_NR_TO_W,        &
                                   PLW_WIN_TO_W, PLW_S_TO_W, PFLX_BLD_WL, PDQS_WL,           &
                                   PABS_LW_WL, PEMIT_LW_WL, PH_WL, PIMB_WL, PRADHT_IN,       &
                                   PRAD_RF_WL, PRAD_WL_WIN, PRAD_WL_FL, PRAD_WL_MA,          &
                                   PCONV_WL_BLD, PLOAD_IN_WL ) 
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_BEM_n, ONLY : BEM_t
IMPLICIT NONE
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(BEM_t), INTENT(INOUT) :: B
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_WL    ! wall layers temperatures
REAL, DIMENSION(:),   INTENT(IN)  :: PTS_WL_B ! opposite wall surface temperature
REAL, DIMENSION(:),   INTENT(IN)  :: PTI_WL_B ! opposite wall internal temperature
REAL,               INTENT(IN)    :: PTSTEP       ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PDN_RD     ! snow-covered fraction on roads
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA        ! rho
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WL     ! aerodynamical conductance [m/s]
REAL, DIMENSION(:), INTENT(IN)    :: PAC_BLD      ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD      ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PPS          ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS        ! surface Exner function
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WL ! absorbed solar radiation
REAL, DIMENSION(:), INTENT(IN)    :: PT_CANYON    ! air canyon temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_RD     ! road surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTSNOW_RD  ! road snow temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GD   ! green area surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_FL    ! floor layers temperatures [K]
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_WB ! LW interactions wall  -> opposite wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_W   ! LW interactions road -> wall ; DEF formulation 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_W   ! LW interactions garden -> wall ; DEF formulation
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_W   ! LW interactions sky   -> wall 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_W  ! LW interactions road(snow) -> wall 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WIN_TO_W ! Radiative heat trasfer coeff wall-window 
REAL, DIMENSION(:), INTENT(OUT)   :: PFLX_BLD_WL! flux from bld to wall
REAL, DIMENSION(:), INTENT(INOUT) :: PDQS_WL    ! heat storage inside the wall 
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_WL ! absorbed infrared rad. [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PEMIT_LW_WL  ! LW flux emitted by the wall [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PH_WL      ! Sensible heat flux from wall to air [W/m2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PIMB_WL    ! wall residual energy imbalance 
REAL, DIMENSION(:),   INTENT(IN)  :: PRADHT_IN     ! Indoor radiant heat transfer coefficient
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_RF_WL ! rad. fluxes from roof to wall [W m-2(roof)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WL_WIN  ! rad. fluxes from wall to win  [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WL_FL! rad. fluxes from wall to floor [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WL_MA ! rad. fluxes from wall to mass [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_WL_BLD ! conv. fluxes from wall to bld [W m-2(wall)]
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_IN_WL  ! LOAD from solar heat gain + rad int. gains  W/m2 [Wall]
END SUBROUTINE WALL_LAYER_E_BUDGET
END INTERFACE
END MODULE MODI_WALL_LAYER_E_BUDGET
