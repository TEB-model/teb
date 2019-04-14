!auto_modi:spll_floor_layer_e_budget.D
MODULE MODI_FLOOR_LAYER_E_BUDGET
INTERFACE
    SUBROUTINE FLOOR_LAYER_E_BUDGET(B, PTSTEP, PFLX_BLD_FL, PDQS_FL, PIMB_FL, PRADHT_IN,  &
                                    PRAD_WL_FL, PRAD_RF_FL, PRAD_WIN_FL, PLOAD_FL,        &
                                    PRAD_FL_MA, PCONV_FL_BLD                  )
USE MODD_BEM_n, ONLY : BEM_t
IMPLICIT NONE
TYPE(BEM_t), INTENT(INOUT) :: B
REAL,                 INTENT(IN)    :: PTSTEP       ! time step
REAL, DIMENSION(:),   INTENT(OUT)  :: PFLX_BLD_FL !flux from building to floor
REAL, DIMENSION(:),   INTENT(OUT) :: PDQS_FL !heat storage inside the floor
REAL, DIMENSION(:),   INTENT(OUT) :: PIMB_FL !floor energy residual imbalance for verification
REAL, DIMENSION(:), INTENT(IN)    :: PRADHT_IN      ! Indoor radiant heat transfer coefficient
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_RF_FL ! rad. fluxes from roof to floor[W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_WL_FL ! rad. fluxes from wall to floor[W m-2(wall)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_WIN_FL  ! rad. fluxes from win to floor[W m-2(win)]
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_FL ! solar and internal load to the floor
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_FL_MA  ! rad. fluxes from floor to mass [W m-2(floor)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_FL_BLD  ! conv. fluxes from floor to bld [W m-2(floor)]
END SUBROUTINE FLOOR_LAYER_E_BUDGET
END INTERFACE
END MODULE MODI_FLOOR_LAYER_E_BUDGET
