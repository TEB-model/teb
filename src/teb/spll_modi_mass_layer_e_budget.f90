!auto_modi:spll_mass_layer_e_budget.D
MODULE MODI_MASS_LAYER_E_BUDGET
INTERFACE
    SUBROUTINE MASS_LAYER_E_BUDGET(B, PTSTEP, PFLX_BLD_MA, PDQS_MA, PIMB_MA, PRADHT_IN, &
                                   PRAD_WL_MA, PRAD_RF_MA, PRAD_WIN_MA, PLOAD_MA,       &
                                   PRAD_FL_MA, PCONV_MA_BLD                  )
USE MODD_BEM_n, ONLY : BEM_t
IMPLICIT NONE
TYPE(BEM_t), INTENT(INOUT) :: B
REAL,                 INTENT(IN)    :: PTSTEP       ! time step
REAL, DIMENSION(:),   INTENT(OUT)  :: PFLX_BLD_MA !flux from building to floor
REAL, DIMENSION(:),   INTENT(OUT) :: PDQS_MA !heat storage inside the floor
REAL, DIMENSION(:),   INTENT(OUT) :: PIMB_MA !floor energy residual imbalance for verification
REAL, DIMENSION(:), INTENT(IN)    :: PRADHT_IN      ! Indoor radiant heat transfer coefficient
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_RF_MA ! rad. fluxes from roof to floor[W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_WL_MA ! rad. fluxes from wall to floor[W m-2(wall)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_WIN_MA  ! rad. fluxes from win to floor[W m-2(win)]
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_MA ! solar and internal load to the floor
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_FL_MA  ! rad. fluxes from floor to mass [W m-2(floor)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_MA_BLD  ! conv. fluxes from floor to bld [W m-2(floor)]
END SUBROUTINE MASS_LAYER_E_BUDGET
END INTERFACE
END MODULE MODI_MASS_LAYER_E_BUDGET
