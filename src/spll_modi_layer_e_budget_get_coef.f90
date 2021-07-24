!auto_modi:spll_layer_e_budget_get_coef.D
MODULE MODI_LAYER_E_BUDGET_GET_COEF
INTERFACE
    SUBROUTINE LAYER_E_BUDGET_GET_COEF(PT, PTSTEP, PIMPL, PHC, PTC, PD, PA, PB, PC, PY)
IMPLICIT NONE
REAL, DIMENSION(:,:), INTENT(IN) :: PT     ! floor layers temperatures
REAL,                 INTENT(IN) :: PTSTEP ! time step
REAL, DIMENSION(:,:), INTENT(IN) :: PHC    ! heat capacity for road layers
REAL, DIMENSION(:,:), INTENT(IN) :: PTC    ! thermal conductivity for road layers
REAL, DIMENSION(:,:), INTENT(IN) :: PD     ! depth of road layers
REAL, DIMENSION(:,:), INTENT(OUT) :: PA 
REAL, DIMENSION(:,:), INTENT(OUT) :: PB 
REAL, DIMENSION(:,:), INTENT(OUT) :: PC
REAL, DIMENSION(:,:), INTENT(OUT)  :: PY
REAL,                 INTENT(IN)   :: PIMPL  ! implicitation coefficient
END SUBROUTINE LAYER_E_BUDGET_GET_COEF
END INTERFACE
END MODULE MODI_LAYER_E_BUDGET_GET_COEF
