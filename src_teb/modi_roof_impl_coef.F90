!auto_modi:spll_roof_impl_coef.D
MODULE MODI_ROOF_IMPL_COEF
INTERFACE
SUBROUTINE ROOF_IMPL_COEF(PTSTEP, KROOF_LAYER, PD_ROOF, PTC_ROOF, PHC_ROOF, PT_ROOF, PTDEEP_A,PTDEEP_B)
IMPLICIT NONE
REAL              ,   INTENT(IN)  :: PTSTEP      ! time step
INTEGER           ,   INTENT(IN)  :: KROOF_LAYER ! number of roof layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_ROOF   ! thickness            of each layer
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_ROOF  ! thermal conductivity of each layer
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_ROOF  ! heat capacity        of each layer
REAL, DIMENSION(:,:), INTENT(IN)  :: PT_ROOF   ! temperature          of each layer
REAL, DIMENSION(:),   INTENT(OUT) :: PTDEEP_A, PTDEEP_B     
END SUBROUTINE ROOF_IMPL_COEF
END INTERFACE
END MODULE MODI_ROOF_IMPL_COEF
