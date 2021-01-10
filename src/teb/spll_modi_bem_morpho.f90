!auto_modi:spll_bem_morpho.D
MODULE MODI_BEM_MORPHO
INTERFACE
SUBROUTINE BEM_MORPHO(PBLD, PWALL_O_HOR, PBLD_HEIGHT, PWALL_O_BLD, B)
USE MODD_BEM_n, ONLY : BEM_t
IMPLICIT NONE
REAL, DIMENSION(:),   INTENT(IN)  :: PBLD         ! Urban horizontal building density
REAL, DIMENSION(:),   INTENT(IN)  :: PWALL_O_HOR  ! Wall to horizontal surface ratio
REAL, DIMENSION(:),   INTENT(IN)  :: PBLD_HEIGHT  ! Average building height [m]
REAL, DIMENSION(:),   INTENT(OUT) :: PWALL_O_BLD  ! wall surface per ground building surface [m2(wall)/m2(bld)]
TYPE(BEM_t), INTENT(INOUT) :: B
END SUBROUTINE BEM_MORPHO
END INTERFACE
END MODULE MODI_BEM_MORPHO
