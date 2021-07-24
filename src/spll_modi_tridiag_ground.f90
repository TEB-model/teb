!auto_modi:spll_tridiag_ground.D
MODULE MODI_TRIDIAG_GROUND
INTERFACE
       SUBROUTINE TRIDIAG_GROUND(PA,PB,PC,PY,PX)
IMPLICIT NONE
REAL,    DIMENSION(:,:), INTENT(IN)  :: PA  ! lower diag. elements of A matrix
REAL,    DIMENSION(:,:), INTENT(IN)  :: PB  ! main  diag. elements of A matrix
REAL,    DIMENSION(:,:), INTENT(IN)  :: PC  ! upper diag. elements of A matrix
REAL,    DIMENSION(:,:), INTENT(IN)  :: PY  ! r.h.s. term   
REAL,    DIMENSION(:,:), INTENT(OUT) :: PX  ! solution of A.X = Y 
END SUBROUTINE TRIDIAG_GROUND
END INTERFACE
END MODULE MODI_TRIDIAG_GROUND
