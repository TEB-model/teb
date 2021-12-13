!auto_modi:spll_surface_cd.D
MODULE MODI_SURFACE_CD
INTERFACE
    SUBROUTINE SURFACE_CD(PRI, PZREF, PUREF, PZ0EFF, PZ0H,   &
                              PCD, PCDN) 
IMPLICIT NONE
REAL, DIMENSION(:), INTENT(IN)    :: PRI      ! Richardson number
REAL, DIMENSION(:), INTENT(IN)    :: PZREF    ! reference height of the first
REAL, DIMENSION(:), INTENT(IN)    :: PUREF    ! reference height of the wind
REAL, DIMENSION(:), INTENT(IN)    :: PZ0EFF   ! roughness length for momentum
REAL, DIMENSION(:), INTENT(IN)    :: PZ0H     ! roughness length for heat
REAL, DIMENSION(:), INTENT(OUT)   :: PCD      ! drag coefficient for momentum
REAL, DIMENSION(:), INTENT(OUT)   :: PCDN     ! neutral drag coefficient for momentum
END SUBROUTINE SURFACE_CD
END INTERFACE
END MODULE MODI_SURFACE_CD
