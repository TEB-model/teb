!auto_modi:spll_surface_aero_cond.D
MODULE MODI_SURFACE_AERO_COND
INTERFACE
    SUBROUTINE SURFACE_AERO_COND(PRI, PZREF, PUREF, PVMOD, PZ0,&
                                     PZ0H, PAC, PRA, PCH           ) 
IMPLICIT NONE
REAL, DIMENSION(:), INTENT(IN)    :: PRI      ! Richardson number
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD    ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PZREF    ! reference height of the first
REAL, DIMENSION(:), INTENT(IN)    :: PUREF    ! reference height of the wind
REAL, DIMENSION(:), INTENT(IN)    :: PZ0      ! roughness length for momentum
REAL, DIMENSION(:), INTENT(IN)    :: PZ0H     ! roughness length for heat
REAL, DIMENSION(:), INTENT(OUT)   :: PAC      ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PRA      ! aerodynamical resistance
REAL, DIMENSION(:), INTENT(OUT)   :: PCH      ! drag coefficient for heat
END SUBROUTINE SURFACE_AERO_COND
END INTERFACE
END MODULE MODI_SURFACE_AERO_COND
