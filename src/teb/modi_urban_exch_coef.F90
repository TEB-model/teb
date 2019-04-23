!auto_modi:spll_urban_exch_coef.D
MODULE MODI_URBAN_EXCH_COEF
INTERFACE
SUBROUTINE URBAN_EXCH_COEF(HZ0H, PZ0_O_Z0H, PTG, PQS, PEXNS, PEXNA, PTA, PQA,   &
                             PZREF, PUREF, PVMOD, PZ0,                            &
                             PRI, PCD, PCDN, PAC, PRA, PCH                        )  
IMPLICIT NONE
 CHARACTER(LEN=6)                  :: HZ0H     ! TEB option for z0h roof & road
REAL,               INTENT(IN)    :: PZ0_O_Z0H! z0/z0h ratio used in Mascart (1995)
REAL, DIMENSION(:), INTENT(IN)    :: PTG      ! surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PQS      ! surface specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS    ! surface exner function
REAL, DIMENSION(:), INTENT(IN)    :: PTA      ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA      ! specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PEXNA    ! exner function
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD    ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PZ0      ! roughness length for momentum
REAL, DIMENSION(:), INTENT(IN)    :: PZREF    ! reference height of the first
REAL, DIMENSION(:), INTENT(IN)    :: PUREF    ! reference height of the wind
REAL, DIMENSION(:), INTENT(OUT)   :: PRI      ! Richardson number
REAL, DIMENSION(:), INTENT(OUT)   :: PCD      ! drag coefficient for momentum
REAL, DIMENSION(:), INTENT(OUT)   :: PCDN     ! neutral drag coefficient for momentum
REAL, DIMENSION(:), INTENT(OUT)   :: PAC      ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PRA      ! aerodynamical resistance
REAL, DIMENSION(:), INTENT(OUT)   :: PCH      ! drag coefficient for heat
END SUBROUTINE URBAN_EXCH_COEF
END INTERFACE
END MODULE MODI_URBAN_EXCH_COEF
