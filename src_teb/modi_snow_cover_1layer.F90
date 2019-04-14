!auto_modi:spll_snow_cover_1layer.D
MODULE MODI_SNOW_COVER_1LAYER
INTERFACE
    SUBROUTINE SNOW_COVER_1LAYER(PTSTEP, PANSMIN, PANSMAX, PTODRY,         &
                                   PRHOSMIN, PRHOSMAX, PRHOFOLD, OALL_MELT,  &
                                   PDRAIN_TIME, PWCRN, PZ0SN, PZ0HSN,        &
                                   PTSNOW, PASNOW, PRSNOW, PWSNOW, PTS_SNOW, &
                                   PESNOW,                                   &
                                   PTG, PTG_COEFA, PTG_COEFB,                &
                                   PABS_SW, PLW1, PLW2,                      &
                                   PTA, PQA, PVMOD, PPS, PRHOA, PSR,         &
                                   PZREF, PUREF,                             &
                                   PRNSNOW, PHSNOW, PLESNOW, PGSNOW, PMELT,  &
                                   PDQS_SNOW, PABS_LW                        )  
IMPLICIT NONE
REAL,                 INTENT(IN)    :: PTSTEP   ! time step
REAL,                 INTENT(IN)    :: PANSMIN  ! minimum snow albedo
REAL,                 INTENT(IN)    :: PANSMAX  ! maximum snow albedo
REAL,                 INTENT(IN)    :: PTODRY   ! snow albedo decreasing constant
REAL,                 INTENT(IN)    :: PRHOSMIN ! minimum snow density
REAL,                 INTENT(IN)    :: PRHOSMAX ! maximum snow density
REAL,                 INTENT(IN)    :: PRHOFOLD ! snow density increasing constant
LOGICAL,              INTENT(IN)    :: OALL_MELT! T --> all snow runs off if
REAL,                 INTENT(IN)    :: PDRAIN_TIME ! drainage folding time (days)
REAL,                 INTENT(IN)    :: PWCRN    ! critical snow amount necessary
REAL,                 INTENT(IN)    :: PZ0SN    ! snow roughness length for momentum
REAL,                 INTENT(IN)    :: PZ0HSN   ! snow roughness length for heat
REAL, DIMENSION(:), INTENT(INOUT) :: PWSNOW   ! snow reservoir (kg/m2)
REAL, DIMENSION(:), INTENT(INOUT) :: PTSNOW   ! snow temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PASNOW   ! snow albedo
REAL, DIMENSION(:), INTENT(INOUT) :: PRSNOW   ! snow density
REAL, DIMENSION(:), INTENT(INOUT) :: PTS_SNOW ! snow surface temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PESNOW   ! snow emissivity
REAL, DIMENSION(:), INTENT(IN)    :: PTG      ! underlying ground temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTG_COEFA! underlying ground temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTG_COEFB! implicit terms
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW  ! absorbed SW energy (Wm-2)
REAL, DIMENSION(:), INTENT(IN)    :: PLW1     ! LW coef independant of TSNOW
REAL, DIMENSION(:), INTENT(IN)    :: PLW2     ! LW coef dependant   of TSNOW
REAL, DIMENSION(:), INTENT(IN)    :: PTA      ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA      ! specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD    ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PPS      ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA    ! air density
REAL, DIMENSION(:), INTENT(IN)    :: PSR      ! snow rate
REAL, DIMENSION(:), INTENT(IN)    :: PZREF    ! reference height of the first
REAL, DIMENSION(:), INTENT(IN)    :: PUREF    ! reference height of the first
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSNOW  ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSNOW   ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESNOW  ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSNOW   ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT    ! snow melting rate (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_SNOW! heat storage inside snow
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW ! absorbed LW rad by snow (W/m2)
END SUBROUTINE SNOW_COVER_1LAYER
END INTERFACE
END MODULE MODI_SNOW_COVER_1LAYER
