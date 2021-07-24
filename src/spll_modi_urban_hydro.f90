!auto_modi:spll_urban_hydro.D
MODULE MODI_URBAN_HYDRO
INTERFACE
    SUBROUTINE URBAN_HYDRO(PWS_ROOF_MAX,PWS_ROAD_MAX, PWS_ROOF, PWS_ROAD,  &
                             PRR, PIRRIG_ROAD, PTSTEP, PBLD, PLE_ROOF,     &
                             PLE_ROAD, PRUNOFF_STRLROOF, PRUNOFF_ROAD      )
IMPLICIT NONE
REAL, DIMENSION(:), INTENT(IN) :: PWS_ROOF_MAX    ! maximum deepness of roof water reservoir
REAL, DIMENSION(:), INTENT(IN) :: PWS_ROAD_MAX    ! maximum deepness of road water reservoir
REAL, DIMENSION(:), INTENT(INOUT) :: PWS_ROOF     ! roof water reservoir
REAL, DIMENSION(:), INTENT(INOUT) :: PWS_ROAD     ! road water reservoir
REAL, DIMENSION(:), INTENT(IN)    :: PRR          ! rain rate
REAL, DIMENSION(:), INTENT(IN)    :: PIRRIG_ROAD  ! watering rate for roads
REAL,               INTENT(IN)    :: PTSTEP       ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PBLD         ! fraction of buildings
REAL, DIMENSION(:), INTENT(IN)    :: PLE_ROOF     ! latent heat flux over roof
REAL, DIMENSION(:), INTENT(IN)    :: PLE_ROAD     ! latent heat flux over road
REAL, DIMENSION(:), INTENT(OUT)   :: PRUNOFF_STRLROOF ! runoff (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT)   :: PRUNOFF_ROAD ! runoff (kg/m2/s)
END SUBROUTINE URBAN_HYDRO
END INTERFACE
END MODULE MODI_URBAN_HYDRO
