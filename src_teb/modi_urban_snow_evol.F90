!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_URBAN_SNOW_EVOL
INTERFACE
    SUBROUTINE URBAN_SNOW_EVOL(                                                 &
                       PT_LOWCAN, PQ_LOWCAN, PU_LOWCAN,                         &
                       PTS_ROOF,PTS_ROAD,PTS_WALL_A, PTS_WALL_B,                &
                       PT_ROOF, PD_ROOF, PTC_ROOF, PHC_ROOF,                    &
                       HSNOW_ROOF,                                              &
                       PWSNOW_ROOF, PTSNOW_ROOF, PRSNOW_ROOF, PASNOW_ROOF,      &
                       PTSSNOW_ROOF, PESNOW_ROOF,                               &
                       HSNOW_ROAD,                                              &
                       PWSNOW_ROAD, PTSNOW_ROAD, PRSNOW_ROAD, PASNOW_ROAD,      &
                       PTSSNOW_ROAD, PESNOW_ROAD,                               &
                       PPS, PTA, PQA, PRHOA,                                    &
                       PLW_RAD,                                                 &
                       PSR, PZREF, PUREF, PVMOD,                                &
                       PTSTEP,                                                  &
                       PZ_LOWCAN,                                               &
                       PDN_ROOF, PABS_SW_SNOW_ROOF, PABS_LW_SNOW_ROOF,          &
                       PDN_ROAD, PABS_SW_SNOW_ROAD, PABS_LW_SNOW_ROAD,          &
                       PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,    &
                       PMELT_ROOF,                                              &
                       PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,    &
                       PMELT_ROAD,                                              &
                       PLW_WA_TO_NR , PLW_WB_TO_NR, PLW_S_TO_NR, PLW_WIN_TO_NR, &
                       PDQS_SNOW_ROOF, PDQS_SNOW_ROAD, PT_WIN1                  )  
REAL, DIMENSION(:),   INTENT(IN)    :: PT_LOWCAN  ! LOWCAN air temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PQ_LOWCAN  ! LOWCAN air specific humidity
REAL, DIMENSION(:),   INTENT(IN)    :: PU_LOWCAN  ! LOWCAN hor. wind
REAL, DIMENSION(:),   INTENT(IN)    :: PTS_ROOF   ! roof surface temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PTS_ROAD   ! road surface temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PTS_WALL_A ! wall surface temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PTS_WALL_B ! wall surface temperature
REAL, DIMENSION(:,:), INTENT(IN)    :: PT_ROOF    ! roof temperature profile
REAL, DIMENSION(:,:), INTENT(IN)    :: PD_ROOF    ! roof layer thickness
REAL, DIMENSION(:,:), INTENT(IN)    :: PTC_ROOF   ! roof layer thermal conductivity
REAL, DIMENSION(:,:), INTENT(IN)    :: PHC_ROOF   ! roof layer heat capacity
 CHARACTER(LEN=*),     INTENT(IN)    :: HSNOW_ROOF ! snow roof scheme
 CHARACTER(LEN=*),     INTENT(IN)    :: HSNOW_ROAD ! snow road scheme
REAL, DIMENSION(:,:), INTENT(INOUT) :: PWSNOW_ROOF ! snow layers reservoir
REAL, DIMENSION(:,:), INTENT(INOUT) :: PTSNOW_ROOF ! snow layers temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRSNOW_ROOF ! snow layers density
REAL, DIMENSION(:),   INTENT(INOUT) :: PASNOW_ROOF ! snow albedo
REAL, DIMENSION(:),   INTENT(INOUT) :: PESNOW_ROOF ! snow emissivity
REAL, DIMENSION(:),   INTENT(INOUT) :: PTSSNOW_ROOF! snow surface temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PWSNOW_ROAD ! snow layers reservoir
REAL, DIMENSION(:,:), INTENT(INOUT) :: PTSNOW_ROAD ! snow layers temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRSNOW_ROAD ! snow layers density
REAL, DIMENSION(:),   INTENT(INOUT) :: PASNOW_ROAD ! snow albedo
REAL, DIMENSION(:),   INTENT(INOUT) :: PESNOW_ROAD ! snow emissivity
REAL, DIMENSION(:),   INTENT(INOUT) :: PTSSNOW_ROAD! snow surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PPS      ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PTA      ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA      ! specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD    ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA    ! air density at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD  ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PSR      ! snow rate
REAL, DIMENSION(:), INTENT(IN)    :: PZREF    ! reference height of the first
REAL, DIMENSION(:), INTENT(IN)    :: PUREF    ! reference height of the first
REAL,               INTENT(IN)    :: PTSTEP   ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PZ_LOWCAN  ! height of forcing
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROOF          ! snow-covered roof frac.
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_SNOW_ROOF ! SW absorbed by roof snow
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_SNOW_ROOF ! absorbed IR rad by snow on roof
REAL, DIMENSION(:), INTENT(INOUT) :: PDN_ROAD          ! snow-covered road frac.
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_SNOW_ROAD ! SW absorbed by road snow
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_SNOW_ROAD ! absorbed IR rad by snow on road
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSNOW_ROOF ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSNOW_ROOF  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESNOW_ROOF ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSNOW_ROOF  ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_ROOF   ! snow melt
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSNOW_ROAD ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSNOW_ROAD  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESNOW_ROAD ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSNOW_ROAD  ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_ROAD   ! snow melt
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_NR        ! LW contrib. wall       -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WB_TO_NR        ! LW contrib. wall       -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_NR         ! LW contrib. sky        -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WIN_TO_NR       ! LW contrib. win       -> road(snow)
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_SNOW_ROOF ! Heat storage in snowpack on roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_SNOW_ROAD ! Heat storage in snowpack on roads
REAL, DIMENSION(:), INTENT(IN)    :: PT_WIN1        ! Window surface temperature
END SUBROUTINE URBAN_SNOW_EVOL
END INTERFACE
END MODULE MODI_URBAN_SNOW_EVOL
