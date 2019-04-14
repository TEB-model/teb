!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_URBAN_DRAG
INTERFACE
    SUBROUTINE URBAN_DRAG(HZ0H, HIMPLICIT_WIND, PTSTEP,                     &
                          PT_CANYON, PQ_CANYON, PU_CANYON,                  &
                          PT_LOWCAN, PQ_LOWCAN, PU_LOWCAN, PZ_LOWCAN,       &
                          PTS_ROOF, PTS_ROAD, PTS_WALL, PTS_GARDEN,         &
                          PDELT_SNOW_ROOF, PDELT_SNOW_ROAD,                 &
                          PEXNS, PEXNA, PTA, PQA, PPS, PRHOA,               &
                          PZREF, PUREF, PVMOD,                              &
                          PZ0_TOWN,                                         &
                          PBLD, PGARDEN, PROAD,                             &
                          PBLD_HEIGHT, PCAN_HW_RATIO,                       &
                          PWALL_O_GRND,                                     &
                          PWS_ROOF, PWS_ROAD,                               &
                          PWS_ROOF_MAX, PWS_ROAD_MAX,                       &
                          PPEW_A_COEF, PPEW_B_COEF,                         &
                          PPEW_A_COEF_LOWCAN, PPEW_B_COEF_LOWCAN,           &
                          PQSAT_ROOF, PQSAT_ROAD, PDELT_ROOF, PDELT_ROAD,   &
                          PCD, PCDN, PAC_ROOF, PAC_ROOF_WAT,                &
                          PAC_WALL, PAC_ROAD, PAC_ROAD_WAT, PAC_TOP,        &
                          PAC_GARDEN, PRI,                                  &
                          PUW_ROAD, PUW_ROOF, PDUWDU_ROAD, PDUWDU_ROOF,     &
                          PUSTAR_TOWN, OCANOPY, PTS_WIN, PAC_WIN, HCH_BEM,  &
                          PROUGH_ROOF, PROUGH_WALL                          ) 
 CHARACTER(LEN=6),   INTENT(IN)    :: HZ0H           ! TEB option for z0h roof & road
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
REAL,               INTENT(IN)    :: PTSTEP         ! time-step
REAL, DIMENSION(:), INTENT(IN)    :: PT_CANYON      ! canyon air temperature
REAL, DIMENSION(:), INTENT(IN)    :: PQ_CANYON      ! canyon air specific humidity.
REAL, DIMENSION(:), INTENT(IN)    :: PU_CANYON      ! hor. wind in canyon
REAL, DIMENSION(:), INTENT(IN)    :: PU_LOWCAN     ! wind near the road
REAL, DIMENSION(:), INTENT(IN)    :: PT_LOWCAN     ! temp. near the road
REAL, DIMENSION(:), INTENT(IN)    :: PQ_LOWCAN     ! hum. near the road
REAL, DIMENSION(:), INTENT(IN)    :: PZ_LOWCAN     ! height of atm. var. near the road
REAL, DIMENSION(:), INTENT(IN)    :: PTS_ROOF       ! surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_ROAD       ! surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_WALL       ! surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GARDEN     ! surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_SNOW_ROOF! fraction of snow on roof
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_SNOW_ROAD! fraction of snow on road
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS          ! surface exner function
REAL, DIMENSION(:), INTENT(IN)    :: PTA            ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA            ! specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD          ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PPS            ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PEXNA          ! exner function
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA          ! air density
REAL, DIMENSION(:), INTENT(IN)    :: PZREF          ! reference height of the first
REAL, DIMENSION(:), INTENT(IN)    :: PUREF          ! reference height of the first
REAL, DIMENSION(:), INTENT(IN)    :: PZ0_TOWN       ! roughness length for momentum
REAL, DIMENSION(:), INTENT(IN)    :: PBLD           ! fraction of buildings
REAL, DIMENSION(:), INTENT(IN)    :: PROAD          ! fraction of roads
REAL, DIMENSION(:), INTENT(IN)    :: PGARDEN        ! fraction of GARDEN areas
REAL, DIMENSION(:), INTENT(IN)    :: PBLD_HEIGHT    ! h
REAL, DIMENSION(:), INTENT(IN)    :: PCAN_HW_RATIO  ! h/W
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_O_GRND   ! wall surf. / (road+GARDEN area) surf.
REAL, DIMENSION(:), INTENT(IN)    :: PWS_ROOF       ! roof water content (kg/m2)
REAL, DIMENSION(:), INTENT(IN)    :: PWS_ROAD       ! road water content (kg/m2)
REAL, DIMENSION(:), INTENT(IN)    :: PWS_ROOF_MAX   ! maximum deepness of roof
REAL, DIMENSION(:), INTENT(IN)    :: PWS_ROAD_MAX   ! and water reservoirs (kg/m2)
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_A_COEF    ! implicit coefficients (m2s/kg)
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_B_COEF    ! for wind coupling     (m/s)
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_A_COEF_LOWCAN ! implicit coefficients for wind coupling (m2s/kg)
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_B_COEF_LOWCAN ! between low canyon wind and road (m/s)
REAL, DIMENSION(:), INTENT(OUT)   :: PQSAT_ROOF     ! qsat(Ts)
REAL, DIMENSION(:), INTENT(OUT)   :: PQSAT_ROAD     ! qsat(Ts)
REAL, DIMENSION(:), INTENT(OUT)   :: PDELT_ROOF     ! water fraction on
REAL, DIMENSION(:), INTENT(OUT)   :: PDELT_ROAD     ! snow-free surfaces
REAL, DIMENSION(:), INTENT(OUT)   :: PCD            ! drag coefficient
REAL, DIMENSION(:), INTENT(OUT)   :: PCDN           ! neutral drag coefficient
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROOF       ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROOF_WAT   ! aerodynamical conductance (for water)
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_WALL       ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROAD       ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROAD_WAT   ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_TOP        ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(IN)    :: PAC_GARDEN     ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PRI            ! Town Richardson number
REAL, DIMENSION(:), INTENT(OUT)   :: PUW_ROAD       ! Momentum flux for roads
REAL, DIMENSION(:), INTENT(OUT)   :: PUW_ROOF       ! Momentum flux for roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PDUWDU_ROAD    ! 
REAL, DIMENSION(:), INTENT(OUT)   :: PDUWDU_ROOF    ! 
REAL, DIMENSION(:), INTENT(OUT)   :: PUSTAR_TOWN    ! Fraction velocity for town
LOGICAL, INTENT(IN)               :: OCANOPY        ! is canopy active
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_WIN        ! aerodynamical conductance for window
REAL, DIMENSION(:), INTENT(IN)    :: PTS_WIN        ! window outdoor surface temp
 CHARACTER(LEN=5), INTENT(IN)      :: HCH_BEM        ! BEM option for convective heat transfer coef.
REAL, DIMENSION(:), INTENT(IN)    :: PROUGH_ROOF    ! roughness coef for the roof
REAL, DIMENSION(:), INTENT(IN)    :: PROUGH_WALL    ! roughness coef for the wall
END SUBROUTINE URBAN_DRAG
END INTERFACE
END MODULE MODI_URBAN_DRAG
