!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_URBAN_FLUXES
INTERFACE
    SUBROUTINE URBAN_FLUXES(HIMPLICIT_WIND, OCANOPY, PT_CANYON,                 &
                       PPEW_A_COEF, PPEW_B_COEF,                                &
                       PEXNS, PRHOA, PVMOD,                                     &
                       PH_TRAFFIC, PLE_TRAFFIC, PH_INDUSTRY, PLE_INDUSTRY,      &
                       PBLD, PROAD, PWALL_O_HOR, PFRAC_GR,                      &
                       PABS_SW_ROOF, PABS_LW_ROOF,                              &
                       PABS_SW_WALL_A, PABS_LW_WALL_A,                          &
                       PABS_SW_WALL_B, PABS_LW_WALL_B,                          &
                       PABS_SW_ROAD, PABS_LW_ROAD,                              &
                       PAC_WALL, PCD,                                           &
                       PDF_ROOF, PDN_ROOF, PDF_ROAD, PDN_ROAD,                  &
                       PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,    &
                       PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,     &
                       PMELT_ROOF, PDQS_ROOF,                                   &
                       PMELT_ROAD, PDQS_ROAD,                                   &
                       PDQS_WALL_A, PDQS_WALL_B,                                &
                       PFLX_BLD_ROOF, PFLX_BLD_WALL_A, PFLX_BLD_WALL_B,         &
                       PFLX_BLD_FLOOR,                                          &
                       PFLX_BLD_MASS, PE_SHADING, HBEM, PGR, PT_WIN1, PH_WASTE, &
                       PLE_WASTE, PF_WASTE_CAN, PABS_SW_WIN, PABS_LW_WIN, PQIN, &
                       PN_FLOOR, PLEW_ROOF, PLEW_ROAD,                          &
                       PRN_GREENROOF, PH_GREENROOF, PLE_GREENROOF,              &
                       PGFLUX_GREENROOF,                                        &
                       PH_ROOF, PH_ROAD, PH_WALL_A, PH_WALL_B,                  &
                       PRN_ROOF, PLE_ROOF,                                      &
                       PGFLUX_ROOF, PRN_ROAD, PLE_ROAD, PGFLUX_ROAD,            &
                       PRN_WALL_A, PLE_WALL_A, PGFLUX_WALL_A,                   &
                       PRN_WALL_B, PLE_WALL_B, PGFLUX_WALL_B,                   &
                       PRN_STRLROOF, PH_STRLROOF, PLE_STRLROOF,                  &
                       PGFLUX_STRLROOF,                                         &
                       PRN_BLT,  PH_BLT,  PLE_BLT,                              &
                       PGFLUX_BLT, PMELT_BLT, PQF_BLD, PFLX_BLD, PDQS_BLT,      &
                       PQF_TOWN, PUSTAR_TOWN, PHVAC_COOL, PHVAC_HEAT            )
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
LOGICAL,            INTENT(IN)    :: OCANOPY      ! is canopy active ?
REAL, DIMENSION(:), INTENT(IN)    :: PT_CANYON    ! canyon air temperature
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_A_COEF  ! implicit coefficients
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_B_COEF  ! for wind coupling
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS        ! surface exner function
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA        ! air density
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD        ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PH_TRAFFIC   ! anthropogenic sensible
REAL, DIMENSION(:), INTENT(IN)    :: PLE_TRAFFIC  ! anthropogenic latent
REAL, DIMENSION(:), INTENT(IN)    :: PH_INDUSTRY  ! anthropogenic sensible
REAL, DIMENSION(:), INTENT(IN)    :: PLE_INDUSTRY ! anthropogenic latent
REAL, DIMENSION(:), INTENT(IN)    :: PBLD         ! fraction of buildings
REAL, DIMENSION(:), INTENT(IN)    :: PROAD        ! fraction of roads
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_O_HOR  ! wall Surf. / (bld+road+green) Surf.
REAL, DIMENSION(:), INTENT(IN)    :: PFRAC_GR     ! Fraction of green roofs
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROOF ! absorbed SW rad. by roof
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_ROOF ! absorbed LW rad. by roof
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WALL_A! absorbed SW rad. by wall A
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_WALL_A! absorbed LW rad. by wall A
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WALL_B! absorbed SW rad. by wall B
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_WALL_B! absorbed LW rad. by wall B
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROAD ! absorbed SW rad. by road
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_ROAD ! absorbed LW rad. by road
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WALL     ! surface conductance
REAL, DIMENSION(:), INTENT(IN)    :: PCD          ! drag coefficient
REAL, DIMENSION(:), INTENT(IN)    :: PDF_ROOF     ! snow-free    roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROOF     ! snow-covered roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDF_ROAD     ! snow-free    road fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROAD     ! snow-covered road fraction
REAL, DIMENSION(:), INTENT(IN)    :: PRNSNOW_ROOF ! net radiation over snow
REAL, DIMENSION(:), INTENT(IN)    :: PHSNOW_ROOF  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PLESNOW_ROOF ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PGSNOW_ROOF  ! flux under the snow
REAL, DIMENSION(:), INTENT(IN)    :: PRNSNOW_ROAD ! net radiation over snow
REAL, DIMENSION(:), INTENT(IN)    :: PHSNOW_ROAD  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PLESNOW_ROAD ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PGSNOW_ROAD  ! flux under the snow
REAL, DIMENSION(:), INTENT(IN)    :: PMELT_ROOF   ! snow melting on roof
REAL, DIMENSION(:), INTENT(IN)    :: PDQS_ROOF    ! storage inside roofs
REAL, DIMENSION(:), INTENT(IN)    :: PMELT_ROAD   ! snow melting on road
REAL, DIMENSION(:), INTENT(IN)    :: PDQS_ROAD    ! storage inside roads
REAL, DIMENSION(:), INTENT(IN)    :: PDQS_WALL_A  ! storage inside walls
REAL, DIMENSION(:), INTENT(IN)    :: PDQS_WALL_B  ! storage inside walls
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_ROOF! heat flx from bld to roof
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_WALL_A! heat flx from bld to wall A
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_WALL_B! heat flx from bld to wall B
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_FLOOR! heat flx from bld to floor
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_MASS! heat flx from bld to mass
REAL, DIMENSION(:), INTENT(IN)    :: PE_SHADING   ! energy not ref., nor absorbed, nor
 CHARACTER(LEN=3), INTENT(IN)      :: HBEM         ! Building Energy model 'DEF' or 'BEM'
REAL, DIMENSION(:), INTENT(IN)    :: PGR          ! glazing ratio
REAL, DIMENSION(:), INTENT(IN)    :: PT_WIN1      ! Outdoor window temperature [K]
REAL, DIMENSION(:), INTENT(IN)    :: PH_WASTE     ! sensible waste heat released by HVAC systems [W m-2(urb)]
REAL, DIMENSION(:), INTENT(IN)    :: PLE_WASTE    ! latent waste heat released by HVAC systems [W m-2(urb)]
REAL, DIMENSION(:), INTENT(IN)    :: PF_WASTE_CAN ! fraction of waste heat released into the canyon
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WIN  ! absorbed SW radiation by windows [W m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_WIN  ! absorbed LW radiation by windows [W m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PQIN         ! Internal heat gains [W m-2(floor)]
REAL, DIMENSION(:), INTENT(IN)    :: PN_FLOOR     ! Number of Building floors [-]
REAL, DIMENSION(:), INTENT(IN)   :: PLEW_ROOF     ! latent heat flux over snow-free roof
REAL, DIMENSION(:), INTENT(IN)   :: PLEW_ROAD     ! latent heat flux of snow-free road
REAL, DIMENSION(:), INTENT(IN)   :: PRN_GREENROOF ! net radiation of green roofs
REAL, DIMENSION(:), INTENT(IN)   :: PH_GREENROOF  ! sensible heat flux of green roofs
REAL, DIMENSION(:), INTENT(IN)   :: PLE_GREENROOF ! latent heat flux of green roofs
REAL, DIMENSION(:), INTENT(IN)   :: PGFLUX_GREENROOF ! heat storage flux of green roofs
REAL, DIMENSION(:), INTENT(INOUT) :: PH_ROOF      ! sensible heat flux over roof
REAL, DIMENSION(:), INTENT(INOUT) :: PH_ROAD      ! sensible heat flux from snow free road
REAL, DIMENSION(:), INTENT(INOUT) :: PH_WALL_A    ! sensible heat flux wall to
REAL, DIMENSION(:), INTENT(INOUT) :: PH_WALL_B    ! sensible heat flux wall to
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_ROOF     ! net radiation over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_ROOF     ! latent heat flux over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_ROOF  ! flux through the roof
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_ROAD     ! net radiation over road
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_ROAD     ! latent heat flux over road
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_ROAD  ! flux through the road
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_WALL_A   ! net radiation over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WALL_A   ! latent heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_WALL_A! flux through the wall
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_WALL_B   ! net radiation over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WALL_B   ! latent heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_WALL_B! flux through the wall
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_STRLROOF  ! net radiation of structural roof
REAL, DIMENSION(:), INTENT(OUT)   :: PH_STRLROOF   ! sensible heat flux of structural roof
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_STRLROOF  ! latent heat flux of structural roof
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_STRLROOF ! heat storage flux of structural roof
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_BLT      ! net radiation over built covers
REAL, DIMENSION(:), INTENT(OUT)   :: PH_BLT       ! sensible heat flux over built covers
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_BLT      ! latent heat flux over built covers
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_BLT   ! flux through the built covers
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_BLT    ! snow melting for town
REAL, DIMENSION(:), INTENT(OUT)   :: PQF_BLD      ! anthropogenic heat due to domestic heating [W m-2(bld)]
REAL, DIMENSION(:), INTENT(OUT)   :: PFLX_BLD     ! heat flx from bld to its structure         [W m-2(bld)]
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_BLT     ! storage inside urban materials
REAL, DIMENSION(:), INTENT(OUT)   :: PQF_TOWN     ! total anthropogenic heat                   [W m-2(urb)]
REAL, DIMENSION(:), INTENT(OUT)   :: PUSTAR_TOWN  ! friction velocity over town
REAL, DIMENSION(:), INTENT(IN)   :: PHVAC_COOL   ! Energy consumption of the cooling system
REAL, DIMENSION(:), INTENT(IN)   :: PHVAC_HEAT   ! Energy consumption of the heating system
END SUBROUTINE URBAN_FLUXES
END INTERFACE
END MODULE MODI_URBAN_FLUXES
