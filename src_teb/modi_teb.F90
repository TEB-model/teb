!auto_modi:spll_teb.D
MODULE MODI_TEB
INTERFACE
    SUBROUTINE TEB  (HZ0H, HIMPLICIT_WIND, HWALL_OPT, HBEM, TPTIME, PTSUN,     &
                     PT_CANYON, PQ_CANYON, PU_CANYON,                          &
                     PT_LOWCAN, PQ_LOWCAN, PU_LOWCAN, PZ_LOWCAN, PTI_BLD,      &
                     PT_ROOF, PT_ROAD, PT_WALL_A, PT_WALL_B,PWS_ROOF, PWS_ROAD,&
                     HSNOW_ROOF, PWSNOW_ROOF, PTSNOW_ROOF, PRSNOW_ROOF,        &
                     PASNOW_ROOF, PTSSNOW_ROOF, PESNOW_ROOF,                   &
                     HSNOW_ROAD, PWSNOW_ROAD, PTSNOW_ROAD, PRSNOW_ROAD,        &
                     PASNOW_ROAD, PTSSNOW_ROAD, PESNOW_ROAD,                   &
                     PPEW_A_COEF, PPEW_B_COEF,                                 &
                     PPEW_A_COEF_LOWCAN, PPEW_B_COEF_LOWCAN,                   &
                     PPS, PPA, PEXNS, PEXNA, PTA, PQA, PRHOA, PLW_RAD,         &
                     PRR, PSR, PZREF, PUREF,                                   &
                     PVMOD, PH_TRAFFIC, PLE_TRAFFIC, PH_INDUSTRY, PLE_INDUSTRY,&
                     PTSTEP, PZ0_TOWN, PBLD, PGARDEN, PROAD, PFRAC_GR,         &
                     PBLD_HEIGHT, PWALL_O_HOR, PCAN_HW_RATIO, PWALL_O_GRND,    &
                     PDF_ROOF, PDN_ROOF, PDF_ROAD,                             &
                     PDN_ROAD, PQSAT_ROOF, PQSAT_ROAD, PDELT_ROOF, PDELT_ROAD, &
                     PEMIS_ROOF, PHC_ROOF, PTC_ROOF, PD_ROOF,                  &
                     PEMIS_ROAD, PHC_ROAD, PTC_ROAD,                           &
                     PD_ROAD, PEMIS_WALL, PTS_GARDEN,                          &
                     PHC_WALL, PTC_WALL, PD_WALL, PRN_ROOF, PH_ROOF, PLE_ROOF, &
                     PLEW_ROOF, PGFLUX_ROOF, PRUNOFF_ROOF,                     &
                     PRN_GREENROOF, PH_GREENROOF, PLE_GREENROOF,               &
                     PGFLUX_GREENROOF, PUW_GREENROOF,                          &
                     PRUNOFF_GREENROOF, PDRAIN_GREENROOF,                      &
                     PRN_STRLROOF, PH_STRLROOF, PLE_STRLROOF, PGFLUX_STRLROOF, &
                     PRUNOFF_STRLROOF,                                         &
                     PRN_ROAD, PH_ROAD,                                        &
                     PLE_ROAD, PLEW_ROAD, PGFLUX_ROAD, PRUNOFF_ROAD,           &
                     PRN_WALL_A, PH_WALL_A, PLE_WALL_A, PGFLUX_WALL_A,         &
                     PRN_WALL_B, PH_WALL_B, PLE_WALL_B, PGFLUX_WALL_B,         &
                     PRN_BLT, PH_BLT, PLE_BLT, PGFLUX_BLT,                     &
                     PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,     &
                     PMELT_ROOF,                                               &
                     PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,     &
                     PMELT_ROAD,                                               &
                     PG_GREENROOF_ROOF,                                        &
                     PUW_ROAD, PUW_ROOF, PDUWDU_ROAD, PDUWDU_ROOF,             &
                     PUSTAR_TOWN, PCD, PCDN, PCH_TOWN, PRI_TOWN,               &
                     PRESA_TOWN, PDQS_TOWN, PQF_TOWN, PQF_BLD, PFLX_BLD,       &
                     PAC_ROOF, PAC_ROAD, PAC_WALL, PAC_TOP, PAC_GARDEN,        &
                     PAC_ROOF_WAT, PAC_ROAD_WAT, PABS_SW_ROOF, PABS_LW_ROOF,   &
                     PABS_SW_SNOW_ROOF, PABS_LW_SNOW_ROOF, PABS_SW_ROAD,       &
                     PABS_LW_ROAD, PABS_SW_SNOW_ROAD, PABS_LW_SNOW_ROAD,       &
                     PABS_SW_WALL_A,PABS_LW_WALL_A,                            &
                     PABS_SW_WALL_B,PABS_LW_WALL_B,                            &
                     PLW_WA_TO_WB,                                             &
                     PLW_WA_TO_R, PLW_WB_TO_R,                                 &
                     PLW_WA_TO_NR, PLW_WB_TO_NR,                               &
                     PLW_R_TO_WA, PLW_R_TO_WB,                                 &
                     PLW_G_TO_WA, PLW_G_TO_WB,                                 &
                     PLW_S_TO_WA, PLW_S_TO_WB, PLW_S_TO_R,                     &
                     PLW_S_TO_NR, PLW_NR_TO_WA, PLW_NR_TO_WB,                  &
                     PLW_NR_TO_WIN, PLW_WA_TO_WIN, PLW_WB_TO_WIN, PLW_G_TO_WIN,&
                     PLW_R_TO_WIN, PLW_S_TO_WIN, PLW_WIN_TO_WA, PLW_WIN_TO_WB, &
                     PLW_WIN_TO_R, PLW_WIN_TO_NR,                              &
                     HNATVENT,                                                 &
                     HCOOL_COIL, PF_WATER_COND, HHEAT_COIL, OAUTOSIZE,         &
                     KDAY, PAUX_MAX, PT_FLOOR, PT_MASS, PH_BLD_COOL,           &
                     PT_BLD_COOL, PH_BLD_HEAT, PLE_BLD_COOL, PLE_BLD_HEAT,     &
                     PH_WASTE, PLE_WASTE, PF_WASTE_CAN, PHVAC_COOL, PHVAC_HEAT,&
                     PQIN, PQIN_FRAD, PQIN_FLAT, PGR, PEFF_HEAT,               &
                     PINF, PTCOOL_TARGET, PTHEAT_TARGET, PHR_TARGET, PT_WIN2,  &
                     PQI_BLD, PV_VENT, PCAP_SYS_HEAT, PCAP_SYS_RAT, PT_ADP,    &
                     PM_SYS_RAT, PCOP_RAT, PCAP_SYS, PM_SYS, PCOP,  PQ_SYS,    &
                     PT_SYS, PTR_SW_WIN, PFAN_POWER, PHC_FLOOR, PTC_FLOOR,     &
                     PD_FLOOR,                                                 &
                     PT_WIN1, PABS_SW_WIN, PABS_LW_WIN, PUGG_WIN, PEMIT_LW_FAC,&
                     PEMIT_LW_ROAD, PT_RAD_IND, PHU_BLD, PTIME, PE_SHADING,    &
                     ONATVENT_NIGHT, PN_FLOOR, PWALL_O_BLD, PGLAZ_O_BLD, PMASS_O_BLD, &
                     PFLOOR_HW_RATIO, PF_FLOOR_MASS, PF_FLOOR_WALL,            &
                     PF_FLOOR_WIN, PF_FLOOR_ROOF, PF_WALL_FLOOR, PF_WALL_MASS, &
                     PF_WALL_WIN, PF_WIN_FLOOR, PF_WIN_MASS, PF_WIN_WALL,      &
                     PF_MASS_FLOOR, PF_MASS_WALL, PF_MASS_WIN, OCANOPY,        &
                     HCH_BEM, PROUGH_ROOF, PROUGH_WALL, PF_WIN_WIN,            &
                     OPAR_RD_IRRIG, PRD_START_MONTH, PRD_END_MONTH,            &
                     PRD_START_HOUR, PRD_END_HOUR, PRD_24H_IRRIG, PIRRIG_ROAD  )
USE MODD_TYPE_DATE_SURF,ONLY: DATE_TIME
IMPLICIT NONE
 CHARACTER(LEN=6), INTENT(IN)      :: HZ0H          ! TEB option for z0h roof & road
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
 CHARACTER(LEN=4), INTENT(IN)      :: HWALL_OPT     ! TEB option for walls
 CHARACTER(LEN=3), INTENT(IN)      :: HBEM          ! Building Energy model 'DEF' or 'BEM'
TYPE(DATE_TIME)     , INTENT(IN)    :: TPTIME             ! current date and time from teb
REAL, DIMENSION(:),   INTENT(IN)    :: PTSUN              ! solar time   (s from midnight)
REAL, DIMENSION(:), INTENT(INOUT) :: PT_CANYON     ! canyon air temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PQ_CANYON     ! canyon air specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PU_CANYON     ! canyon hor. wind
REAL, DIMENSION(:), INTENT(IN)    :: PU_LOWCAN     ! wind near the road
REAL, DIMENSION(:), INTENT(IN)    :: PT_LOWCAN     ! temp. near the road
REAL, DIMENSION(:), INTENT(IN)    :: PQ_LOWCAN     ! hum. near the road
REAL, DIMENSION(:), INTENT(IN)    :: PZ_LOWCAN     ! height of atm. var. near the road
REAL, DIMENSION(:), INTENT(INOUT) :: PTI_BLD       ! inside building temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_ROOF     ! roof layers temperatures
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_ROAD     ! road layers temperatures
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_WALL_A   ! wall A layers temperatures
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_WALL_B   ! wall B layers temperatures
REAL, DIMENSION(:), INTENT(INOUT) :: PWS_ROOF      ! roof water reservoir
REAL, DIMENSION(:), INTENT(INOUT) :: PWS_ROAD      ! road water reservoir
 CHARACTER(LEN=*),   INTENT(IN)    :: HSNOW_ROOF    ! snow roof scheme
 CHARACTER(LEN=*),   INTENT(IN)    :: HSNOW_ROAD    ! snow road scheme
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
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_A_COEF   ! implicit coefficients
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_B_COEF   ! for wind coupling
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_A_COEF_LOWCAN ! implicit coefficients for wind coupling
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_B_COEF_LOWCAN ! between low canyon wind and road
REAL, DIMENSION(:), INTENT(IN)    :: PPS           ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PPA           ! pressure at the first atmospheric level
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS         ! surface exner function
REAL, DIMENSION(:), INTENT(IN)    :: PTA           ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA           ! specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD         ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PEXNA         ! exner function
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA         ! air density
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD       ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PRR           ! rain rate
REAL, DIMENSION(:), INTENT(IN)    :: PSR           ! snow rate
REAL, DIMENSION(:), INTENT(IN)    :: PH_TRAFFIC    ! anthropogenic sensible
REAL, DIMENSION(:), INTENT(IN)    :: PLE_TRAFFIC   ! anthropogenic latent
REAL, DIMENSION(:), INTENT(IN)    :: PH_INDUSTRY   ! anthropogenic sensible
REAL, DIMENSION(:), INTENT(IN)    :: PLE_INDUSTRY  ! anthropogenic latent
REAL, DIMENSION(:), INTENT(IN)    :: PZREF         ! reference height of the first
REAL, DIMENSION(:), INTENT(IN)    :: PUREF         ! reference height of the first
REAL,               INTENT(IN)    :: PTSTEP        ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PZ0_TOWN      ! town roughness length
REAL, DIMENSION(:), INTENT(IN)    :: PBLD          ! fraction of buildings
REAL, DIMENSION(:), INTENT(IN)    :: PGARDEN       ! fraction of GARDEN areas
REAL, DIMENSION(:), INTENT(IN)    :: PROAD         ! fraction of roads
REAL, DIMENSION(:), INTENT(IN)    :: PFRAC_GR      ! fraction of green roofs
REAL, DIMENSION(:), INTENT(IN)    :: PBLD_HEIGHT   ! buildings h
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_O_HOR   ! wall surf. / hor. surf.
REAL, DIMENSION(:), INTENT(IN)    :: PCAN_HW_RATIO ! canyon    h/W
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_O_GRND  ! wall surf. / (road + gard.  surf.)
REAL, DIMENSION(:), INTENT(INOUT) :: PDF_ROOF      ! snow-free    fraction on roofs
REAL, DIMENSION(:), INTENT(INOUT) :: PDN_ROOF      ! snow-covered fraction on roofs
REAL, DIMENSION(:), INTENT(INOUT) :: PDF_ROAD      ! snow-free    fraction on roads
REAL, DIMENSION(:), INTENT(INOUT) :: PDN_ROAD      ! snow-covered fraction on roads
REAL, DIMENSION(:), INTENT(OUT)   :: PQSAT_ROOF    ! hum at saturation over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PQSAT_ROAD    ! hum at saturation over road
REAL, DIMENSION(:), INTENT(OUT)   :: PDELT_ROOF    ! water fraction on roof
REAL, DIMENSION(:), INTENT(OUT)   :: PDELT_ROAD    ! water fraction on road
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_ROOF    ! roof emissivity
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_ROOF      ! heat capacity for roof layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_ROOF      ! thermal conductivity for roof layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_ROOF       ! depth of roof layers
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_ROAD    ! road emissivity
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_ROAD      ! heat capacity for road layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_ROAD      ! thermal conductivity for road layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_ROAD       ! depth of road layers
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_WALL    ! wall emissivity
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_WALL      ! heat capacity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_WALL      ! thermal conductivity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_WALL       ! depth of wall layers
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GARDEN    ! GARDEN area surf temp.
REAL, DIMENSION(:), INTENT(IN)    :: PG_GREENROOF_ROOF ! snow melt
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_ROOF     ! net radiation over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PH_ROOF      ! sensible heat flux over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_ROOF     ! latent heat flux over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PLEW_ROOF    ! latent heat flux over roof (snow)
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_ROOF  ! flux through the roof
REAL, DIMENSION(:), INTENT(OUT)   :: PRUNOFF_STRLROOF ! runoff over the ground
REAL, DIMENSION(:), INTENT(OUT)   :: PRUNOFF_ROOF ! runoff over roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_ROAD     ! net radiation over road
REAL, DIMENSION(:), INTENT(OUT)   :: PH_ROAD      ! sensible heat flux over road
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_ROAD     ! latent heat flux over road
REAL, DIMENSION(:), INTENT(OUT)   :: PLEW_ROAD    ! latent heat flux over road (snow)
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_ROAD  ! flux through the road
REAL, DIMENSION(:), INTENT(OUT)   :: PRUNOFF_ROAD ! runoff over the ground
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_WALL_A   ! net radiation over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PH_WALL_A    ! sensible heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WALL_A   ! latent heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_WALL_A! flux through the wall
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_WALL_B   ! net radiation over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PH_WALL_B    ! sensible heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WALL_B   ! latent heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_WALL_B! flux through the wall
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_BLT       ! net radiation over built surf 
REAL, DIMENSION(:), INTENT(OUT)   :: PH_BLT        ! sensible heat flux over built surf 
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_BLT       ! latent heat flux over built surf 
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_BLT    ! flux through the built surf 
REAL, DIMENSION(:), INTENT(IN)    :: PRN_GREENROOF     ! net radiation over greenroof
REAL, DIMENSION(:), INTENT(IN)    :: PH_GREENROOF      ! sensible heat flux over greenroof
REAL, DIMENSION(:), INTENT(IN)    :: PLE_GREENROOF     ! latent heat flux over greenroof
REAL, DIMENSION(:), INTENT(IN)    :: PGFLUX_GREENROOF  ! flux through the greenroof
REAL, DIMENSION(:), INTENT(IN)    :: PUW_GREENROOF     ! Momentum flux for greenroofs
REAL, DIMENSION(:), INTENT(IN)    :: PRUNOFF_GREENROOF ! runoff over green roofs
REAL, DIMENSION(:), INTENT(IN)    :: PDRAIN_GREENROOF  ! outlet drainage at base of green roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_STRLROOF      ! net radiation over structural roof
REAL, DIMENSION(:), INTENT(OUT)   :: PH_STRLROOF       ! sensible heat flux over structural roof
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_STRLROOF      ! latent heat flux over structural roof
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_STRLROOF   ! flux through the structural roof
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
REAL, DIMENSION(:), INTENT(OUT)   :: PUW_ROAD     ! Momentum flux for roads
REAL, DIMENSION(:), INTENT(OUT)   :: PUW_ROOF     ! Momentum flux for roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PDUWDU_ROAD  !
REAL, DIMENSION(:), INTENT(OUT)   :: PDUWDU_ROOF  !
REAL, DIMENSION(:), INTENT(OUT)   :: PUSTAR_TOWN  ! friciton velocity over town
REAL, DIMENSION(:), INTENT(OUT)   :: PCD          ! town averaged drag coefficient
REAL, DIMENSION(:), INTENT(OUT)   :: PCDN         ! town averaged neutral drag coefficient
REAL, DIMENSION(:), INTENT(OUT)   :: PCH_TOWN     ! town averaged heat transfer
REAL, DIMENSION(:), INTENT(OUT)   :: PRI_TOWN      ! town averaged Richardson number
REAL, DIMENSION(:), INTENT(OUT)   :: PRESA_TOWN    ! town aerodynamical resistance
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_TOWN     ! heat storage inside town
REAL, DIMENSION(:), INTENT(OUT)   :: PQF_TOWN      ! total anthropogenic heat
REAL, DIMENSION(:), INTENT(OUT)   :: PQF_BLD       ! anthropogenic heat flux of
REAL, DIMENSION(:), INTENT(OUT)   :: PFLX_BLD      ! heat flx from inside bld through its structure
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROOF      ! roof conductance
REAL, DIMENSION(:), INTENT(INOUT) :: PAC_ROAD      ! road conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_WALL      ! wall conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_TOP       ! top conductance
REAL, DIMENSION(:), INTENT(IN)    :: PAC_GARDEN    ! garden conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROOF_WAT  ! roof water conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROAD_WAT  ! roof water conductance
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROOF      ! absorbed solar rad by roof
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_SNOW_ROOF ! absorbed solar rad by snow on roof
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROAD      ! absorbed solar rad by road
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_SNOW_ROAD ! absorbed solar rad by snow on road
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WALL_A    ! absorbed solar rad by wall
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WALL_B    ! absorbed solar rad by wall
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_ROOF      ! absorbed IR rad by roof
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_SNOW_ROOF ! absorbed IR rad by snow on roof
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_ROAD      ! absorbed IR rad by road
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_SNOW_ROAD ! absorbed IR rad by snow on road
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_WALL_A    ! absorbed IR rad by wall
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_WALL_B    ! absorbed IR rad by wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_WB      ! LW contrib. wall A (orB) -> wall B (or A)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_R         ! LW contrib. wall       -> road
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WB_TO_R         ! LW contrib. wall       -> road
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_NR        ! LW contrib. wall       -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WB_TO_NR        ! LW contrib. wall       -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_WA         ! LW contrib. road       -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_WB         ! LW contrib. road       -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_WA         ! LW contrib. GARDEN     -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_WB         ! LW contrib. GARDEN     -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_WA        ! LW contrib. road(snow) -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_WB        ! LW contrib. road(snow) -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_WA         ! LW contrib. sky        -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_WB         ! LW contrib. sky        -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_R          ! LW contrib. sky        -> road
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_NR         ! LW contrib. sky        -> road(snow)
 CHARACTER(LEN=4), DIMENSION(:), INTENT(IN) :: HNATVENT
 CHARACTER(LEN=6), INTENT(IN)       :: HCOOL_COIL
REAL, DIMENSION(:), INTENT(IN)     :: PF_WATER_COND!fraction of evaporation for the condensers
 CHARACTER(LEN=6), INTENT(IN)       :: HHEAT_COIL
LOGICAL,            INTENT(IN)     :: OAUTOSIZE    ! Flag to activate autosize calculations
INTEGER,            INTENT(IN)     :: KDAY         ! Simulation day
REAL, DIMENSION(:), INTENT(INOUT)  :: PAUX_MAX     ! Auxiliar variable for autosize calcs
REAL, DIMENSION(:,:), INTENT(INOUT):: PT_FLOOR     ! Floor layers temperatures [K]
REAL, DIMENSION(:,:), INTENT(INOUT):: PT_MASS      ! Internal mass layers temperatures [K]
REAL, DIMENSION(:), INTENT(OUT)    :: PH_BLD_COOL  ! Sensible cooling energy demand  
REAL, DIMENSION(:), INTENT(OUT)   :: PT_BLD_COOL   ! Total cooling energy demand  
REAL, DIMENSION(:), INTENT(OUT)   :: PH_BLD_HEAT   ! Heating energy demand       
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_BLD_COOL  ! Latent cooling energy demand 
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_BLD_HEAT  ! Latent heating energy demand 
REAL, DIMENSION(:), INTENT(OUT)   :: PH_WASTE      ! Sensible waste heat from HVAC system
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WASTE     ! Latent waste heat from HVAC system
REAL, DIMENSION(:), INTENT(IN)    :: PF_WASTE_CAN  ! fraction of waste heat released into the canyon
REAL, DIMENSION(:), INTENT(OUT)   :: PHVAC_COOL    ! Energy consumption of the cooling system
REAL, DIMENSION(:), INTENT(OUT)  :: PHVAC_HEAT     ! Energy consumption of the heating system
REAL, DIMENSION(:), INTENT(IN)   :: PQIN          ! Internal heat gains [W m-2(floor)]
REAL, DIMENSION(:), INTENT(IN)   :: PQIN_FRAD     ! Radiant fraction of internal heat gains
REAL, DIMENSION(:), INTENT(IN)   :: PQIN_FLAT     ! Latent franction of internal heat gains
REAL, DIMENSION(:), INTENT(IN)   :: PGR           ! Glazing ratio
REAL, DIMENSION(:), INTENT(IN)   :: PEFF_HEAT     ! Efficiency of the heating system
REAL, DIMENSION(:), INTENT(IN)   :: PINF          ! Infiltration flow rate [AC/H]
REAL, DIMENSION(:), INTENT(IN)   :: PTCOOL_TARGET ! Cooling setpoint of HVAC system [K]
REAL, DIMENSION(:), INTENT(IN)   :: PTHEAT_TARGET ! Heating setpoint of HVAC system [K]
REAL, DIMENSION(:), INTENT(IN)   :: PHR_TARGET    ! Relative humidity setpoint
REAL, DIMENSION(:), INTENT(INOUT):: PT_WIN2       ! Indoor window temperature [K]
REAL, DIMENSION(:), INTENT(INOUT):: PQI_BLD       ! Indoor air specific humidity [kg kg-1]
REAL, DIMENSION(:), INTENT(IN)   :: PV_VENT       ! Ventilation flow rate [AC/H]
REAL, DIMENSION(:), INTENT(IN)   :: PCAP_SYS_HEAT ! Capacity of the heating system 
REAL, DIMENSION(:), INTENT(INOUT):: PCAP_SYS_RAT  ! Rated capacity of the cooling system
REAL, DIMENSION(:), INTENT(IN)   :: PT_ADP        ! Apparatus dewpoint temperature of the
REAL, DIMENSION(:), INTENT(INOUT):: PM_SYS_RAT    ! Rated HVAC mass flow rate 
REAL, DIMENSION(:), INTENT(IN)   :: PCOP_RAT      ! Rated COP of the cooling system
REAL, DIMENSION(:), INTENT(OUT)  :: PCAP_SYS      ! Actual capacity of the cooling system
REAL, DIMENSION(:), INTENT(OUT)  :: PM_SYS        ! Actual HVAC mass flow rate 
REAL, DIMENSION(:), INTENT(OUT)  :: PCOP          ! COP of the cooling system
REAL, DIMENSION(:), INTENT(OUT)  :: PQ_SYS        ! Supply air specific humidity [kg kg-1]
REAL, DIMENSION(:), INTENT(OUT)  :: PT_SYS        ! Supply air temperature [K]
REAL, DIMENSION(:), INTENT(IN)   :: PTR_SW_WIN    ! Solar radiation transmitted throught
REAL, DIMENSION(:), INTENT(OUT)  :: PFAN_POWER    ! HVAC fan power
REAL, DIMENSION(:,:),INTENT(IN)  :: PHC_FLOOR     ! heat capacity for road layers
REAL, DIMENSION(:,:),INTENT(IN)  :: PTC_FLOOR     ! thermal conductivity for 
REAL, DIMENSION(:,:),INTENT(IN)  :: PD_FLOOR      ! depth of road layers
REAL, DIMENSION(:), INTENT(IN)   :: PLW_WA_TO_WIN ! Radiative heat trasfer coeff wall-window 
REAL, DIMENSION(:), INTENT(IN)   :: PLW_WB_TO_WIN ! Radiative heat trasfer coeff wall-window 
REAL, DIMENSION(:), INTENT(IN)   :: PLW_G_TO_WIN  ! Radiative heat trasfer coeff garden-window 
REAL, DIMENSION(:), INTENT(IN)   :: PLW_R_TO_WIN  ! Radiative heat trasfer coeff road-window 
REAL, DIMENSION(:), INTENT(IN)   :: PLW_S_TO_WIN ! Radiative heat trasfer coeff window-sky 
REAL, DIMENSION(:), INTENT(IN)   :: PLW_WIN_TO_WA! Radiative heat trasfer coeff window-wall
REAL, DIMENSION(:), INTENT(IN)   :: PLW_WIN_TO_WB! Radiative heat trasfer coeff window-wall
REAL, DIMENSION(:), INTENT(IN)   :: PLW_WIN_TO_R ! Radiative heat trasfer coeff window-road 
REAL, DIMENSION(:), INTENT(IN)   :: PLW_NR_TO_WIN! Radiative heat trasfer coeff road(snow)-win 
REAL, DIMENSION(:), INTENT(IN)   :: PLW_WIN_TO_NR! Radiative heat trasfer coeff win-road(snow) 
REAL, DIMENSION(:), INTENT(INOUT):: PT_WIN1       ! outdoor window temperature [K]
REAL, DIMENSION(:), INTENT(IN)   :: PABS_SW_WIN   ! window absorbed shortwave radiation [W m-2] 
REAL, DIMENSION(:), INTENT(OUT)  :: PABS_LW_WIN   ! absorbed infrared rad. [W m-2]
REAL, DIMENSION(:), INTENT(IN)   :: PUGG_WIN      ! window glass-to-glass U-factro [W m-2 K-1]
REAL, DIMENSION(:), INTENT(OUT) :: PEMIT_LW_ROAD ! LW fluxes emitted by road (W/m2 surf road)
REAL, DIMENSION(:), INTENT(OUT) :: PEMIT_LW_FAC  ! LW fluxes emitted by wall (W/m2 surf wall)
REAL, DIMENSION(:), INTENT(OUT) :: PT_RAD_IND    ! Indoor mean radiant temperature [K]
REAL, DIMENSION(:), INTENT(OUT) :: PHU_BLD       ! Indoor relative humidity 0 < (-) < 1
REAL,                INTENT(IN)  :: PTIME        ! current time since midnight (UTC, s)
REAL, DIMENSION(:), INTENT(IN)  :: PE_SHADING    !energy not ref., nor absorbed, nor
LOGICAL, DIMENSION(:),INTENT(INOUT) :: ONATVENT_NIGHT ! has natural ventilation
REAL, DIMENSION(:)  , INTENT(IN)  :: PN_FLOOR        ! Number of floors     
REAL, DIMENSION(:)  , INTENT(IN)  :: PWALL_O_BLD     ! Wall area [m2_wall/m2_bld]
REAL, DIMENSION(:)  , INTENT(IN)  :: PGLAZ_O_BLD     ! Window area [m2_win/m2_bld]
REAL, DIMENSION(:)  , INTENT(IN)  :: PMASS_O_BLD     ! Mass area [m2_mass/m2_bld]
REAL, DIMENSION(:)  , INTENT(IN)  :: PFLOOR_HW_RATIO ! H/W ratio of 1 floor level
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_FLOOR_MASS   ! View factor floor-mass
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_FLOOR_WALL   ! View factor floor-wall
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_FLOOR_WIN    ! View factor floor-window
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_FLOOR_ROOF   ! View factor floor-roof
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_WALL_FLOOR   ! View factor wall-floor
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_WALL_MASS    ! View factor wall-mass
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_WALL_WIN     ! View factor wall-win
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_WIN_FLOOR    ! View factor win-floor
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_WIN_MASS     ! View factor win-mass
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_WIN_WALL     ! View factor win-wall
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_MASS_FLOOR   ! View factor mass-floor
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_MASS_WALL    ! View factor mass-wall
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_MASS_WIN     ! View factor mass-window
LOGICAL,              INTENT(IN)  :: OCANOPY         ! is canopy active ?
 CHARACTER(LEN=5), INTENT(IN)      :: HCH_BEM         ! TEB option for outside building coef
REAL, DIMENSION(:)  , INTENT(IN)  :: PROUGH_ROOF     ! roof roughness coef
REAL, DIMENSION(:)  , INTENT(IN)  :: PROUGH_WALL     ! wall roughness coef
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_WIN_WIN      ! indoor win to win view factor
LOGICAL,              INTENT(IN)  :: OPAR_RD_IRRIG  ! True if there is potentially irrigation
REAL, DIMENSION(:),   INTENT(IN)  :: PRD_START_MONTH! First month of irrigation
REAL, DIMENSION(:),   INTENT(IN)  :: PRD_END_MONTH  ! Last month of irrigation
REAL, DIMENSION(:),   INTENT(IN)  :: PRD_START_HOUR ! First hour of irrigation (included)
REAL, DIMENSION(:),   INTENT(IN)  :: PRD_END_HOUR   ! Last hour of irrigation  (excluded)
REAL, DIMENSION(:),   INTENT(IN)  :: PRD_24H_IRRIG  ! diurnal averaged irrigation flux (kg/m2/s)
REAL, DIMENSION(:),   INTENT(OUT) :: PIRRIG_ROAD    ! man-made watering for road (kg/m2/s)
END SUBROUTINE TEB
END INTERFACE
END MODULE MODI_TEB
