!auto_modi:spll_bem.D
MODULE MODI_BEM
INTERFACE
        SUBROUTINE BEM(PTSTEP, PSUNTIME, HCOOL_COIL, HHEAT_COIL,              &
                OAUTOSIZE, KDAY,  HNATVENT,                                   &
                PPS, PRHOA, PT_CANYON, PQ_CANYON, PU_CANYON,                  &
                PT_ROOF, PT_WALL_A, PT_WALL_B, PBLD, PBLD_HEIGHT, PWALL_O_HOR,&
                PHC_FLOOR, PTC_FLOOR, PD_FLOOR,                               &
                PT_WIN2, PGR, PQIN, PQIN_FRAD, PQIN_FLAT,                     &
                PEFF_HEAT, PINF, PTCOOL_TARGET, PTHEAT_TARGET, PHR_TARGET,    &
                PF_WATER_COND, PV_VENT, PCAP_SYS_HEAT, PT_ADP, PCOP_RAT,      &
                PAUX_MAX, PT_FLOOR, PT_MASS, PTI_BLD, PQI_BLD,                &
                PCAP_SYS_RAT, PM_SYS_RAT, ONATVENT_NIGHT,                     &
                PH_BLD_COOL, PH_BLD_HEAT, PLE_BLD_COOL, PLE_BLD_HEAT,         &
                PT_BLD_COOL, PHVAC_COOL, PT_SYS, PQ_SYS, PH_WASTE, PLE_WASTE, &
                PFAN_POWER, PHVAC_HEAT, PM_SYS, PCOP, PCAP_SYS, PHU_BLD,      &
                PTR_SW_WIN, PT_RAD_IND, PFLX_BLD_FLOOR, PFLX_BLD_MASS,        &
                PF_FLOOR_MASS, PF_FLOOR_WALL, PF_FLOOR_WIN,                   &
                PRADHT_IN, PN_FLOOR, PWALL_O_BLD, PGLAZ_O_BLD, PMASS_O_BLD,   &
                PFLOOR_HW_RATIO, PF_FLOOR_ROOF, PF_MASS_FLOOR, PF_MASS_WALL,  &
                PF_MASS_WIN, PRAD_ROOF_MASS, PRAD_ROOF_FLOOR, PRAD_WALL_MASS, &
                PRAD_WALL_FLOOR, PRAD_WIN_MASS, PRAD_WIN_FLOOR,               &
                PCONV_ROOF_BLD, PCONV_WALL_BLD, PCONV_WIN_BLD, &
                PLOAD_IN_FLOOR, PLOAD_IN_MASS                )
IMPLICIT NONE
REAL,                INTENT(IN)   :: PTSTEP        ! Time step
REAL, DIMENSION(:),  INTENT(IN)   :: PSUNTIME       ! current solar time since midnight (solar time, s)
 CHARACTER(LEN=6),    INTENT(IN)   :: HCOOL_COIL    ! type of cooling system IDEAL/DX_COOL
 CHARACTER(LEN=6),    INTENT(IN)   :: HHEAT_COIL    ! type of heating system IDEAL/FIN_CAP
LOGICAL,             INTENT(IN)   :: OAUTOSIZE     ! Flag to activate autosize calculations
INTEGER,             INTENT(IN)   :: KDAY          ! Simulation day
 CHARACTER(LEN=4), DIMENSION(:), INTENT(IN) :: HNATVENT
REAL, DIMENSION(:),   INTENT(IN)  :: PPS          ! Canyon air pressure [Pa]
REAL, DIMENSION(:),   INTENT(IN)  :: PRHOA        ! Air density at the lowest level [kg m-3]
REAL, DIMENSION(:),   INTENT(IN)  :: PT_CANYON    ! Canyon air temperature [K]
REAL, DIMENSION(:),   INTENT(IN)  :: PQ_CANYON    ! Canyon air specific humidity [kg kg-1]
REAL, DIMENSION(:),   INTENT(IN)  :: PU_CANYON    ! Canyon wind speed (m s-1)
REAL, DIMENSION(:,:), INTENT(IN)  :: PT_ROOF      ! Roof layers temperatures [K]
REAL, DIMENSION(:,:), INTENT(IN)  :: PT_WALL_A    ! Wall A layers temperatures [K]
REAL, DIMENSION(:,:), INTENT(IN)  :: PT_WALL_B    ! Wall B layers temperatures [K]
REAL, DIMENSION(:),   INTENT(IN)  :: PBLD         ! Urban horizontal building density
REAL, DIMENSION(:),   INTENT(IN)  :: PBLD_HEIGHT  ! Average building height [m]
REAL, DIMENSION(:),   INTENT(IN)  :: PWALL_O_HOR  ! Wall to horizontal surface ratio
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_FLOOR    ! heat capacity for road layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_FLOOR    ! thermal conductivity for 
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_FLOOR     ! depth of road layers
REAL, DIMENSION(:),   INTENT(IN)   :: PT_WIN2       ! Indoor window temperature [K]
REAL, DIMENSION(:),   INTENT(IN)   :: PGR           ! Glazing ratio
REAL, DIMENSION(:),   INTENT(IN)   :: PQIN          ! Internal heat gains [W m-2(floor)]
REAL, DIMENSION(:),   INTENT(IN)   :: PQIN_FRAD     ! Radiant fraction of internal heat gains
REAL, DIMENSION(:),   INTENT(IN)   :: PQIN_FLAT     ! Latent franction of internal heat gains
REAL, DIMENSION(:),   INTENT(IN)   :: PEFF_HEAT     ! Efficiency of the heating system
REAL, DIMENSION(:),   INTENT(IN)   :: PINF          ! Infiltration flow rate [AC/H]
REAL, DIMENSION(:),   INTENT(IN)   :: PTCOOL_TARGET ! Cooling setpoint of HVAC system [K]
REAL, DIMENSION(:),   INTENT(IN)   :: PTHEAT_TARGET ! Heating setpoint of HVAC system [K]
REAL, DIMENSION(:),   INTENT(IN)   :: PHR_TARGET    ! Relative humidity setpoint
REAL,DIMENSION(:),    INTENT(IN)   :: PF_WATER_COND !fraction of evaporation for the condensers
REAL, DIMENSION(:),   INTENT(IN)   :: PV_VENT       ! Ventilation flow rate [AC/H]
REAL, DIMENSION(:),   INTENT(IN)   :: PCAP_SYS_HEAT ! Capacity of the heating system 
REAL, DIMENSION(:),   INTENT(IN)   :: PT_ADP        ! Apparatus dewpoint temperature of the
REAL, DIMENSION(:),   INTENT(IN)   :: PCOP_RAT      ! Rated COP of the cooling system
REAL, DIMENSION(:),   INTENT(INOUT):: PAUX_MAX      ! Auxiliar variable for autosize calcs
REAL, DIMENSION(:,:), INTENT(INOUT):: PT_FLOOR      ! Floor layers temperatures [K]
REAL, DIMENSION(:,:), INTENT(INOUT):: PT_MASS       ! Internal mass layers temperatures [K]
REAL, DIMENSION(:),   INTENT(INOUT):: PTI_BLD       ! Indoor air temperature [K]
REAL, DIMENSION(:),   INTENT(INOUT):: PQI_BLD       ! Indoor air specific humidity [kg kg-1]
REAL, DIMENSION(:),   INTENT(INOUT):: PCAP_SYS_RAT  ! Rated capacity of the cooling system
REAL, DIMENSION(:),   INTENT(INOUT):: PM_SYS_RAT    ! Rated HVAC mass flow rate 
LOGICAL, DIMENSION(:),INTENT(INOUT):: ONATVENT_NIGHT ! has natural ventilation
REAL, DIMENSION(:),   INTENT(OUT)  :: PH_BLD_COOL   ! Sensible cooling energy demand  
REAL, DIMENSION(:),   INTENT(OUT)  :: PH_BLD_HEAT   ! Heating energy demand       
REAL, DIMENSION(:),   INTENT(OUT)  :: PLE_BLD_COOL  ! Latent cooling energy demand 
REAL, DIMENSION(:),   INTENT(OUT)  :: PLE_BLD_HEAT  ! Latent heating energy demand 
REAL, DIMENSION(:),   INTENT(OUT)  :: PT_BLD_COOL   ! Total cooling energy demand  
REAL, DIMENSION(:),   INTENT(OUT)  :: PHVAC_COOL    ! Energy consumption of the cooling system
REAL, DIMENSION(:),   INTENT(OUT)  :: PT_SYS        ! Supply air temperature [K]                         
REAL, DIMENSION(:),   INTENT(OUT)  :: PQ_SYS        ! Supply air specific humidity [kg kg-1]
REAL, DIMENSION(:),   INTENT(OUT)  :: PH_WASTE      ! Sensible waste heat from HVAC system
REAL, DIMENSION(:),   INTENT(OUT)  :: PLE_WASTE     ! Latent waste heat from HVAC system
REAL, DIMENSION(:),   INTENT(OUT)  :: PFAN_POWER    ! HVAC fan power
REAL, DIMENSION(:),   INTENT(OUT)  :: PHVAC_HEAT    ! Energy consumption of the heating system
REAL, DIMENSION(:),   INTENT(OUT)  :: PM_SYS        ! Actual HVAC mass flow rate 
REAL, DIMENSION(:),   INTENT(OUT)  :: PCOP          ! COP of the cooling system
REAL, DIMENSION(:),   INTENT(OUT)  :: PCAP_SYS      ! Actual capacity of the cooling system
REAL, DIMENSION(:),   INTENT(OUT)  :: PHU_BLD       ! Indoor relative humidity 0 < (-) < 1
REAL, DIMENSION(:),   INTENT(IN)   :: PTR_SW_WIN    ! Solar radiation transmitted throught
REAL, DIMENSION(:),   INTENT(IN)  :: PT_RAD_IND    ! Indoor mean radiant temperature [K]
REAL, DIMENSION(:),   INTENT(OUT)  :: PFLX_BLD_FLOOR! Heat flux from indoor air to floor 
REAL, DIMENSION(:),   INTENT(OUT)  :: PFLX_BLD_MASS ! Heat flux from indoor air to mass 
REAL, DIMENSION(:),   INTENT(IN)  :: PF_FLOOR_MASS ! View factor floor-mass
REAL, DIMENSION(:),   INTENT(IN)  :: PF_FLOOR_WALL ! View factor floor-wall
REAL, DIMENSION(:),   INTENT(IN)  :: PF_FLOOR_WIN  ! View factor floor-window
REAL, DIMENSION(:),   INTENT(IN) :: PRADHT_IN     ! Indoor radiant heat transfer coefficient
REAL, DIMENSION(:)  , INTENT(IN)  :: PN_FLOOR        ! Number of floors     
REAL, DIMENSION(:)  , INTENT(IN)  :: PWALL_O_BLD         ! Wall area [m2_wall/m2_bld]
REAL, DIMENSION(:)  , INTENT(IN)  :: PGLAZ_O_BLD          ! Window area [m2_win/m2_bld]
REAL, DIMENSION(:)  , INTENT(IN)  :: PMASS_O_BLD         ! Mass area [m2_mass/m2_bld]
REAL, DIMENSION(:)  , INTENT(IN)  :: PFLOOR_HW_RATIO ! H/W ratio of 1 floor level
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_FLOOR_ROOF   ! View factor floor-roof
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_MASS_FLOOR   ! View factor mass-floor
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_MASS_WALL    ! View factor mass-wall
REAL, DIMENSION(:)  , INTENT(IN)  :: PF_MASS_WIN     ! View factor mass-window
REAL, DIMENSION(:)  , INTENT(IN)  :: PRAD_ROOF_MASS  ! Rad. fluxes between roof and mass
REAL, DIMENSION(:)  , INTENT(IN)  :: PRAD_ROOF_FLOOR ! Rad. fluxes between roof and floor
REAL, DIMENSION(:)  , INTENT(IN)  :: PRAD_WALL_MASS  ! Rad. fluxes between wall and mass
REAL, DIMENSION(:)  , INTENT(IN)  :: PRAD_WALL_FLOOR ! Rad. fluxes between wall and floor
REAL, DIMENSION(:)  , INTENT(IN)  :: PRAD_WIN_MASS   ! Rad. fluxes between wind. and mass
REAL, DIMENSION(:)  , INTENT(IN)  :: PRAD_WIN_FLOOR  ! Rad. fluxes between wind. and floor
REAL, DIMENSION(:)  , INTENT(IN)  :: PCONV_ROOF_BLD  ! Conv. fluxes between roof and indoor air
REAL, DIMENSION(:)  , INTENT(IN)  :: PCONV_WALL_BLD  ! Conv. fluxes between wall and indoor air
REAL, DIMENSION(:)  , INTENT(IN)  :: PCONV_WIN_BLD   ! Conv. fluxes between wind. and indoor air
REAL, DIMENSION(:)  , INTENT(IN)  :: PLOAD_IN_FLOOR  ! solar + int heat gain on floor W/m2 [floor]
REAL, DIMENSION(:)  , INTENT(IN)  :: PLOAD_IN_MASS   ! solar + int heat gain on floor W/m2 [mass]
END SUBROUTINE BEM
END INTERFACE
END MODULE MODI_BEM
