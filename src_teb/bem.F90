!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!   ##########################################################################
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
!   ##########################################################################
!
!!****  *BEM*
!!
!!    PURPOSE
!!    -------
!
!     Computes the temperature and humidity evolution of indoor air, 
!     building energy demand, HVAC energy consumption, 
!     waste heat from HVAC systems, and heat fluxes from indoor to building surfaces.
!
!
!!**  METHOD
!     ------
!
!              NOMENCLATURE: bld  - refers to building plant area; 
!                            floor- refers to building plant area multiplied 
!                                   by the number of floors;
!                            wall - refers to wall area (excluding windows).
!                            win  - refers to window area. 
!                            mass - refers to internal mass area. 
!
!
!        solar radiation transmitted through windows
!        *******************************************
!
!     Qsol_tr_win = Qsol_facade * tr_win * GR 
!
!
!        indoor wall conv/rad heat transfer coefficients
!        ***********************************************
!
!     The calculation of CHTC accounts for favorable or unfavorable convection 
!     depending on the relative position between the hot layer and cold layer
!
! 
!        building energy demand
!        **********************
!
!     Calculation of the cooling and heating, sensible and latent building energy demand.
!     The sensible demand includes the convective heat transfer from indoor surfaces, the 
!     convective fraction of internal heat gains, and sensible infiltration/ventilation heat
!     gains. The latent demand includes the latent fraction of internal heat gains and latent
!     infiltration/ventilation heat gains.  
!
!        surface areas and volummes (referred to m2_bld)
!        ***********************************************
!
!     Awall   =  WALL_O_HOR * (1 - GR) / BLD [m2_wall/m2_bld]
!     Awin    =  WALL_O_HOR * GR / BLD       [m2_win/m2_bld]   
!     Amass   =  2 * N_FLOOR                  [m2_mass/m2_bld]  
!     N_FLOOR  =  BLD_HEIGHT / FLOOR_HEIGHT   [#]
!     Aroof   =  1                           [m2_roof/m2_bld]  
!     Afloor  =  1                           [m2_floor/m2_bld]   
!     Vol_air =  BLD_HEIGHT                  [m3_bld/m2_bld]
!
!
!        evolution of the internal temperature
!        *************************************
!

!                                  dTin  
!     Vol_air * ro_air * cp_air * ---- = h_wall * Awall * (Twall - Tin)
!                                   dt    + h_roof * Aroof * (Troof -Tin)
!                                         + h_floor * Afloor *(Tfloor - Tin)
!                                         + h_mass * Amass * (Tmass - Tin)  
!                                         + h_win * Awin * (Twin - Tin)
!                                         + Qig * (1 - fig_rad) * (1-fig_lat)
!                                         + Vinf * ro_air * cp_air * (Tout - Tin) 
!                                         + Vsys * ro_air * cp_air * (Tsys - Tin) 
!
!
!        evolution of the internal specific humidity
!        *******************************************
!
!                                  dQin  
!      Vol_air * ro_air * lv_air * ---- = Qig * fig_lat
!                                   dt    + Vinf * ro_air * lv_air * (Qout - Qin) 
!                                         + Vsys * ro_air * lv_air * (Qsys - Qin) 
!
!
!        heat fluxes from indoor to surfaces
!        ***********************************
!
!      Qin_wall  = h_wall  * (Tin - Twall)  [W/m2_wall]
!      Qin_roof  = h_roof  * (Tin - Troof)  [W/m2_roof] 
!      Qin_floor = h_floor * (Tin - Tfloor) [W/m2_floor] 
!      Qin_mass  = h_wall  * (Tin - Tmass)  
!                + Qig * fig_rad * (1-fig_lat)/ 2  
!                + Qsol_tr_win              [W/m2_mass]
!
!
!        energy consumption and waste heat from cooling system
!        *****************************************************
!
!      Qhvac  = Qbld / COP
!      Qwaste = Qbld + Qhvac
!
!
!        energy consumption and waste heat from heating system
!        *****************************************************
!
!      Qhvac  = Qbld / Eff
!      Qwaste = Qhvac - Qbld
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      B. Bueno           * Meteo-France *
!!
!!!    MODIFICATIONS
!!    -------------
!!     Original 2010
!!     G. Pigeon nov. 2011: inclusion floor/mass budget inside
!!                          add automatic/manual ventilation
!!                          conserve exchanges with the different surfaces inside 1 time step
!!    G. Pigeon sept. 2012: use of TARP/DOE coef for indoor convection
!!                          use of both PT_WALL_A and PT_WALL_B for calculations
!!                          the internal mass depth is 1/2 of the floor depth
!!                          add the option of no atmospheric heat releases by HVAC system (PF_WATER_COND < 0)
!!    G. Pigeon oct. 2012:  use indoor air density + new solar heat gain distribution
!!    V. Masson May  2013   implicitation of internal building temperature evolution
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XCPD,XSTEFAN,XLVTT,XG, XRV, XRD
 
USE MODE_THERMOS
USE MODE_PSYCHRO
USE MODI_DX_AIR_COOLING_COIL_CV
USE MODI_FLOOR_LAYER_E_BUDGET
USE MODI_MASS_LAYER_E_BUDGET
USE MODE_CONV_DOE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    Declarations of arguments
!
!
REAL,                INTENT(IN)   :: PTSTEP        ! Time step
REAL, DIMENSION(:),  INTENT(IN)   :: PSUNTIME       ! current solar time since midnight (solar time, s)
 CHARACTER(LEN=6),    INTENT(IN)   :: HCOOL_COIL    ! type of cooling system IDEAL/DX_COOL
 CHARACTER(LEN=6),    INTENT(IN)   :: HHEAT_COIL    ! type of heating system IDEAL/FIN_CAP
LOGICAL,             INTENT(IN)   :: OAUTOSIZE     ! Flag to activate autosize calculations
INTEGER,             INTENT(IN)   :: KDAY          ! Simulation day
!
 CHARACTER(LEN=4), DIMENSION(:), INTENT(IN) :: HNATVENT
!
REAL, DIMENSION(:),   INTENT(IN)  :: PPS          ! Canyon air pressure [Pa]
REAL, DIMENSION(:),   INTENT(IN)  :: PRHOA        ! Air density at the lowest level [kg m-3]
REAL, DIMENSION(:),   INTENT(IN)  :: PT_CANYON    ! Canyon air temperature [K]
REAL, DIMENSION(:),   INTENT(IN)  :: PQ_CANYON    ! Canyon air specific humidity [kg kg-1]
REAL, DIMENSION(:),   INTENT(IN)  :: PU_CANYON    ! Canyon wind speed (m s-1)
REAL, DIMENSION(:,:), INTENT(IN)  :: PT_ROOF      ! Roof layers temperatures [K]
REAL, DIMENSION(:,:), INTENT(IN)  :: PT_WALL_A    ! Wall A layers temperatures [K]
REAL, DIMENSION(:,:), INTENT(IN)  :: PT_WALL_B    ! Wall B layers temperatures [K]
!
REAL, DIMENSION(:),   INTENT(IN)  :: PBLD         ! Urban horizontal building density
REAL, DIMENSION(:),   INTENT(IN)  :: PBLD_HEIGHT  ! Average building height [m]
REAL, DIMENSION(:),   INTENT(IN)  :: PWALL_O_HOR  ! Wall to horizontal surface ratio
!
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_FLOOR    ! heat capacity for road layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_FLOOR    ! thermal conductivity for 
                                                  ! road layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_FLOOR     ! depth of road layers
!
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
                                                    ! [W m-2(bld)]
REAL, DIMENSION(:),   INTENT(IN)   :: PT_ADP        ! Apparatus dewpoint temperature of the
                                                    ! cooling coil [K]          
REAL, DIMENSION(:),   INTENT(IN)   :: PCOP_RAT      ! Rated COP of the cooling system
!
REAL, DIMENSION(:),   INTENT(INOUT):: PAUX_MAX      ! Auxiliar variable for autosize calcs
REAL, DIMENSION(:,:), INTENT(INOUT):: PT_FLOOR      ! Floor layers temperatures [K]
REAL, DIMENSION(:,:), INTENT(INOUT):: PT_MASS       ! Internal mass layers temperatures [K]
REAL, DIMENSION(:),   INTENT(INOUT):: PTI_BLD       ! Indoor air temperature [K]
REAL, DIMENSION(:),   INTENT(INOUT):: PQI_BLD       ! Indoor air specific humidity [kg kg-1]
                 
REAL, DIMENSION(:),   INTENT(INOUT):: PCAP_SYS_RAT  ! Rated capacity of the cooling system
                                                    ! [W m-2(bld)]
REAL, DIMENSION(:),   INTENT(INOUT):: PM_SYS_RAT    ! Rated HVAC mass flow rate 
                                                    ! [kg s-1 m-2(bld)]
LOGICAL, DIMENSION(:),INTENT(INOUT):: ONATVENT_NIGHT ! has natural ventilation
                                                     ! been necessary/possible during the night
!
REAL, DIMENSION(:),   INTENT(OUT)  :: PH_BLD_COOL   ! Sensible cooling energy demand  
                                                    ! of the building [W m-2(bld)]
REAL, DIMENSION(:),   INTENT(OUT)  :: PH_BLD_HEAT   ! Heating energy demand       
                                                    ! of the building [W m-2(bld)]
REAL, DIMENSION(:),   INTENT(OUT)  :: PLE_BLD_COOL  ! Latent cooling energy demand 
                                                    ! of the building [W m-2(bld)]
REAL, DIMENSION(:),   INTENT(OUT)  :: PLE_BLD_HEAT  ! Latent heating energy demand 
                                                    ! of the building [W m-2(bld)]
REAL, DIMENSION(:),   INTENT(OUT)  :: PT_BLD_COOL   ! Total cooling energy demand  
                                                    ! of the building [W m-2(bld)]
REAL, DIMENSION(:),   INTENT(OUT)  :: PHVAC_COOL    ! Energy consumption of the cooling system
                                                    ! [W m-2(bld)]
REAL, DIMENSION(:),   INTENT(OUT)  :: PT_SYS        ! Supply air temperature [K]                         
REAL, DIMENSION(:),   INTENT(OUT)  :: PQ_SYS        ! Supply air specific humidity [kg kg-1]
REAL, DIMENSION(:),   INTENT(OUT)  :: PH_WASTE      ! Sensible waste heat from HVAC system
                                                    ! [W m-2(bld)]
REAL, DIMENSION(:),   INTENT(OUT)  :: PLE_WASTE     ! Latent waste heat from HVAC system
                                                    ! [W m-2(bld)]
REAL, DIMENSION(:),   INTENT(OUT)  :: PFAN_POWER    ! HVAC fan power
REAL, DIMENSION(:),   INTENT(OUT)  :: PHVAC_HEAT    ! Energy consumption of the heating system
                                                    ! [W m-2(bld)]
!
REAL, DIMENSION(:),   INTENT(OUT)  :: PM_SYS        ! Actual HVAC mass flow rate 
                                                    ! [kg s-1 m-2(bld)]
REAL, DIMENSION(:),   INTENT(OUT)  :: PCOP          ! COP of the cooling system
REAL, DIMENSION(:),   INTENT(OUT)  :: PCAP_SYS      ! Actual capacity of the cooling system
                                                    ! [W m-2(bld)] 
!
REAL, DIMENSION(:),   INTENT(OUT)  :: PHU_BLD       ! Indoor relative humidity 0 < (-) < 1
REAL, DIMENSION(:),   INTENT(IN)   :: PTR_SW_WIN    ! Solar radiation transmitted throught
                                                    ! windows [W m-2(bld)]
REAL, DIMENSION(:),   INTENT(IN)  :: PT_RAD_IND    ! Indoor mean radiant temperature [K]
!
REAL, DIMENSION(:),   INTENT(OUT)  :: PFLX_BLD_FLOOR! Heat flux from indoor air to floor 
                                                    ! [W m-2(bld)]
REAL, DIMENSION(:),   INTENT(OUT)  :: PFLX_BLD_MASS ! Heat flux from indoor air to mass 
                                                    ! [W m-2(bld)]
REAL, DIMENSION(:),   INTENT(IN)  :: PF_FLOOR_MASS ! View factor floor-mass
REAL, DIMENSION(:),   INTENT(IN)  :: PF_FLOOR_WALL ! View factor floor-wall
REAL, DIMENSION(:),   INTENT(IN)  :: PF_FLOOR_WIN  ! View factor floor-window
REAL, DIMENSION(:),   INTENT(IN) :: PRADHT_IN     ! Indoor radiant heat transfer coefficient
                                                    ! [W K-1 m-2]
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
!
!*      0.2    Declarations of local variables 
!
INTEGER                        :: IROOF        ! Number of roof layers
INTEGER                        :: IWALL        ! Number of wall layers
!REAL                           :: ZTCOMF_MAX   ! Maximum comfort temperature for nat.vent [K]
!
REAL, DIMENSION(SIZE(PTI_BLD)) :: ZFAN_AP      ! Fan design pressure increase [Pa]
REAL, DIMENSION(SIZE(PTI_BLD)) :: ZFAN_EFF     ! Fan total efficiency
!
LOGICAL, DIMENSION(SIZE(PTI_BLD)):: GSCHED     ! Day-night schedule flag 
                                               ! *to be transported to inputs*
!
REAL, DIMENSION(SIZE(PTI_BLD)) :: ZF_NIGHT     ! Reduction factor of int.gains at night
REAL, DIMENSION(SIZE(PTI_BLD)) :: ZF_DAY       ! Amplification factor of int.gains at daytime
!
REAL, DIMENSION(SIZE(PTI_BLD)):: ZAC_IN_MASS_COOL, ZAC_IN_FLOOR_COOL, &
                                 ZAC_IN_ROOF_COOL, ZAC_IN_WALL_A_COOL, &
                                 ZAC_IN_WALL_B_COOL, ZAC_IN_WIN_COOL   
REAL, DIMENSION(SIZE(PTI_BLD)):: ZAC_IN_MASS_HEAT, ZAC_IN_FLOOR_HEAT, &
                                 ZAC_IN_ROOF_HEAT, ZAC_IN_WALL_A_HEAT, &
                                 ZAC_IN_WALL_B_HEAT, ZAC_IN_WIN_HEAT   
!
REAL, DIMENSION(SIZE(PTI_BLD)):: ZQIN          ! Internal heat gains [W m-2(bld)]
!
REAL, DIMENSION(SIZE(PTI_BLD)):: ZV_VENT       ! Ventilation flow rate [m3 s-1 m-2(bld)]
REAL, DIMENSION(SIZE(PTI_BLD)):: ZINF          ! Infiltration flow rate [m3 s-1 m-2(bld)]
!
LOGICAL, DIMENSION(SIZE(PTI_BLD)):: GNAT_VENT  ! Is Natural ventilation active ? 
REAL, DIMENSION(SIZE(PTI_BLD)):: ZNAT_VENT     ! Nat.vent airflow rate [m3 s-1 m-2(bld)]
REAL,DIMENSION(SIZE(PTI_BLD)) :: ZTI_BLD       ! Indoor air temperature at time step t + dt [K]
REAL,DIMENSION(SIZE(PTI_BLD)) :: ZTI_BLD_OPEN  ! Indoor air temperature if windows opened
REAL,DIMENSION(SIZE(PTI_BLD)) :: ZTI_BLD_CLOSED! Indoor air temperature if windows closed
!
REAL, DIMENSION(SIZE(PTI_BLD)):: ZQCOOL_TARGET ! Specific humidity cooling setpoing [kg kg-1]
REAL, DIMENSION(SIZE(PTI_BLD)):: ZQHEAT_TARGET ! Specific humidity heating setpoing [kg kg-1]
!
REAL, DIMENSION(SIZE(PTI_BLD)):: ZSHR          ! Rated sensible heat rate
REAL, DIMENSION(SIZE(PTI_BLD)):: ZM_SYS_RAT    ! Auxiliar mass flow rate [kg s-1 m-2(bld)]
!
REAL, DIMENSION(SIZE(PTI_BLD)):: ZXMIX         ! Outdoor mixing fraction
REAL, DIMENSION(SIZE(PTI_BLD)):: ZT_MIX        ! Mixing air temperature [K]
REAL, DIMENSION(SIZE(PTI_BLD)):: ZQ_MIX        ! Mixing air specific humidity [kg kg-1]
!
REAL,DIMENSION(SIZE(PTI_BLD)) :: ZQI_BLD       ! Indoor air humidity at time step t + dt [K}
REAL, DIMENSION(SIZE(PTI_BLD)):: ZWASTE
!
REAL, DIMENSION(SIZE(PTI_BLD)):: ZDQS_FLOOR
REAL, DIMENSION(SIZE(PTI_BLD)):: ZIMB_FLOOR
REAL, DIMENSION(SIZE(PTI_BLD)):: ZDQS_MASS
REAL, DIMENSION(SIZE(PTI_BLD)):: ZIMB_MASS
REAL, DIMENSION(SIZE(PTI_BLD)):: ZLOAD_FLOOR   ! sum of solar and internal loads on floor
REAL, DIMENSION(SIZE(PTI_BLD)):: ZLOAD_MASS    ! sum of solar and internal loads on mass
REAL, DIMENSION(SIZE(PTI_BLD)):: ZRAD_FLOOR_MASS ! Rad. fluxes from floor to mass
REAL, DIMENSION(SIZE(PTI_BLD)):: ZCONV_FLOOR_BLD ! Conv. fluxes from floor to indoor air
REAL, DIMENSION(SIZE(PTI_BLD)):: ZCONV_MASS_BLD  ! Conv. fluxes from mass to indoor air
REAL, DIMENSION(SIZE(PTI_BLD)):: ZRHOI  ! indoor air density
!
INTEGER :: JJ                                  ! Loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!!REAL :: ZEXPL = 0.5 !explicit coefficient for internal temperature evol.
!!REAL :: ZIMPL = 0.5 !implicit coef..
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('BEM',0,ZHOOK_HANDLE)
!
!*      1.   Initializations
!            ---------------
!
ZRHOI(:) = PPS(:) / (XRD * PTI_BLD(:) * ( 1.+((XRV/XRD)-1.)*PQI_BLD(:) ) )
! *Temperal definitions for nat.vent*
!ZTCOMF_MAX  = 26. + 273.16
!
! *Definitions
ZFAN_AP (:) = 600.0
ZFAN_EFF(:) = 0.7
!
! *Other calcs
IROOF  = SIZE(PT_ROOF,2)
IWALL  = SIZE(PT_WALL_A,2)
!
!
! initial condition of QI_BLD equivalent to 50% RH
IF (ANY(PQI_BLD(:) <= 1E-6)) PQI_BLD = 0.5 * QSAT(PTI_BLD, PPS)
!
! *Temperal definitions for shedule*
GSCHED(:) = .FALSE.
WHERE (GSCHED(:))
 ZF_NIGHT(:)  = 0.8
 ZF_DAY(:)    = 1.2
ELSE WHERE
 ZF_NIGHT(:)  = 1.
 ZF_DAY(:)    = 1.
END WHERE
!
! *Int.gains schedule
!
ZQIN = PQIN * PN_FLOOR
WHERE (PSUNTIME(:) > 0. .AND. PSUNTIME(:) < 25200.) ! night between 0000 and 0700
  ZQIN(:) = ZQIN(:) * ZF_NIGHT(:)
ELSEWHERE
  ZQIN(:) = ZQIN(:) * ZF_DAY(:)
END WHERE

! *Change of units AC/H -> [m3 s-1 m-2(bld)]
ZV_VENT(:) = PV_VENT(:) * PBLD_HEIGHT(:) / 3600.
ZINF   (:) = PINF   (:) * PBLD_HEIGHT(:) / 3600.  
!
!*      2.   heat balance for building floor and mass
!            ----------------------------------------
!
!*      2.1 total load on the internal mass or floor
ZLOAD_FLOOR(:) = (ZQIN(:) * PQIN_FRAD(:) * (1.-PQIN_FLAT(:)) + PTR_SW_WIN(:)) / (PMASS_O_BLD(:)+1.)
WHERE (PN_FLOOR(:) > 1.)
   ZLOAD_MASS(:) = ZLOAD_FLOOR(:)
ELSEWHERE
   ZLOAD_MASS(:) = 0.
ENDWHERE
!
!*      2.2 FLOOR HEAT BALANCE
!
 CALL FLOOR_LAYER_E_BUDGET(PT_FLOOR, PTSTEP,               &
                          PHC_FLOOR, PTC_FLOOR, PD_FLOOR, &
                          PFLX_BLD_FLOOR, ZDQS_FLOOR, ZIMB_FLOOR,&
                          PF_FLOOR_MASS, PF_FLOOR_WALL, PF_FLOOR_WIN,&
                          PF_FLOOR_ROOF, PRADHT_IN, &
                          PT_MASS(:,1), PRAD_WALL_FLOOR, PRAD_ROOF_FLOOR, &
                          PRAD_WIN_FLOOR, PLOAD_IN_FLOOR, PTI_BLD, ZRAD_FLOOR_MASS, ZCONV_FLOOR_BLD)
!
!*      2.3 MASS HEAT BALANCE
!
 CALL MASS_LAYER_E_BUDGET(PT_MASS, PTSTEP,                &
                          PHC_FLOOR, PTC_FLOOR, PD_FLOOR/2., &
                          PFLX_BLD_MASS, ZDQS_MASS, ZIMB_MASS,&
                          PF_MASS_WALL, PF_MASS_WIN,&
                          PF_MASS_FLOOR, PRADHT_IN, &
                          PRAD_WALL_MASS, PRAD_ROOF_MASS,       &
                          PRAD_WIN_MASS, PLOAD_IN_MASS, PTI_BLD,             &
                          ZRAD_FLOOR_MASS, ZCONV_MASS_BLD  )
!
!
ZAC_IN_WALL_A_COOL = CHTC_VERT_DOE(PT_WALL_A(:,IWALL), PTCOOL_TARGET)
ZAC_IN_WALL_B_COOL = CHTC_VERT_DOE(PT_WALL_B(:,IWALL), PTCOOL_TARGET)
ZAC_IN_WIN_COOL = CHTC_VERT_DOE(PT_WIN2, PTCOOL_TARGET)
ZAC_IN_MASS_COOL = CHTC_VERT_DOE(PT_MASS(:,1), PTCOOL_TARGET)
ZAC_IN_ROOF_COOL = CHTC_DOWN_DOE(PT_ROOF(:,IROOF),PTCOOL_TARGET)
ZAC_IN_FLOOR_COOL = CHTC_UP_DOE(PT_FLOOR(:,1),PTCOOL_TARGET)

ZAC_IN_WALL_A_HEAT = CHTC_VERT_DOE(PT_WALL_A(:,IWALL), PTHEAT_TARGET)
ZAC_IN_WALL_B_HEAT = CHTC_VERT_DOE(PT_WALL_B(:,IWALL), PTHEAT_TARGET)
ZAC_IN_WIN_HEAT = CHTC_VERT_DOE(PT_WIN2, PTHEAT_TARGET)
ZAC_IN_MASS_HEAT = CHTC_VERT_DOE(PT_MASS(:,1), PTHEAT_TARGET)
ZAC_IN_ROOF_HEAT = CHTC_DOWN_DOE(PT_ROOF(:,IROOF),PTHEAT_TARGET)
ZAC_IN_FLOOR_HEAT = CHTC_UP_DOE(PT_FLOOR(:,1),PTHEAT_TARGET)

DO JJ=1,SIZE(ZAC_IN_WIN_COOL)
   ZAC_IN_WALL_A_COOL(JJ) = MAX(1.,ZAC_IN_WALL_A_COOL(JJ))
   ZAC_IN_WALL_B_COOL(JJ) = MAX(1., ZAC_IN_WALL_B_COOL(JJ))
   ZAC_IN_WIN_COOL(JJ)    = MAX(1., ZAC_IN_WIN_COOL(JJ))
   ZAC_IN_MASS_COOL(JJ)   = MAX(1., ZAC_IN_MASS_COOL(JJ))
   ZAC_IN_ROOF_COOL(JJ)   = MAX(1., ZAC_IN_ROOF_COOL(JJ))
   ZAC_IN_FLOOR_COOL(JJ)  = MAX(1., ZAC_IN_FLOOR_COOL(JJ))
   
   ZAC_IN_WALL_A_HEAT(JJ) = MAX(1.,ZAC_IN_WALL_A_HEAT(JJ))
   ZAC_IN_WALL_B_HEAT(JJ) = MAX(1., ZAC_IN_WALL_B_HEAT(JJ))
   ZAC_IN_WIN_HEAT(JJ)    = MAX(1., ZAC_IN_WIN_HEAT(JJ))
   ZAC_IN_MASS_HEAT(JJ)   = MAX(1., ZAC_IN_MASS_HEAT(JJ))
   ZAC_IN_ROOF_HEAT(JJ)   = MAX(1., ZAC_IN_ROOF_HEAT(JJ))
   ZAC_IN_FLOOR_HEAT(JJ)  = MAX(1., ZAC_IN_FLOOR_HEAT(JJ))
ENDDO

!*      4.   Indoor energy balance calculation
!            ---------------------------------
!
DO JJ=1,SIZE(PT_CANYON)
  ! *first guess of indoor temperature

  ZTI_BLD(JJ) = PTI_BLD(JJ)                                          &
          + PTSTEP/(ZRHOI(JJ) * XCPD * PBLD_HEIGHT(JJ))              & 
          * (  PWALL_O_BLD(JJ) * PCONV_WALL_BLD(JJ)                        &
             + PGLAZ_O_BLD (JJ) * PCONV_WIN_BLD(JJ)                        &
             + PMASS_O_BLD(JJ) * ZCONV_MASS_BLD(JJ)                        &
             + PCONV_ROOF_BLD(JJ)                                    &
             + ZCONV_FLOOR_BLD(JJ)                                   &
             + ZQIN(JJ) * (1 - PQIN_FRAD(JJ))  * (1 - PQIN_FLAT(JJ)) )
  !
  !################################################################################
  ! *is natural surventilation active at the current time step ?
  !---------------------------------------------------------------------------------
  !
  !    *no surventilation possible

  IF (HNATVENT(JJ)=='NONE') THEN
    !
    GNAT_VENT(JJ) = .FALSE.
    !
  !    *automatic management of surventilation
  ELSEIF (HNATVENT(JJ)=='AUTO' .OR. HNATVENT(JJ)=='MECH') THEN
    !
    IF (MOD(PSUNTIME(JJ), 3600.) .LT. PTSTEP) THEN
      !
      IF ( PTI_BLD(JJ).GT. PT_CANYON(JJ) + 1 ) THEN ! condition to enable the
        IF (HNATVENT(JJ)=='AUTO') THEN
        ! natural surventilation rate calculation (window opening)
          CALL GET_NAT_VENT(PTI_BLD(JJ), PT_CANYON(JJ), PU_CANYON(JJ), PGR(JJ), &
                            PFLOOR_HW_RATIO(JJ), PBLD_HEIGHT(JJ), ZNAT_VENT(JJ))
        ELSE IF (HNATVENT(JJ)=='MECH') THEN
        ! mechanical surventilation rate calculation : 5 volumes/hour
          ZNAT_VENT(JJ) =  5.0*PBLD_HEIGHT(JJ)/3600.
        END IF
        !
        ZTI_BLD_OPEN  (JJ) = ZTI_BLD(JJ) &
                + ZNAT_VENT(JJ)            * PTSTEP/PBLD_HEIGHT(JJ) * (PT_CANYON(JJ) - PTI_BLD(JJ)) 
        ZTI_BLD_CLOSED(JJ) = ZTI_BLD(JJ) &
                + (ZINF(JJ) + ZV_VENT(JJ)) * PTSTEP/PBLD_HEIGHT(JJ) * (PT_CANYON(JJ) - PTI_BLD(JJ)) 
        !
        GNAT_VENT(JJ) = (ZTI_BLD_OPEN(JJ) <= PTCOOL_TARGET (JJ) .AND. &            
                         ZTI_BLD_OPEN(JJ) <  ZTI_BLD_CLOSED(JJ) .AND. &
                         ZTI_BLD_OPEN(JJ) >  PTHEAT_TARGET (JJ) + 4.)
        !
      ELSE
        GNAT_VENT(JJ) = .FALSE.
      ENDIF
      ONATVENT_NIGHT(JJ) = GNAT_VENT(JJ)
    ELSE 
      GNAT_VENT(JJ) = ONATVENT_NIGHT(JJ)
    ENDIF
    !
  !    *manual management of surventilation
  ELSEIF (HNATVENT(JJ)=='MANU') THEN
    !
    ONATVENT_NIGHT(JJ) = ONATVENT_NIGHT(JJ) .AND. &
                         .NOT. ( PSUNTIME(JJ) > 5.*3600 .AND. PSUNTIME(JJ) < 18.*3600 )
    !
    GNAT_VENT(JJ) = ( PSUNTIME(JJ) > 18.*3600. .AND. PSUNTIME(JJ) < 21.*3600.  &
                      .AND. PT_CANYON(JJ) < PTI_BLD(JJ)+2.       &
                      .AND. PT_CANYON(JJ) > PTHEAT_TARGET(JJ)    & 
                      .AND. ( PTI_BLD(JJ) > PTHEAT_TARGET(JJ)+5. &
                       .OR. PTI_BLD(JJ) == PTCOOL_TARGET(JJ) )   ) 
    GNAT_VENT(JJ) = GNAT_VENT(JJ) .OR. ONATVENT_NIGHT(JJ)
    !
  ENDIF
  !
  ! Decicion about natural surventilation OK
  !################################################################################
  !
  !
  !################################################################################
  ! COMPUTE ENERGY DEMAND
  !---------------------------------------------------------------------------------

  ! *If natural surventilation ACTIVE
  IF (GNAT_VENT(JJ)) THEN
     ! 
     CALL GET_NAT_VENT(PTI_BLD(JJ), PT_CANYON(JJ), PU_CANYON(JJ), PGR(JJ), &
                       PFLOOR_HW_RATIO(JJ), PBLD_HEIGHT(JJ), ZNAT_VENT(JJ)         )
     !
     ZV_VENT     (JJ) = 0.
     ZINF        (JJ) = 0.
     !
     PH_BLD_COOL (JJ) = 0.0         ! No HVAC consumption
     PH_BLD_HEAT (JJ) = 0.0    
     PLE_BLD_COOL(JJ) = 0.0         ! No HVAC consumption
     PLE_BLD_HEAT(JJ) = 0.0   
     !    
     PT_BLD_COOL (JJ) = 0.0         ! No HVAC consumption
     PHVAC_COOL  (JJ) = 0.0    
     PT_SYS      (JJ) = PTI_BLD(JJ) ! No mechanical ventilation
     PQ_SYS      (JJ) = PQI_BLD(JJ) ! 
     PH_WASTE    (JJ) = 0.0
     PLE_WASTE   (JJ) = 0.0     
     PFAN_POWER  (JJ) = 0.0    
     PHVAC_HEAT  (JJ) = 0.0
     !
     PM_SYS  (JJ) = 0.0
     PCOP    (JJ) = 0.0
     PCAP_SYS(JJ) = 0.0
     !
  ! *If natural surventilation INACTIVE
  ELSE 
     !
     ZNAT_VENT(JJ) = 0.
     !
     ! ------------------------------------------------
     ! * Building energy demand for heating and cooling
     ! ------------------------------------------------
     !
     PH_BLD_COOL(JJ) = PWALL_O_BLD(JJ)/2. * (ZAC_IN_WALL_A_COOL(JJ) * (PT_WALL_A(JJ,IWALL) - PTCOOL_TARGET(JJ))  &
                                           + ZAC_IN_WALL_B_COOL(JJ) * (PT_WALL_B(JJ,IWALL) - PTCOOL_TARGET(JJ))) &
                     + PGLAZ_O_BLD(JJ) * ZAC_IN_WIN_COOL(JJ) * (PT_WIN2(JJ)       - PTCOOL_TARGET(JJ)) &    
                + ZAC_IN_MASS_COOL (JJ)* PMASS_O_BLD(JJ) * (PT_MASS(JJ,1)     - PTCOOL_TARGET(JJ))  &
                + ZAC_IN_ROOF_COOL (JJ)              * (PT_ROOF(JJ,IROOF) - PTCOOL_TARGET(JJ))  &
                + ZAC_IN_FLOOR_COOL(JJ)              * (PT_FLOOR(JJ,1)    - PTCOOL_TARGET(JJ))  &
                + ZQIN(JJ) * (1 - PQIN_FRAD(JJ)) * (1 - PQIN_FLAT(JJ))                     &
                + (ZINF(JJ) + ZV_VENT(JJ)) * ZRHOI(JJ) * XCPD * (PT_CANYON(JJ) - PTCOOL_TARGET(JJ))
     !
     PH_BLD_HEAT(JJ) = - ( PWALL_O_BLD(JJ)/2. * (ZAC_IN_WALL_A_HEAT(JJ) * (PT_WALL_A(JJ,IWALL) - PTHEAT_TARGET(JJ))  &
                                                +ZAC_IN_WALL_B_HEAT(JJ) * (PT_WALL_B(JJ,IWALL) - PTHEAT_TARGET(JJ)))  &
                         + PGLAZ_O_BLD(JJ) * ZAC_IN_WIN_HEAT(JJ) * (PT_WIN2(JJ)       - PTHEAT_TARGET(JJ)) &    
                +  ZAC_IN_MASS_HEAT(JJ)* PMASS_O_BLD(JJ) * (PT_MASS(JJ,1)     - PTHEAT_TARGET(JJ))  &
                +  ZAC_IN_ROOF_HEAT(JJ)              * (PT_ROOF(JJ,IROOF) - PTHEAT_TARGET(JJ))  &
                + ZAC_IN_FLOOR_HEAT(JJ)              * (PT_FLOOR(JJ,1)    - PTHEAT_TARGET(JJ))  &
                + ZQIN(JJ) * (1 - PQIN_FRAD(JJ))* (1 - PQIN_FLAT(JJ))                      &
                + (ZINF(JJ) + ZV_VENT(JJ)) * ZRHOI(JJ) * XCPD * (PT_CANYON(JJ) - PTHEAT_TARGET(JJ)))
     !

     ZQCOOL_TARGET(JJ) = 0.62198 * PHR_TARGET(JJ) * PSAT(PTCOOL_TARGET(JJ)) / &
                         (PPS(JJ)- PHR_TARGET(JJ) * PSAT(PTCOOL_TARGET(JJ)))    
     !
     PLE_BLD_COOL(JJ) = ZQIN(JJ) * PQIN_FLAT(JJ)                                           &
                + (ZINF(JJ) + ZV_VENT(JJ)) * ZRHOI(JJ) * XLVTT * (PQ_CANYON(JJ) - ZQCOOL_TARGET(JJ)) 
     !

     ZQHEAT_TARGET(JJ) = 0.62198 * PHR_TARGET(JJ) * PSAT(PTHEAT_TARGET(JJ)) / &
                         (PPS(JJ)- PHR_TARGET(JJ) * PSAT(PTHEAT_TARGET(JJ)))    
     !
     PLE_BLD_HEAT(JJ) = ZQIN(JJ) * PQIN_FLAT(JJ)                                           &
                + (ZINF(JJ) + ZV_VENT(JJ)) * ZRHOI(JJ) * XLVTT * (PQ_CANYON(JJ) - ZQHEAT_TARGET(JJ))       
     !
     ! * Autosize calculations
     !
     IF (OAUTOSIZE .AND. KDAY==15) THEN
        !
        IF (PH_BLD_COOL(JJ) > PAUX_MAX(JJ))  THEN
          !  
          PAUX_MAX    (JJ) = PH_BLD_COOL(JJ)
          !
          ! Cooling coil sensible heat rate 
          ZSHR        (JJ) = MIN(XCPD * (PTCOOL_TARGET(JJ) - PT_ADP(JJ)) /           &
                                (ENTH_FN_T_Q(PTCOOL_TARGET(JJ),ZQCOOL_TARGET(JJ)) -  &
                                 ENTH_FN_T_Q(PT_ADP(JJ),QSAT(PT_ADP(JJ),PPS(JJ)))), 1.)
          ! Cooling Coil Capacity [W m-2(bld)]
          PCAP_SYS_RAT(JJ) = PH_BLD_COOL(JJ) / ZSHR(JJ) 
          !
          ! Cooling rated air flow rate [kg s-1 m-2(bld)]
          ZM_SYS_RAT  (JJ) = PH_BLD_COOL(JJ) / XCPD / (PTCOOL_TARGET(JJ)-(14.0+273.16))
          IF (ZM_SYS_RAT(JJ) > PM_SYS_RAT(JJ)) PM_SYS_RAT(JJ) = ZM_SYS_RAT(JJ)
          !
          ! Impose condition 
          IF (PM_SYS_RAT(JJ)/ZRHOI(JJ)/PCAP_SYS_RAT(JJ) < 0.00004027) THEN
            PCAP_SYS_RAT(JJ) = PM_SYS_RAT(JJ)/ZRHOI(JJ)/0.00004027
          ELSE IF (PM_SYS_RAT(JJ)/ZRHOI(JJ)/PCAP_SYS_RAT(JJ) > 0.00006041) THEN
            PCAP_SYS_RAT(JJ) = PM_SYS_RAT(JJ)/ZRHOI(JJ)/0.00006041
          END IF
          !
        END IF
        !
     END IF
     !
     ! * END Autosize calculations
     !
     ! * system efficiency
     ! ...................
     !
     PM_SYS  (JJ) = PM_SYS_RAT  (JJ)
     PCOP    (JJ) = PCOP_RAT    (JJ)
     PCAP_SYS(JJ) = PCAP_SYS_RAT(JJ)
     !
     ! * Mixing conditions
     ! .................
     !
     ZXMIX (JJ) = ZV_VENT(JJ) * ZRHOI(JJ) / PM_SYS(JJ)
     ZT_MIX(JJ) = ZXMIX(JJ) * PT_CANYON(JJ) + (1.-ZXMIX(JJ)) * PTI_BLD(JJ)
     ZQ_MIX(JJ) = ZXMIX(JJ) * PQ_CANYON(JJ) + (1.-ZXMIX(JJ)) * PQI_BLD(JJ)
     ! 
     ! ---------------------------------------------
     ! * COOLING system : Performance and Waste heat
     ! ---------------------------------------------
     !
     IF (PH_BLD_COOL(JJ) >= 0.0) THEN
        !
        ! *ideal system
        IF (HCOOL_COIL=='IDEAL') THEN
           !
           PT_BLD_COOL(JJ) = PH_BLD_COOL(JJ) + PLE_BLD_COOL(JJ)
           !desactivation of LE_BLD_COOL impact on HVAC_COOL calculation
           !following too much impact in VURCA simulation (23/01/2012)
           !this would be the case for a vaporization system !
           !PHVAC_COOL (JJ) = PT_BLD_COOL(JJ) / PCOP_RAT(JJ)
           PHVAC_COOL (JJ) = PH_BLD_COOL(JJ) / PCOP_RAT(JJ)
           IF (PHVAC_COOL(JJ) < 0.0) PHVAC_COOL(JJ) = 0.0
           !
           PT_SYS(JJ) = ZT_MIX(JJ) - PH_BLD_COOL (JJ)  /PM_SYS(JJ) / XCPD
           !PQ_SYS(JJ) = ZQ_MIX(JJ) - PLE_BLD_COOL(JJ) / PM_SYS(JJ)/ XLVTT
           !desactivation following too much impact in VURCA simulation
           !(23/01/2012)
           PQ_SYS(JJ) = ZQ_MIX(JJ)
           !
           PH_WASTE(JJ)  = PHVAC_COOL(JJ) * (1.+PCOP_RAT(JJ)) * (1. - PF_WATER_COND(JJ))
           PLE_WASTE(JJ) = PHVAC_COOL(JJ) * (1.+PCOP_RAT(JJ)) * PF_WATER_COND(JJ)
           !
        ! *real system
        ELSEIF (HCOOL_COIL=='DXCOIL') THEN
           !
           CALL DX_AIR_COOLING_COIL_CV(PT_CANYON(JJ), PQ_CANYON(JJ), PPS(JJ),  &
                     ZRHOI(JJ), ZT_MIX(JJ), ZQ_MIX(JJ), PCOP_RAT(JJ),          &
                     PCAP_SYS_RAT(JJ), PT_ADP(JJ), PF_WATER_COND(JJ),          &
                     PM_SYS(JJ), PH_BLD_COOL(JJ), PH_WASTE(JJ), PLE_WASTE(JJ), &
                     PCOP(JJ), PCAP_SYS(JJ), PT_SYS(JJ), PQ_SYS(JJ),           &
                     PHVAC_COOL(JJ), PT_BLD_COOL(JJ)                           )
           !
        ENDIF !end type of cooling system

        !!! case of system without atmospheric releases. I-e releases in soil/water F_WATER_COND < 0 
        IF (PF_WATER_COND(JJ) < 0) THEN
          PH_WASTE(JJ) = 0. 
          PLE_WASTE(JJ) = 0. 
        ENDIF
        !!!!
        !
        !         From EP Engineering Reference (p. 647)
        PFAN_POWER(JJ) = PM_SYS(JJ) * ZFAN_AP(JJ) * ZFAN_EFF(JJ) * ZRHOI(JJ)
        !
        PH_BLD_HEAT (JJ) = 0.0
        PLE_BLD_HEAT(JJ) = 0.0
        PHVAC_HEAT  (JJ) = 0.0
     !
     ! ---------------------------------------------
     ! * HEATING system : Performance and Waste heat
     ! ---------------------------------------------
     !
     ELSE IF (PH_BLD_HEAT(JJ) > 0.0) THEN
        !
        ! *specific computation for real heating system
        IF  (HHEAT_COIL .EQ. 'FINCAP') THEN
          IF (PH_BLD_HEAT(JJ) > PCAP_SYS_HEAT(JJ)) PH_BLD_HEAT(JJ) =  PCAP_SYS_HEAT(JJ)
        END IF
        !
        PT_SYS(JJ) = ZT_MIX(JJ) + PH_BLD_HEAT(JJ) / PM_SYS(JJ) / XCPD
        PQ_SYS(JJ) = ZQ_MIX(JJ)
        !
        PHVAC_HEAT  (JJ) = PH_BLD_HEAT(JJ) / PEFF_HEAT(JJ)
        PH_WASTE    (JJ) = PHVAC_HEAT(JJ) - PH_BLD_HEAT(JJ)  
        PLE_WASTE   (JJ) = 0.0
        PH_BLD_COOL (JJ) = 0.0
        PLE_BLD_COOL(JJ) = 0.0
        PT_BLD_COOL (JJ) = 0.0
        PHVAC_COOL  (JJ) = 0.0
!       From EP Engineering Reference (p. 647)
        PFAN_POWER(JJ) = PM_SYS(JJ)*ZFAN_AP(JJ)*(ZFAN_EFF(JJ)*ZRHOI(JJ))
     !
     ! ------------------------------
     ! * NEITHEIR COOLING NOR HEATING
     ! ------------------------------
     ! 
     ELSE
        !
        PH_BLD_COOL (JJ) = 0.0
        PH_BLD_HEAT (JJ) = 0.0 
        PLE_BLD_COOL(JJ) = 0.0
        PLE_BLD_HEAT(JJ) = 0.0
        !
        PT_BLD_COOL (JJ) = 0.0 
        PHVAC_COOL  (JJ) = 0.0
        PT_SYS      (JJ) = ZT_MIX(JJ)
        PQ_SYS      (JJ) = ZQ_MIX(JJ)
        PH_WASTE    (JJ) = 0.0
        PLE_WASTE   (JJ) = 0.0
        PFAN_POWER  (JJ) = 0.0     
        PHVAC_HEAT  (JJ) = 0.0
        !
     END IF !end for heating/cooling sytem
     !
  END IF
  !
  !---------------------------------------------------------------------------------
  ! ENERGY DEMAND COMPUTED
  !################################################################################
ENDDO
!
!---------------------------------------------------
! EVOLUTION OF THE INTERNAL TEMPERATURE AND HUMIDITY
!###################################################
!
ZTI_BLD(:) = ( ZTI_BLD(:) + PTSTEP/PBLD_HEIGHT(:) *                   & 
            ((ZINF(:) + ZNAT_VENT(:)) * (PT_CANYON(:)             )   &
            + PM_SYS(:) / ZRHOI(:)    * (PT_SYS   (:)             ) ))&
          / (1. + PTSTEP/PBLD_HEIGHT(:)*                              &
                (ZINF(:) + ZNAT_VENT(:) + PM_SYS(:) / ZRHOI(:))      )
ZQI_BLD(:) = ( PQI_BLD(:) +  PTSTEP/PBLD_HEIGHT(:) *                    & 
             ( ZQIN(:) * PQIN_FLAT(:) / ( ZRHOI(:) * XLVTT)             &
              + (ZINF(:) + ZNAT_VENT(:)) * (PQ_CANYON(:)            )   &
              + PM_SYS(:) / ZRHOI(:)     * (PQ_SYS   (:)            ) ))&
          / (1. + PTSTEP/PBLD_HEIGHT(:)*                              &
                (ZINF(:) + ZNAT_VENT(:) + PM_SYS(:) / ZRHOI(:))      )
!
!
! Update variables
PTI_BLD(:) = ZTI_BLD(:)
PQI_BLD(:) = ZQI_BLD(:)


!
! Waste heat due to infiltration/ventilation
ZWASTE   (:) = (ZINF(:)+ZV_VENT(:)+ZNAT_VENT(:)) * ZRHOI(:) 
PH_WASTE (:) = PH_WASTE (:) + ZWASTE(:) * XCPD  * (PTI_BLD(:) - PT_CANYON(:))
PLE_WASTE(:) = PLE_WASTE(:) + ZWASTE(:) * XLVTT * (PQI_BLD(:) - PQ_CANYON(:)) 
!
!
IF (LHOOK) CALL DR_HOOK('BEM',1,ZHOOK_HANDLE)
!
 CONTAINS
!
SUBROUTINE GET_NAT_VENT(PPTI_BLD, PPT_CANYON, PPU_CANYON, PPGR, &
                        PF_AUX, PPBLD_HEIGHT, PNAT_VENT)
!
IMPLICIT NONE
!
REAL, INTENT(IN) :: PPTI_BLD
REAL, INTENT(IN) :: PPT_CANYON
REAL, INTENT(IN) :: PPU_CANYON
REAL, INTENT(IN) :: PPGR
REAL, INTENT(IN) :: PF_AUX
REAL, INTENT(IN) :: PPBLD_HEIGHT
REAL, INTENT(OUT) :: PNAT_VENT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('BEM:GET_NAT_VENT',0,ZHOOK_HANDLE)
!
PNAT_VENT = XG * (PPTI_BLD - PPT_CANYON)
IF (PNAT_VENT .LT. 0.) THEN ! exceptional case with MANU ventilation system
   PNAT_VENT= PPBLD_HEIGHT/3600. !minimum value
ELSE
   PNAT_VENT = 1./3. * (PNAT_VENT/PPT_CANYON)**(1./2.)                             &
               * (1.5 + PPTI_BLD/PNAT_VENT * 1./2. * PPU_CANYON**2*0.1)**(3./2.) &
               * PPGR * PF_AUX / 1.5 / 2.
   PNAT_VENT = MIN(PNAT_VENT, 5.0*PPBLD_HEIGHT/3600.)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('BEM:GET_NAT_VENT',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_NAT_VENT
!
END SUBROUTINE BEM
