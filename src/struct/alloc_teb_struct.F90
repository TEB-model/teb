!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL licence
!SFX_LIC version 2.1. See Licence_CeCILL_V2.1-en.txt and Licence_CeCILL_V2.1-fr.txt  
!SFX_LIC for details.
!   ##########################################################################
    SUBROUTINE ALLOC_TEB_STRUCT(KROOF_LAYER,KROAD_LAYER,KWALL_LAYER,KFLOOR_LAYER,&
                               DTCO, G, TOP, T, BOP, B, TPN, TIR, DMT, GDM, GRM )
!   ##########################################################################
!
!!****  *ALLOC_TEB_STRUCT*
!!
!!    PURPOSE
!!    -------
!
!       To allocate all the structures used in Teb_garden and subsequent routines.
!
!!**  METHOD
!     ------
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
!!      V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!    Original    03/2018
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t, DATA_COVER_INIT
USE MODD_SFX_GRID_n, ONLY : GRID_t, GRID_INIT
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t, TEB_OPTIONS_INIT
USE MODD_TEB_n, ONLY : TEB_t, TEB_INIT
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t, BEM_OPTIONS_INIT
USE MODD_BEM_n, ONLY : BEM_t, BEM_INIT
USE MODD_TEB_PANEL_n, ONLY : TEB_PANEL_t, TEB_PANEL_INIT
USE MODD_TEB_IRRIG_n, ONLY : TEB_IRRIG_t, TEB_IRRIG_INIT
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_t, DIAG_MISC_TEB_INIT
USE MODD_SURFEX_n, ONLY : TEB_GARDEN_MODEL_t
USE MODD_SURFEX_n, ONLY : TEB_GREENROOF_MODEL_t
!
!
IMPLICIT NONE
!
!*      0.1    Declarations of arguments
!
!
!
!
INTEGER, INTENT(IN) :: KROOF_LAYER   ! number of roof layers
INTEGER, INTENT(IN) :: KROAD_LAYER   ! number of road layers
INTEGER, INTENT(IN) :: KWALL_LAYER   ! number of wall layers
INTEGER, INTENT(IN) :: KFLOOR_LAYER  ! number of floor layers
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(GRID_t), INTENT(INOUT) :: G
TYPE(TEB_OPTIONS_t),   INTENT(INOUT) :: TOP
TYPE(TEB_t),           INTENT(INOUT) :: T
TYPE(BEM_OPTIONS_t),   INTENT(INOUT) :: BOP
TYPE(BEM_t),           INTENT(INOUT) :: B
TYPE(TEB_PANEL_t),     INTENT(INOUT) :: TPN
TYPE(TEB_IRRIG_t),     INTENT(INOUT) :: TIR
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: DMT
TYPE(TEB_GARDEN_MODEL_t), INTENT(INOUT) :: GDM
TYPE(TEB_GREENROOF_MODEL_t), INTENT(INOUT) :: GRM
!
!
!*      0.2    Declarations of local variables
!
!
!-------------------------------------------------------------------------------
!
CALL GRID_INIT(G)
CALL TEB_INIT(T)
CALL TEB_PANEL_INIT(TPN)
CALL TEB_OPTIONS_INIT(TOP)
CALL TEB_IRRIG_INIT(TIR)
CALL DIAG_MISC_TEB_INIT(DMT)
CALL BEM_INIT(B)
CALL BEM_OPTIONS_INIT(BOP)
CALL DATA_COVER_INIT(DTCO)
!
!-------------------------------------------------------------------------------
!
! BEM structure
!
!
! Floor parameters
!
ALLOCATE(B%XHC_FLOOR(1,KFLOOR_LAYER)) ! floor layers heat capacity        (J/K/m3)
ALLOCATE(B%XTC_FLOOR(1,KFLOOR_LAYER)) ! floor layers thermal conductivity (W/K/m)
ALLOCATE(B%XD_FLOOR (1,KFLOOR_LAYER)) ! depth of floor layers             (m)
!
! HVAC parameters
!
ALLOCATE(B%XTCOOL_TARGET(1)) ! cooling setpoint of indoor air
ALLOCATE(B%XTHEAT_TARGET(1)) ! heating setpoint of indoor air
ALLOCATE(B%XF_WASTE_CAN (1)) ! fraction of waste heat released into the canyon
ALLOCATE(B%XEFF_HEAT    (1)) ! efficiency of the heating system
!
! Indoor parameters
!
ALLOCATE(B%XTI_BLD(1))                    ! building interior temperature    (K)
ALLOCATE(B%XT_FLOOR (1,KFLOOR_LAYER)) ! floor layer temperatures         (K)
ALLOCATE(B%XT_MASS  (1,KFLOOR_LAYER)) ! Air cooled building internal th. mass temperature (K)
!
ALLOCATE(B%XQIN(1))           ! internal heat gains [W m-2(floor)]
ALLOCATE(B%XQIN_FRAD(1))      ! radiant fraction of internal heat gains
ALLOCATE(B%XSHGC(1))          ! solar heat gain coef. of windows
ALLOCATE(B%XSHGC_SH(1))       ! solar heat gain coef. of windows + shading
ALLOCATE(B%XU_WIN(1))         ! window U-factor [K m W-2]
ALLOCATE(B%XTRAN_WIN(1))      ! window transmittance (-)
ALLOCATE(B%XGR(1))            ! glazing ratio
ALLOCATE(B%XFLOOR_HEIGHT(1))  ! building floor height [m]
ALLOCATE(B%XINF(1))           ! infiltration/ventilation flow rate [AC/H]
!
! New parameters
!
ALLOCATE(B%XF_WATER_COND(1)) ! fraction of evaporation for condensers (cooling system)
ALLOCATE(B%XAUX_MAX(1))      ! Auxiliar variable for autosize calcs
ALLOCATE(B%XQIN_FLAT(1))     ! Latent franction of internal heat gains
ALLOCATE(B%XHR_TARGET(1))    ! Relative humidity setpoint
ALLOCATE(B%XT_WIN2(1))       ! Indoor window temperature [K]
ALLOCATE(B%XQI_BLD(1))       ! Indoor air specific humidity [kg kg-1]
ALLOCATE(B%XV_VENT(1))       ! Ventilation flow rate [AC/H]
ALLOCATE(B%XCAP_SYS_HEAT(1)) ! Capacity of the heating system 
                             ! [W m-2(bld)]
ALLOCATE(B%XCAP_SYS_RAT(1))  ! Rated capacity of the cooling system
                             ! [W m-2(bld)]
ALLOCATE(B%XT_ADP(1))        ! Apparatus dewpoint temperature of the
                             ! cooling coil [K]
ALLOCATE(B%XM_SYS_RAT(1))    ! Rated HVAC mass flow rate 
                             ! [kg s-1 m-2(bld)]
ALLOCATE(B%XCOP_RAT(1))      ! Rated COP of the cooling system
ALLOCATE(B%XT_WIN1(1))       ! outdoor window temperature [K]
ALLOCATE(B%XALB_WIN(1))      ! window albedo
ALLOCATE(B%XABS_WIN(1))      ! window absortance
ALLOCATE(B%XT_SIZE_MAX(1))   ! Maximum outdoor air temperature for
                             ! HVAC sizing [K]
ALLOCATE(B%XT_SIZE_MIN(1))   ! Minimum outdoor air temperature for
                             ! HVAC sizing [K]
ALLOCATE(B%XUGG_WIN(1))      ! Window glass-to-glass U-factor [K m W-2]
ALLOCATE(B%LSHADE(1))        ! flag to activate shading devices -> LOGICAL in the code
ALLOCATE(B%XSHADE(1))        ! flag to activate shading devices -> REAL for i/o 0. or 1.
ALLOCATE(B%CNATVENT(1))      ! flag to activate natural ventilation 'NONE', 'MANU', 'AUTO'
ALLOCATE(B%XNATVENT(1))      ! flag to describe surventilation system for i/o 
                             ! 0 for NONE, 1 for MANU and 2 for AUTO
ALLOCATE(B%LSHAD_DAY(1))     !Has shading been necessary this day ?
ALLOCATE(B%LNATVENT_NIGHT(1))!Has nocturnal surventilation been necessary and possible this night ?
  !
  !indoor relative surfaces and view factors
ALLOCATE(B%XN_FLOOR(1))        ! Number of floors     
ALLOCATE(B%XGLAZ_O_BLD(1))     ! Window area [m2_win/m2_bld]
ALLOCATE(B%XMASS_O_BLD(1))     ! Mass area [m2_mass/m2_bld]
ALLOCATE(B%XFLOOR_HW_RATIO(1)) ! H/W ratio of 1 floor level
ALLOCATE(B%XF_FLOOR_MASS(1))   ! View factor floor-mass
ALLOCATE(B%XF_FLOOR_WALL(1))   ! View factor floor-wall
ALLOCATE(B%XF_FLOOR_WIN(1))    ! View factor floor-window
ALLOCATE(B%XF_FLOOR_ROOF(1))   ! View factor floor-roof
ALLOCATE(B%XF_WALL_FLOOR(1))   ! View factor wall-floor
ALLOCATE(B%XF_WALL_MASS(1))    ! View factor wall-mass
ALLOCATE(B%XF_WALL_WIN(1))     ! View factor wall-win
ALLOCATE(B%XF_WIN_FLOOR(1))    ! View factor win-floor
ALLOCATE(B%XF_WIN_MASS(1))     ! View factor win-mass
ALLOCATE(B%XF_WIN_WALL(1))     ! View factor win-wall
ALLOCATE(B%XF_WIN_WIN(1))      ! indoor View factor win-win
ALLOCATE(B%XF_MASS_FLOOR(1))   ! View factor mass-floor
ALLOCATE(B%XF_MASS_WALL(1))    ! View factor mass-wall
ALLOCATE(B%XF_MASS_WIN(1))     ! View factor mass-window
!
!-------------------------------------------------------------------------------
!
! TEB structure
!
! Geometric Parameters:
!
ALLOCATE(T%XROAD_DIR(1))     ! Road direction (deg from North, clockwise)
ALLOCATE(T%XGARDEN(1))       ! fraction of veg in the streets   (-)
ALLOCATE(T%XGREENROOF(1))    ! fraction of greenroofs on roofs  (-)
ALLOCATE(T%XBLD(1))          ! fraction of buildings            (-)
ALLOCATE(T%XROAD(1))         ! fraction of roads                (-)
ALLOCATE(T%XCAN_HW_RATIO(1)) ! canyon    h/W                    (-)
ALLOCATE(T%XBLD_HEIGHT(1))   ! buildings height 'h'             (m)
ALLOCATE(T%XWALL_O_HOR(1))   ! wall surf. / hor. surf.          (-)
ALLOCATE(T%XROAD_O_GRND(1))  ! road surf. / (road + garden surf.) (-)
ALLOCATE(T%XGARDEN_O_GRND(1))! gard. surf. / (road + garden surf.)(-)
ALLOCATE(T%XWALL_O_GRND(1))  ! wall surf. / (road + garden surf.) (-)
ALLOCATE(T%XWALL_O_BLD(1))   ! wall surf. / bld surf. (-)
ALLOCATE(T%XZ0_TOWN(1))      ! roughness length for momentum    (m)
ALLOCATE(T%XSVF_ROAD(1))     ! road sky view factor             (-)
ALLOCATE(T%XSVF_GARDEN(1))   ! green area sky view factor       (-)
ALLOCATE(T%XSVF_WALL(1))     ! wall sky view factor             (-)
!
! Roof parameters
!
ALLOCATE(T%XALB_ROOF(1))     ! roof albedo                      (-)
ALLOCATE(T%XEMIS_ROOF(1))    ! roof emissivity                  (-)
ALLOCATE(T%XHC_ROOF(1,KROOF_LAYER))      ! roof layers heat capacity        (J/K/m3)
ALLOCATE(T%XTC_ROOF(1,KROOF_LAYER))      ! roof layers thermal conductivity (W/K/m)
ALLOCATE(T%XD_ROOF(1,KROOF_LAYER))       ! depth of roof layers             (m)
ALLOCATE(T%XROUGH_ROOF(1))   ! roof roughness coef
!
!
! Road parameters
!
ALLOCATE(T%XALB_ROAD(1))     ! road albedo                      (-)
ALLOCATE(T%XEMIS_ROAD(1))    ! road emissivity                  (-)
ALLOCATE(T%XHC_ROAD(1,KROAD_LAYER))      ! road layers heat capacity        (J/K/m3)
ALLOCATE(T%XTC_ROAD(1,KROAD_LAYER))      ! road layers thermal conductivity (W/K/m)
ALLOCATE(T%XD_ROAD(1,KROAD_LAYER))       ! depth of road layers             (m)
!
! Wall parameters
!
ALLOCATE(T%XALB_WALL(1))     ! wall albedo                      (-)
ALLOCATE(T%XEMIS_WALL(1))    ! wall emissivity                  (-)
ALLOCATE(T%XHC_WALL(1,KWALL_LAYER))      ! wall layers heat capacity        (J/K/m3)
ALLOCATE(T%XTC_WALL(1,KWALL_LAYER))      ! wall layers thermal conductivity (W/K/m)
ALLOCATE(T%XD_WALL(1,KWALL_LAYER))       ! depth of wall layers             (m)
ALLOCATE(T%XROUGH_WALL(1))   ! wall roughness coef
!
! Building's use type
!
ALLOCATE(T%XRESIDENTIAL(1))  ! fraction of Residential use      (-)
  
!
! anthropogenic fluxes
!
ALLOCATE(T%XH_TRAFFIC(1))    ! anthropogenic sensible
!                                                  ! heat fluxes due to traffic       (W/m2)
ALLOCATE(T%XLE_TRAFFIC(1))   ! anthropogenic latent
!                                                  ! heat fluxes due to traffic       (W/m2)
ALLOCATE(T%XH_INDUSTRY(1))   ! anthropogenic sensible                   
!                                                  ! heat fluxes due to factories     (W/m2)
ALLOCATE(T%XLE_INDUSTRY(1))  ! anthropogenic latent
!                                                  ! heat fluxes due to factories     (W/m2)
!
! temperatures for boundary conditions
!
ALLOCATE(T%XTI_ROAD(1))      ! road interior temperature        (K)
!
! Prognostic variables:
!
ALLOCATE(T%XWS_ROOF(1))      ! roof water reservoir             (kg/m2)
ALLOCATE(T%XWS_ROAD(1))      ! road water reservoir             (kg/m2)
ALLOCATE(T%XT_ROOF(1,KROOF_LAYER))       ! roof layer temperatures          (K)
ALLOCATE(T%XT_ROAD(1,KROAD_LAYER))       ! road layer temperatures          (K)
ALLOCATE(T%XT_WALL_A(1,KWALL_LAYER))     ! wall layer temperatures          (K)
ALLOCATE(T%XT_WALL_B(1,KWALL_LAYER))     ! wall layer temperatures          (K)
!
ALLOCATE(T%XAC_ROOF(1))      ! roof aerodynamic conductance     ()
ALLOCATE(T%XAC_ROAD(1))      ! road aerodynamic conductance     ()
ALLOCATE(T%XAC_WALL(1))      ! wall aerodynamic conductance     ()
ALLOCATE(T%XAC_TOP(1))       ! top  aerodynamic conductance     ()
ALLOCATE(T%XAC_ROOF_WAT(1))  ! water aerodynamic conductance    ()
ALLOCATE(T%XAC_ROAD_WAT(1))  ! water aerodynamic conductance    ()
!
ALLOCATE(T%XQSAT_ROOF(1))    ! humidity of saturation for roofs (kg/kg)
ALLOCATE(T%XQSAT_ROAD(1))    ! humidity of saturation for roads (kg/kg)
!
ALLOCATE(T%XDELT_ROOF(1))    ! humidity of saturation for roofs (-)
ALLOCATE(T%XDELT_ROAD(1))    ! humidity of saturation for roads (-)
!
! Semi-prognostic variables:
!
ALLOCATE(T%XT_CANYON(1))     ! canyon air temperature           (K)
ALLOCATE(T%XQ_CANYON(1))     ! canyon air specific humidity     (kg/kg)
!
!
! Prognostic snow:
!
ALLOCATE(T%TSNOW_ROOF%WSNOW(1,1))     ! snow (& liq. water) content      (kg/m2) 
ALLOCATE(T%TSNOW_ROOF%T    (1,1))     ! temperature                      (K)
ALLOCATE(T%TSNOW_ROOF%RHO  (1,1))     ! density                          (kg m-3)
ALLOCATE(T%TSNOW_ROOF%ALB  (1))     ! albedo                           (-)
ALLOCATE(T%TSNOW_ROOF%EMIS (1))     ! emissivity                       (-)
ALLOCATE(T%TSNOW_ROOF%TS   (1))     ! Surface temperature              (K)

ALLOCATE(T%TSNOW_ROAD%WSNOW(1,1))     ! snow (& liq. water) content      (kg/m2) 
ALLOCATE(T%TSNOW_ROAD%T    (1,1))     ! temperature                      (K)
ALLOCATE(T%TSNOW_ROAD%RHO  (1,1))     ! density                          (kg m-3)
ALLOCATE(T%TSNOW_ROAD%ALB  (1))     ! albedo                           (-)
ALLOCATE(T%TSNOW_ROAD%EMIS (1))     ! emissivity                       (-)
ALLOCATE(T%TSNOW_ROAD%TS   (1))     ! Surface temperature              (K)

!-------------------------------------------------------------------------------
!
! solar panels
!
ALLOCATE(TPN%XEMIS_PANEL(1))       ! emissivity of solar panels               (-)
ALLOCATE(TPN%XALB_PANEL(1))        ! albedo of soolar panels                  (-)
ALLOCATE(TPN%XEFF_PANEL(1))        ! Efficiency coefficient of solar panels   (-)
ALLOCATE(TPN%XFRAC_PANEL(1))       ! Fraction of SOLAR panels on roofs        (-)
ALLOCATE(TPN%XTHER_PRODC_DAY(1))   ! Present day integrated thermal production of energy    (J/m2 panel)


!-------------------------------------------------------------------------------
!
! Irrigation
!
ALLOCATE(TIR%XRD_START_MONTH(1))   ! roads : start month for irrigation (included)
ALLOCATE(TIR%XRD_END_MONTH(1))     ! roads : end   month for irrigation (included)
ALLOCATE(TIR%XRD_START_HOUR(1))    ! roads : start solar hour for irrigation (included, hour)
ALLOCATE(TIR%XRD_END_HOUR(1))      ! roads : end   solar hour for irrigation (excluded, hour)
ALLOCATE(TIR%XRD_24H_IRRIG(1))     ! roads : total irrigation over 24 hours (kg/m2)
!
!-------------------------------------------------------------------------------
!
! TEB diagnostics structure
!
ALLOCATE(DMT%XZ0_TOWN(1))  ! town roughness length
ALLOCATE(DMT%XQF_BLD(1))   ! domestic heating
ALLOCATE(DMT%XFLX_BLD(1)) ! heat flux from bld
ALLOCATE(DMT%XQF_TOWN(1))  ! total anthropogenic heat
ALLOCATE(DMT%XDQS_TOWN(1)) ! storage inside building
!
ALLOCATE(DMT%XH_WALL_A(1))   ! wall sensible heat flux          (W/m2)
ALLOCATE(DMT%XH_WALL_B(1))   ! wall sensible heat flux          (W/m2)
ALLOCATE(DMT%XH_ROOF(1))     ! roof sensible heat flux          (W/m2)
ALLOCATE(DMT%XH_ROAD(1))     ! road sensible heat flux          (W/m2)
ALLOCATE(DMT%XH_STRLROOF(1)) ! structural roof sens. heat flux  (W/m2)
ALLOCATE(DMT%XH_BLT(1))      ! built surf sensible heat flux    (W/m2)
ALLOCATE(DMT%XRN_WALL_A(1))  ! net radiation at wall            (W/m2)
ALLOCATE(DMT%XRN_WALL_B(1))  ! net radiation at wall            (W/m2)
ALLOCATE(DMT%XRN_ROOF(1))    ! net radiation at roof            (W/m2)
ALLOCATE(DMT%XRN_ROAD(1))    ! net radiation at road            (W/m2)
ALLOCATE(DMT%XRN_STRLROOF(1)) !net radiation at structural roofs(W/m2)
ALLOCATE(DMT%XRN_BLT(1))     ! net radiation at built surf      (W/m2)
ALLOCATE(DMT%XGFLUX_WALL_A(1)) !net wall conduction flux        (W/m2)
ALLOCATE(DMT%XGFLUX_WALL_B(1)) !net wall conduction flux        (W/m2)
ALLOCATE(DMT%XGFLUX_ROOF(1)) ! net roof conduction flux         (W/m2)                                         
ALLOCATE(DMT%XGFLUX_ROAD(1)) ! net road conduction flux         (W/m2)
ALLOCATE(DMT%XGFLUX_STRLROOF(1)) !net structural roof cond flux (W/m2)
ALLOCATE(DMT%XGFLUX_BLT(1))  ! net built surf conduction flux   (W/m2)
ALLOCATE(DMT%XLE_ROOF(1))    ! roof latent heat flux            (W/m2)
ALLOCATE(DMT%XLE_ROAD(1))    ! road latent heat flux            (W/m2)
ALLOCATE(DMT%XLE_STRLROOF(1)) !structural roof latent heat flux (W/m2)
ALLOCATE(DMT%XLE_BLT(1))     ! built surf latent heat flux      (W/m2)
!
ALLOCATE(DMT%XRUNOFF_TOWN(1))      ! aggregated water runoff for town      (kg/m2/s)
ALLOCATE(DMT%XIRRIG_GARDEN(1))     ! summer ground irrigation rate         (kg/m2/s)
ALLOCATE(DMT%XRUNOFF_ROAD(1))      ! water runoff for roads                (kg/m2/s)
ALLOCATE(DMT%XIRRIG_ROAD(1))       ! road man-made watering rate           (kg/m2/s)
ALLOCATE(DMT%XRUNOFF_ROOF(1))      ! aggregated water runoff for roofs     (kg/m2/s)
ALLOCATE(DMT%XRUNOFF_STRLROOF(1))  ! water runoff for structural roofs     (kg/m2/s)
ALLOCATE(DMT%XIRRIG_GREENROOF(1))  ! summer ground irrigation rate         (kg/m2/s)
!
ALLOCATE(DMT%XABS_SW_ROOF(1))      ! absorbed shortwave radiation over roofs
ALLOCATE(DMT%XABS_SW_SNOW_ROOF(1)) ! absorbed longwave radiation over roofs
ALLOCATE(DMT%XABS_LW_ROOF(1))      ! absorbed shortwave radiation over roofs
ALLOCATE(DMT%XABS_LW_SNOW_ROOF(1)) ! absorbed longwave radiation over roofs
ALLOCATE(DMT%XABS_SW_ROAD(1))      ! absorbed shortwave radiation over roads
ALLOCATE(DMT%XABS_SW_SNOW_ROAD(1)) ! absorbed longwave radiation over roads
ALLOCATE(DMT%XABS_LW_ROAD(1))      ! absorbed shortwave radiation over roads
ALLOCATE(DMT%XABS_LW_SNOW_ROAD(1)) ! absorbed longwave radiation over roads
ALLOCATE(DMT%XABS_SW_WALL_A(1))    ! absorbed shortwave radiation over walls
ALLOCATE(DMT%XABS_SW_WALL_B(1))    ! absorbed shortwave radiation over walls
ALLOCATE(DMT%XABS_LW_WALL_A(1))    ! absorbed shortwave radiation over walls
ALLOCATE(DMT%XABS_LW_WALL_B(1))    ! absorbed shortwave radiation over walls
ALLOCATE(DMT%XABS_SW_GARDEN(1))    ! absorbed shortwave radiation over green areas
ALLOCATE(DMT%XABS_LW_GARDEN(1))    ! absorbed shortwave radiation over green areas
ALLOCATE(DMT%XABS_SW_GREENROOF(1)) ! absorbed shortwave radiation over green roofs
ALLOCATE(DMT%XABS_LW_GREENROOF(1)) ! absorbed shortwave radiation over green roofs
ALLOCATE(DMT%XG_GREENROOF_ROOF(1)) ! Heat flux between green roof and structural roof
ALLOCATE(DMT%XABS_SW_PANEL(1))     ! absorbed shortwave radiation over solar panels
ALLOCATE(DMT%XABS_LW_PANEL(1))     ! absorbed longwave  radiation over solar panels
!
ALLOCATE(DMT%XRN_PANEL(1))         ! net radiation           over solar panels (W/m2)
ALLOCATE(DMT%XH_PANEL(1))          ! sensible heat flux      over solar panels (W/m2)
ALLOCATE(DMT%XTHER_PROD_PANEL(1))  ! thermal      production of   solar panels (W/m2 thermal panel)
ALLOCATE(DMT%XPHOT_PROD_PANEL(1))  ! photovoltaic production of   solar panels (W/m2 photovoltaic panel)
ALLOCATE(DMT%XPROD_PANEL(1))       !              production of   solar panels (W/m2 panel)
ALLOCATE(DMT%XTHER_PROD_BLD(1))    ! thermal      production of   solar panels (W/m2 bld)
ALLOCATE(DMT%XPHOT_PROD_BLD(1))    ! photovoltaic production of   solar panels (W/m2 bld)

ALLOCATE(DMT%XH_BLD_COOL(1))       ! Sensible cooling energy demand  
                                                     ! of the building [W m-2(bld)]
ALLOCATE(DMT%XT_BLD_COOL(1))       ! Total cooling energy demand  
                                                     ! of the building [W m-2(bld)]
ALLOCATE(DMT%XH_BLD_HEAT(1))       ! Heating energy demand       
                                                     ! of the building [W m-2(bld)]
ALLOCATE(DMT%XLE_BLD_COOL(1))      ! Latent cooling energy demand 
                                                     ! of the building [W m-2(bld)]
ALLOCATE(DMT%XLE_BLD_HEAT(1))      ! Latent heating energy demand 
                                                     ! of the building [W m-2(bld)]
ALLOCATE(DMT%XH_WASTE(1))          ! Sensible waste heat from HVAC system
                                                     ! [W m-2(tot)]
ALLOCATE(DMT%XLE_WASTE(1))         ! Latent waste heat from HVAC system
                                                     ! [W m-2(tot)]
ALLOCATE(DMT%XHVAC_COOL(1))        ! Energy consumption of the cooling system
                                                     ! [W m-2(bld)]
ALLOCATE(DMT%XHVAC_HEAT(1))        ! Energy consumption of the heating system
                                                     ! [W m-2(bld)]
ALLOCATE(DMT%XCAP_SYS(1))          ! Actual capacity of the cooling system
                                                     ! [W m-2(bld)] 
ALLOCATE(DMT%XM_SYS(1))            ! Actual HVAC mass flow rate 
                                                     ! [kg s-1 m-2(bld)]
ALLOCATE(DMT%XCOP(1))              ! COP of the cooling system
ALLOCATE(DMT%XQ_SYS(1))            ! Supply air specific humidity [kg kg-1]
ALLOCATE(DMT%XT_SYS(1))            ! Supply air temperature [K]
ALLOCATE(DMT%XTR_SW_WIN(1))        ! Solar radiation transmitted throught
                                                     ! windows [W m-2(bld)]
ALLOCATE(DMT%XFAN_POWER(1))        ! HVAC fan power
ALLOCATE(DMT%XABS_SW_WIN(1))       ! window absorbed shortwave radiation [W m-2] 
ALLOCATE(DMT%XABS_LW_WIN(1))       ! absorbed infrared rad. [W m-2]

ALLOCATE(DMT%XEMIT_LW_FAC(1))      ! LW flux emitted by the facade (W/m2 facade)
ALLOCATE(DMT%XEMIT_LW_GRND(1))     ! LW flux emitted by the ground (W/m2 ground = road + garden)
ALLOCATE(DMT%XT_RAD_IND(1))        !Indoor mean radiant temperature [K]
ALLOCATE(DMT%XREF_SW_GRND(1))      ! total solar rad reflected by ground
ALLOCATE(DMT%XREF_SW_FAC(1))       ! total solar rad reflected by facade
ALLOCATE(DMT%XHU_BLD(1))           !Indoor relative humidity
!
ALLOCATE(DMT%XTCOOL_TARGET(1))     ! current cooling setpoint of indoor air
ALLOCATE(DMT%XTHEAT_TARGET(1))     ! current heating setpoint of indoor air
ALLOCATE(DMT%XQIN(1))              ! current internal heat gains [W m-2(floor)]

! New fields
ALLOCATE(GDM%VD%ND%AL(1))
ALLOCATE(GDM%VD%NDE%AL(1))
ALLOCATE(GDM%VD%NDM%AL(1))
ALLOCATE(GDM%NPE%AL(1))

ALLOCATE(GRM%VD%ND%AL(1))
ALLOCATE(GRM%VD%NDE%AL(1))
ALLOCATE(GRM%VD%NDM%AL(1))
ALLOCATE(GRM%NPE%AL(1))

!-------------------------------------------------------------------------------
!
!
END SUBROUTINE ALLOC_TEB_STRUCT
