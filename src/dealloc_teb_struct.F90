!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL licence
!SFX_LIC version 2.1. See Licence_CeCILL_V2.1-en.txt and Licence_CeCILL_V2.1-fr.txt  
!SFX_LIC for details.
!   ##########################################################################
    SUBROUTINE DEALLOC_TEB_STRUCT(TOP,T,BOP,B,TPN,TIR,DMT,GDM,GRM)
!   ##########################################################################
!
!!****  *DEALLOC_TEB_STRUCT*
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
! BEM structure
!
!
! Floor parameters
!
DEALLOCATE(B%XHC_FLOOR) ! floor layers heat capacity        (J/K/m3)
DEALLOCATE(B%XTC_FLOOR) ! floor layers thermal conductivity (W/K/m)
DEALLOCATE(B%XD_FLOOR ) ! depth of floor layers             (m)
!
! HVAC parameters
!
DEALLOCATE(B%XTCOOL_TARGET) ! cooling setpoint of indoor air
DEALLOCATE(B%XTHEAT_TARGET) ! heating setpoint of indoor air
DEALLOCATE(B%XF_WASTE_CAN ) ! fraction of waste heat released into the canyon
DEALLOCATE(B%XEFF_HEAT    ) ! efficiency of the heating system
!
! Indoor parameters
!
DEALLOCATE(B%XTI_BLD)                    ! building interior temperature    (K)
DEALLOCATE(B%XT_FLOOR ) ! floor layer temperatures         (K)
DEALLOCATE(B%XT_MASS  ) ! Air cooled building internal th. mass temperature (K)
!
DEALLOCATE(B%XQIN)           ! internal heat gains [W m-2(floor)]
DEALLOCATE(B%XQIN_FRAD)      ! radiant fraction of internal heat gains
DEALLOCATE(B%XSHGC)          ! solar heat gain coef. of windows
DEALLOCATE(B%XSHGC_SH)       ! solar heat gain coef. of windows + shading
DEALLOCATE(B%XU_WIN)         ! window U-factor [K m W-2]
DEALLOCATE(B%XTRAN_WIN)      ! window transmittance (-)
DEALLOCATE(B%XGR)            ! glazing ratio
DEALLOCATE(B%XFLOOR_HEIGHT)  ! building floor height [m]
DEALLOCATE(B%XINF)           ! infiltration/ventilation flow rate [AC/H]
!
! New parameters
!
DEALLOCATE(B%XF_WATER_COND) ! fraction of evaporation for condensers (cooling system)
DEALLOCATE(B%XAUX_MAX)      ! Auxiliar variable for autosize calcs
DEALLOCATE(B%XQIN_FLAT)     ! Latent franction of internal heat gains
DEALLOCATE(B%XHR_TARGET)    ! Relative humidity setpoint
DEALLOCATE(B%XT_WIN2)       ! Indoor window temperature [K]
DEALLOCATE(B%XQI_BLD)       ! Indoor air specific humidity [kg kg-1]
DEALLOCATE(B%XV_VENT)       ! Ventilation flow rate [AC/H]
DEALLOCATE(B%XCAP_SYS_HEAT) ! Capacity of the heating system 
                             ! [W m-2(bld)]
DEALLOCATE(B%XCAP_SYS_RAT)  ! Rated capacity of the cooling system
                             ! [W m-2(bld)]
DEALLOCATE(B%XT_ADP)        ! Apparatus dewpoint temperature of the
                             ! cooling coil [K]
DEALLOCATE(B%XM_SYS_RAT)    ! Rated HVAC mass flow rate 
                             ! [kg s-1 m-2(bld)]
DEALLOCATE(B%XCOP_RAT)      ! Rated COP of the cooling system
DEALLOCATE(B%XT_WIN1)       ! outdoor window temperature [K]
DEALLOCATE(B%XALB_WIN)      ! window albedo
DEALLOCATE(B%XABS_WIN)      ! window absortance
DEALLOCATE(B%XT_SIZE_MAX)   ! Maximum outdoor air temperature for
                             ! HVAC sizing [K]
DEALLOCATE(B%XT_SIZE_MIN)   ! Minimum outdoor air temperature for
                             ! HVAC sizing [K]
DEALLOCATE(B%XUGG_WIN)      ! Window glass-to-glass U-factor [K m W-2]
DEALLOCATE(B%LSHADE)        ! flag to activate shading devices -> LOGICAL in the code
DEALLOCATE(B%XSHADE)        ! flag to activate shading devices -> REAL for i/o 0. or 1.
DEALLOCATE(B%CNATVENT)      ! flag to activate natural ventilation 'NONE', 'MANU', 'AUTO'
DEALLOCATE(B%XNATVENT)      ! flag to describe surventilation system for i/o 
                             ! 0 for NONE, 1 for MANU and 2 for AUTO
DEALLOCATE(B%LSHAD_DAY)     !Has shading been necessary this day ?
DEALLOCATE(B%LNATVENT_NIGHT)!Has nocturnal surventilation been necessary and possible this night ?
  !
  !indoor relative surfaces and view factors
DEALLOCATE(B%XN_FLOOR)        ! Number of floors     
DEALLOCATE(B%XGLAZ_O_BLD)     ! Window area [m2_win/m2_bld]
DEALLOCATE(B%XMASS_O_BLD)     ! Mass area [m2_mass/m2_bld]
DEALLOCATE(B%XFLOOR_HW_RATIO) ! H/W ratio of 1 floor level
DEALLOCATE(B%XF_FLOOR_MASS)   ! View factor floor-mass
DEALLOCATE(B%XF_FLOOR_WALL)   ! View factor floor-wall
DEALLOCATE(B%XF_FLOOR_WIN)    ! View factor floor-window
DEALLOCATE(B%XF_FLOOR_ROOF)   ! View factor floor-roof
DEALLOCATE(B%XF_WALL_FLOOR)   ! View factor wall-floor
DEALLOCATE(B%XF_WALL_MASS)    ! View factor wall-mass
DEALLOCATE(B%XF_WALL_WIN)     ! View factor wall-win
DEALLOCATE(B%XF_WIN_FLOOR)    ! View factor win-floor
DEALLOCATE(B%XF_WIN_MASS)     ! View factor win-mass
DEALLOCATE(B%XF_WIN_WALL)     ! View factor win-wall
DEALLOCATE(B%XF_WIN_WIN)      ! indoor View factor win-win
DEALLOCATE(B%XF_MASS_FLOOR)   ! View factor mass-floor
DEALLOCATE(B%XF_MASS_WALL)    ! View factor mass-wall
DEALLOCATE(B%XF_MASS_WIN)     ! View factor mass-window
!
!-------------------------------------------------------------------------------
!
! TEB structure
!
! Geometric Parameters:
!
DEALLOCATE(T%XROAD_DIR)     ! Road direction (deg from North, clockwise)
DEALLOCATE(T%XGARDEN)       ! fraction of veg in the streets   (-)
DEALLOCATE(T%XGREENROOF)    ! fraction of greenroofs on roofs  (-)
DEALLOCATE(T%XBLD)          ! fraction of buildings            (-)
DEALLOCATE(T%XROAD)         ! fraction of roads                (-)
DEALLOCATE(T%XCAN_HW_RATIO) ! canyon    h/W                    (-)
DEALLOCATE(T%XBLD_HEIGHT)   ! buildings height 'h'             (m)
DEALLOCATE(T%XWALL_O_HOR)   ! wall surf. / hor. surf.          (-)
DEALLOCATE(T%XROAD_O_GRND)  ! road surf. / (road + garden surf.) (-)
DEALLOCATE(T%XGARDEN_O_GRND)! gard. surf. / (road + garden surf.)(-)
DEALLOCATE(T%XWALL_O_GRND)  ! wall surf. / (road + garden surf.) (-)
DEALLOCATE(T%XWALL_O_BLD)   ! wall surf. / bld surf. (-)
DEALLOCATE(T%XZ0_TOWN)      ! roughness length for momentum    (m)
DEALLOCATE(T%XSVF_ROAD)     ! road sky view factor             (-)
DEALLOCATE(T%XSVF_GARDEN)   ! green area sky view factor       (-)
DEALLOCATE(T%XSVF_WALL)     ! wall sky view factor             (-)
!
! Roof parameters
!
DEALLOCATE(T%XALB_ROOF)     ! roof albedo                      (-)
DEALLOCATE(T%XEMIS_ROOF)    ! roof emissivity                  (-)
DEALLOCATE(T%XHC_ROOF)      ! roof layers heat capacity        (J/K/m3)
DEALLOCATE(T%XTC_ROOF)      ! roof layers thermal conductivity (W/K/m)
DEALLOCATE(T%XD_ROOF)       ! depth of roof layers             (m)
DEALLOCATE(T%XROUGH_ROOF)   ! roof roughness coef
!
!
! Road parameters
!
DEALLOCATE(T%XALB_ROAD)     ! road albedo                      (-)
DEALLOCATE(T%XEMIS_ROAD)    ! road emissivity                  (-)
DEALLOCATE(T%XHC_ROAD)      ! road layers heat capacity        (J/K/m3)
DEALLOCATE(T%XTC_ROAD)      ! road layers thermal conductivity (W/K/m)
DEALLOCATE(T%XD_ROAD)       ! depth of road layers             (m)
!
! Wall parameters
!
DEALLOCATE(T%XALB_WALL)     ! wall albedo                      (-)
DEALLOCATE(T%XEMIS_WALL)    ! wall emissivity                  (-)
DEALLOCATE(T%XHC_WALL)      ! wall layers heat capacity        (J/K/m3)
DEALLOCATE(T%XTC_WALL)      ! wall layers thermal conductivity (W/K/m)
DEALLOCATE(T%XD_WALL)       ! depth of wall layers             (m)
DEALLOCATE(T%XROUGH_WALL)   ! wall roughness coef
!
! Building's use type
!
DEALLOCATE(T%XRESIDENTIAL)  ! fraction of Residential use      (-)
  
!
! anthropogenic fluxes
!
DEALLOCATE(T%XH_TRAFFIC)    ! anthropogenic sensible
!                                                  ! heat fluxes due to traffic       (W/m2)
DEALLOCATE(T%XLE_TRAFFIC)   ! anthropogenic latent
!                                                  ! heat fluxes due to traffic       (W/m2)
DEALLOCATE(T%XH_INDUSTRY)   ! anthropogenic sensible                   
!                                                  ! heat fluxes due to factories     (W/m2)
DEALLOCATE(T%XLE_INDUSTRY)  ! anthropogenic latent
!                                                  ! heat fluxes due to factories     (W/m2)
!
! temperatures for boundary conditions
!
DEALLOCATE(T%XTI_ROAD)      ! road interior temperature        (K)
!
! Prognostic variables:
!
DEALLOCATE(T%XWS_ROOF)      ! roof water reservoir             (kg/m2)
DEALLOCATE(T%XWS_ROAD)      ! road water reservoir             (kg/m2)
DEALLOCATE(T%XT_ROOF)       ! roof layer temperatures          (K)
DEALLOCATE(T%XT_ROAD)       ! road layer temperatures          (K)
DEALLOCATE(T%XT_WALL_A)     ! wall layer temperatures          (K)
DEALLOCATE(T%XT_WALL_B)     ! wall layer temperatures          (K)
!
DEALLOCATE(T%XAC_ROOF)      ! roof aerodynamic conductance     ()
DEALLOCATE(T%XAC_ROAD)      ! road aerodynamic conductance     ()
DEALLOCATE(T%XAC_WALL)      ! wall aerodynamic conductance     ()
DEALLOCATE(T%XAC_TOP)       ! top  aerodynamic conductance     ()
DEALLOCATE(T%XAC_ROOF_WAT)  ! water aerodynamic conductance    ()
DEALLOCATE(T%XAC_ROAD_WAT)  ! water aerodynamic conductance    ()
!
DEALLOCATE(T%XQSAT_ROOF)    ! humidity of saturation for roofs (kg/kg)
DEALLOCATE(T%XQSAT_ROAD)    ! humidity of saturation for roads (kg/kg)
!
DEALLOCATE(T%XDELT_ROOF)    ! humidity of saturation for roofs (-)
DEALLOCATE(T%XDELT_ROAD)    ! humidity of saturation for roads (-)
!
! Semi-prognostic variables:
!
DEALLOCATE(T%XT_CANYON)     ! canyon air temperature           (K)
DEALLOCATE(T%XQ_CANYON)     ! canyon air specific humidity     (kg/kg)
!
!
! Prognostic snow:
!
DEALLOCATE(T%TSNOW_ROOF%WSNOW)     ! snow (& liq. water) content      (kg/m2) 
DEALLOCATE(T%TSNOW_ROOF%T    )     ! temperature                      (K)
DEALLOCATE(T%TSNOW_ROOF%RHO  )     ! density                          (kg m-3)
DEALLOCATE(T%TSNOW_ROOF%ALB  )     ! albedo                           (-)
DEALLOCATE(T%TSNOW_ROOF%EMIS )     ! emissivity                       (-)
DEALLOCATE(T%TSNOW_ROOF%TS   )     ! Surface temperature              (K)

DEALLOCATE(T%TSNOW_ROAD%WSNOW)     ! snow (& liq. water) content      (kg/m2) 
DEALLOCATE(T%TSNOW_ROAD%T    )     ! temperature                      (K)
DEALLOCATE(T%TSNOW_ROAD%RHO  )     ! density                          (kg m-3)
DEALLOCATE(T%TSNOW_ROAD%ALB  )     ! albedo                           (-)
DEALLOCATE(T%TSNOW_ROAD%EMIS )     ! emissivity                       (-)
DEALLOCATE(T%TSNOW_ROAD%TS   )     ! Surface temperature              (K)

!-------------------------------------------------------------------------------
!
! solar panels
!
DEALLOCATE(TPN%XEMIS_PANEL)       ! emissivity of solar panels               (-)
DEALLOCATE(TPN%XALB_PANEL)        ! albedo of soolar panels                  (-)
DEALLOCATE(TPN%XEFF_PANEL)        ! Efficiency coefficient of solar panels   (-)
DEALLOCATE(TPN%XFRAC_PANEL)       ! Fraction of SOLAR panels on roofs        (-)
DEALLOCATE(TPN%XTHER_PRODC_DAY)   ! Present day integrated thermal production of energy    (J/m2 panel)


!-------------------------------------------------------------------------------
!
! Irrigation
!
DEALLOCATE(TIR%XRD_START_MONTH)   ! roads : start month for irrigation (included)
DEALLOCATE(TIR%XRD_END_MONTH)     ! roads : end   month for irrigation (included)
DEALLOCATE(TIR%XRD_START_HOUR)    ! roads : start solar hour for irrigation (included, hour)
DEALLOCATE(TIR%XRD_END_HOUR)      ! roads : end   solar hour for irrigation (excluded, hour)
DEALLOCATE(TIR%XRD_24H_IRRIG)     ! roads : total irrigation over 24 hours (kg/m2)
!
!-------------------------------------------------------------------------------
!
! TEB diagnostics structure
!
DEALLOCATE(DMT%XZ0_TOWN)  ! town roughness length
DEALLOCATE(DMT%XQF_BLD)   ! domestic heating
DEALLOCATE(DMT%XFLX_BLD) ! heat flux from bld
DEALLOCATE(DMT%XQF_TOWN)  ! total anthropogenic heat
DEALLOCATE(DMT%XDQS_TOWN) ! storage inside building
!
DEALLOCATE(DMT%XH_WALL_A)   ! wall sensible heat flux          (W/m2)
DEALLOCATE(DMT%XH_WALL_B)   ! wall sensible heat flux          (W/m2)
DEALLOCATE(DMT%XH_ROOF)     ! roof sensible heat flux          (W/m2)
DEALLOCATE(DMT%XH_ROAD)     ! road sensible heat flux          (W/m2)
DEALLOCATE(DMT%XH_STRLROOF) ! structural roof sens. heat flux  (W/m2)
DEALLOCATE(DMT%XH_BLT)      ! built surf sensible heat flux    (W/m2)
DEALLOCATE(DMT%XRN_WALL_A)  ! net radiation at wall            (W/m2)
DEALLOCATE(DMT%XRN_WALL_B)  ! net radiation at wall            (W/m2)
DEALLOCATE(DMT%XRN_ROOF)    ! net radiation at roof            (W/m2)
DEALLOCATE(DMT%XRN_ROAD)    ! net radiation at road            (W/m2)
DEALLOCATE(DMT%XRN_STRLROOF) !net radiation at structural roofs(W/m2)
DEALLOCATE(DMT%XRN_BLT)     ! net radiation at built surf      (W/m2)
DEALLOCATE(DMT%XGFLUX_WALL_A) !net wall conduction flux        (W/m2)
DEALLOCATE(DMT%XGFLUX_WALL_B) !net wall conduction flux        (W/m2)
DEALLOCATE(DMT%XGFLUX_ROOF) ! net roof conduction flux         (W/m2)                                         
DEALLOCATE(DMT%XGFLUX_ROAD) ! net road conduction flux         (W/m2)
DEALLOCATE(DMT%XGFLUX_STRLROOF) !net structural roof cond flux (W/m2)
DEALLOCATE(DMT%XGFLUX_BLT)  ! net built surf conduction flux   (W/m2)
DEALLOCATE(DMT%XLE_ROOF)    ! roof latent heat flux            (W/m2)
DEALLOCATE(DMT%XLE_ROAD)    ! road latent heat flux            (W/m2)
DEALLOCATE(DMT%XLE_STRLROOF) !structural roof latent heat flux (W/m2)
DEALLOCATE(DMT%XLE_BLT)     ! built surf latent heat flux      (W/m2)
!
DEALLOCATE(DMT%XRUNOFF_TOWN)      ! aggregated water runoff for town      (kg/m2/s)
DEALLOCATE(DMT%XIRRIG_GARDEN)     ! summer ground irrigation rate         (kg/m2/s)
DEALLOCATE(DMT%XRUNOFF_ROAD)      ! water runoff for roads                (kg/m2/s)
DEALLOCATE(DMT%XIRRIG_ROAD)       ! road man-made watering rate           (kg/m2/s)
DEALLOCATE(DMT%XRUNOFF_ROOF)      ! aggregated water runoff for roofs     (kg/m2/s)
DEALLOCATE(DMT%XRUNOFF_STRLROOF)  ! water runoff for structural roofs     (kg/m2/s)
DEALLOCATE(DMT%XIRRIG_GREENROOF)  ! summer ground irrigation rate         (kg/m2/s)
!
DEALLOCATE(DMT%XABS_SW_ROOF)      ! absorbed shortwave radiation over roofs
DEALLOCATE(DMT%XABS_SW_SNOW_ROOF) ! absorbed longwave radiation over roofs
DEALLOCATE(DMT%XABS_LW_ROOF)      ! absorbed shortwave radiation over roofs
DEALLOCATE(DMT%XABS_LW_SNOW_ROOF) ! absorbed longwave radiation over roofs
DEALLOCATE(DMT%XABS_SW_ROAD)      ! absorbed shortwave radiation over roads
DEALLOCATE(DMT%XABS_SW_SNOW_ROAD) ! absorbed longwave radiation over roads
DEALLOCATE(DMT%XABS_LW_ROAD)      ! absorbed shortwave radiation over roads
DEALLOCATE(DMT%XABS_LW_SNOW_ROAD) ! absorbed longwave radiation over roads
DEALLOCATE(DMT%XABS_SW_WALL_A)    ! absorbed shortwave radiation over walls
DEALLOCATE(DMT%XABS_SW_WALL_B)    ! absorbed shortwave radiation over walls
DEALLOCATE(DMT%XABS_LW_WALL_A)    ! absorbed shortwave radiation over walls
DEALLOCATE(DMT%XABS_LW_WALL_B)    ! absorbed shortwave radiation over walls
DEALLOCATE(DMT%XABS_SW_GARDEN)    ! absorbed shortwave radiation over green areas
DEALLOCATE(DMT%XABS_LW_GARDEN)    ! absorbed shortwave radiation over green areas
DEALLOCATE(DMT%XABS_SW_GREENROOF) ! absorbed shortwave radiation over green roofs
DEALLOCATE(DMT%XABS_LW_GREENROOF) ! absorbed shortwave radiation over green roofs
DEALLOCATE(DMT%XG_GREENROOF_ROOF) ! Heat flux between green roof and structural roof
DEALLOCATE(DMT%XABS_SW_PANEL)     ! absorbed shortwave radiation over solar panels
DEALLOCATE(DMT%XABS_LW_PANEL)     ! absorbed longwave  radiation over solar panels
!
DEALLOCATE(DMT%XRN_PANEL)         ! net radiation           over solar panels (W/m2)
DEALLOCATE(DMT%XH_PANEL)          ! sensible heat flux      over solar panels (W/m2)
DEALLOCATE(DMT%XTHER_PROD_PANEL)  ! thermal      production of   solar panels (W/m2 thermal panel)
DEALLOCATE(DMT%XPHOT_PROD_PANEL)  ! photovoltaic production of   solar panels (W/m2 photovoltaic panel)
DEALLOCATE(DMT%XPROD_PANEL)       !              production of   solar panels (W/m2 panel)
DEALLOCATE(DMT%XTHER_PROD_BLD)    ! thermal      production of   solar panels (W/m2 bld)
DEALLOCATE(DMT%XPHOT_PROD_BLD)    ! photovoltaic production of   solar panels (W/m2 bld)

DEALLOCATE(DMT%XH_BLD_COOL)       ! Sensible cooling energy demand  
                                                     ! of the building [W m-2(bld)]
DEALLOCATE(DMT%XT_BLD_COOL)       ! Total cooling energy demand  
                                                     ! of the building [W m-2(bld)]
DEALLOCATE(DMT%XH_BLD_HEAT)       ! Heating energy demand       
                                                     ! of the building [W m-2(bld)]
DEALLOCATE(DMT%XLE_BLD_COOL)      ! Latent cooling energy demand 
                                                     ! of the building [W m-2(bld)]
DEALLOCATE(DMT%XLE_BLD_HEAT)      ! Latent heating energy demand 
                                                     ! of the building [W m-2(bld)]
DEALLOCATE(DMT%XH_WASTE)          ! Sensible waste heat from HVAC system
                                                     ! [W m-2(tot)]
DEALLOCATE(DMT%XLE_WASTE)         ! Latent waste heat from HVAC system
                                                     ! [W m-2(tot)]
DEALLOCATE(DMT%XHVAC_COOL)        ! Energy consumption of the cooling system
                                                     ! [W m-2(bld)]
DEALLOCATE(DMT%XHVAC_HEAT)        ! Energy consumption of the heating system
                                                     ! [W m-2(bld)]
DEALLOCATE(DMT%XCAP_SYS)          ! Actual capacity of the cooling system
                                                     ! [W m-2(bld)] 
DEALLOCATE(DMT%XM_SYS)            ! Actual HVAC mass flow rate 
                                                     ! [kg s-1 m-2(bld)]
DEALLOCATE(DMT%XCOP)              ! COP of the cooling system
DEALLOCATE(DMT%XQ_SYS)            ! Supply air specific humidity [kg kg-1]
DEALLOCATE(DMT%XT_SYS)            ! Supply air temperature [K]
DEALLOCATE(DMT%XTR_SW_WIN)        ! Solar radiation transmitted throught
                                                     ! windows [W m-2(bld)]
DEALLOCATE(DMT%XFAN_POWER)        ! HVAC fan power
DEALLOCATE(DMT%XABS_SW_WIN)       ! window absorbed shortwave radiation [W m-2] 
DEALLOCATE(DMT%XABS_LW_WIN)       ! absorbed infrared rad. [W m-2]

DEALLOCATE(DMT%XEMIT_LW_FAC)      ! LW flux emitted by the facade (W/m2 facade)
DEALLOCATE(DMT%XEMIT_LW_GRND)     ! LW flux emitted by the ground (W/m2 ground = road + garden)
DEALLOCATE(DMT%XT_RAD_IND)        !Indoor mean radiant temperature [K]
DEALLOCATE(DMT%XREF_SW_GRND)      ! total solar rad reflected by ground
DEALLOCATE(DMT%XREF_SW_FAC)       ! total solar rad reflected by facade
DEALLOCATE(DMT%XHU_BLD)           !Indoor relative humidity
!
DEALLOCATE(DMT%XTCOOL_TARGET)     ! current cooling setpoint of indoor air
DEALLOCATE(DMT%XTHEAT_TARGET)     ! current heating setpoint of indoor air
DEALLOCATE(DMT%XQIN)              ! current internal heat gains [W m-2(floor)]

! New fields
DEALLOCATE(GDM%VD%ND%AL)
DEALLOCATE(GDM%VD%NDE%AL)
DEALLOCATE(GDM%VD%NDM%AL)
DEALLOCATE(GDM%NPE%AL)

DEALLOCATE(GRM%VD%ND%AL)
DEALLOCATE(GRM%VD%NDE%AL)
DEALLOCATE(GRM%VD%NDM%AL)
DEALLOCATE(GRM%NPE%AL)

!-------------------------------------------------------------------------------
!
!
END SUBROUTINE DEALLOC_TEB_STRUCT
