!auto_modi:spll_facade_e_budget.D
MODULE MODI_FACADE_E_BUDGET
INTERFACE
SUBROUTINE FACADE_E_BUDGET(HWALL_OPT, HBEM,                                    &
                           PT_WALL_A, PT_WALL_B, PTSTEP,                       &
                           PHC_WALL, PTC_WALL, PD_WALL,                        &
                           PDN_ROAD, PRHOA, PAC_WALL,                          &
                           PAC_BLD,                                            &
                           PTI_BLD, PLW_RAD, PPS, PEXNS,                       &
                           PABS_SW_WALL_A, PABS_SW_WALL_B, PT_CANYON,          &
                           PTS_ROAD, PTSNOW_ROAD, PTS_GARDEN, PU_WIN, PGR,     &
                           PABS_SW_WIN, PF_WIN_MASS, PF_WIN_WALL,              &
                           PF_WIN_FLOOR, PTS_MASS, PTS_FLOOR, PEMIS_WALL,      &
                           PLW_WA_TO_WB,                                       &
                           PLW_R_TO_WA, PLW_R_TO_WB, PLW_G_TO_WA, PLW_G_TO_WB, &
                           PLW_S_TO_WA, PLW_S_TO_WB, PLW_NR_TO_WA,PLW_NR_TO_WB,&
                           PLW_WIN_TO_WA, PLW_WIN_TO_WB,                       &
                           PLW_S_TO_WIN, PLW_WA_TO_WIN, PLW_WB_TO_WIN,         &
                           PLW_R_TO_WIN, PLW_G_TO_WIN, PLW_NR_TO_WIN, PT_WIN1, &
                           PFLX_BLD_WALL_A, PDQS_WALL_A,                       &
                           PFLX_BLD_WALL_B, PDQS_WALL_B,                       &
                           PT_WIN2, PEMIT_LW_FAC,                              &
                           PABS_LW_WALL_A, PABS_LW_WALL_B, PABS_LW_WIN,        &
                           PH_WALL_A, PH_WALL_B,                               &
                           PIMB_WALL, PF_WALL_MASS, PF_WALL_FLOOR,             &
                           PF_WALL_WIN, PRADHT_IN, PRAD_ROOF_WALL,             &
                           PRAD_ROOF_WIN, PRAD_WALL_FLOOR, PRAD_WALL_MASS,     &
                           PRAD_WIN_FLOOR, PRAD_WIN_MASS, PCONV_WALL_BLD,      &
                           PCONV_WIN_BLD, PAC_WIN,  PLOAD_IN_WALL, PLOAD_IN_WIN )
IMPLICIT NONE
 CHARACTER(LEN=4), INTENT(IN)      :: HWALL_OPT     ! TEB option for walls
 CHARACTER(LEN=3), INTENT(IN)      :: HBEM         ! Building Energy model 'DEF' or 'BEM'
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_WALL_A  ! wall A layers temperatures
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_WALL_B  ! wall B layers temperatures
REAL,               INTENT(IN)    :: PTSTEP       ! time step
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_WALL     ! heat capacity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_WALL     ! thermal conductivity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_WALL      ! depth of wall layers
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROAD     ! snow-covered fraction on roads
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA        ! rho
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WALL     ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(IN)    :: PAC_BLD      ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(IN)    :: PTI_BLD      ! inside building temperature
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD      ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PPS          ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS        ! surface Exner function
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WALL_A! absorbed solar radiation
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WALL_B! absorbed solar radiation
REAL, DIMENSION(:), INTENT(IN)    :: PT_CANYON    ! air canyon temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_ROAD     ! road surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTSNOW_ROAD  ! road snow temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GARDEN   ! green area surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PU_WIN       ! window U-factor [W m-2 K-1]
REAL, DIMENSION(:), INTENT(IN)    :: PGR          ! glazing ratio
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WIN  ! window absorbed shortwave radiation [W m-2] 
REAL, DIMENSION(:),   INTENT(IN)  :: PF_WIN_MASS  ! window mass view factor
REAL, DIMENSION(:),   INTENT(IN)  :: PF_WIN_WALL  ! window wall view factor
REAL, DIMENSION(:),   INTENT(IN)  :: PF_WIN_FLOOR ! window floor view factor
REAL, DIMENSION(:), INTENT(IN)    :: PTS_MASS     ! surface mass temperature  [K]
REAL, DIMENSION(:), INTENT(IN)    :: PTS_FLOOR    ! floor layers temperatures [K]
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_WALL   ! wall emissivity
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_WB ! LW interactions wall  -> wall (both versions)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_WA  ! LW interactions road -> wall ; DEF formulation 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_WB  ! LW interactions road -> wall ; DEF formulation 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_WA  ! LW interactions garden -> wall ; DEF formulation
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_WB  ! LW interactions garden -> wall ; DEF formulation
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_WA  ! LW interactions sky   -> wall 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_WB  ! LW interactions sky   -> wall 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_WA ! LW interactions road(snow) -> wall 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_WB ! LW interactions road(snow) -> wall 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_WIN! Radiative heat trasfer coeff wall-window 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WB_TO_WIN! Radiative heat trasfer coeff wall-window 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_WIN ! Radiative heat trasfer coeff window-sky 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WIN_TO_WA! Radiative heat trasfer coeff window-wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WIN_TO_WB! Radiative heat trasfer coeff window-wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_WIN ! Radiative heat trasfer coeff window-road 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_WIN ! Radiative heat trasfer coeff window-garden
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_WIN! Radiative heat trasfer coeff window-snow
REAL, DIMENSION(:), INTENT(INOUT) :: PT_WIN1      ! outdoor window temperature [K]
REAL, DIMENSION(:), INTENT(OUT)   :: PFLX_BLD_WALL_A! flux from bld to wall
REAL, DIMENSION(:), INTENT(INOUT) :: PDQS_WALL_A    ! heat storage inside the wall 
REAL, DIMENSION(:), INTENT(OUT)   :: PFLX_BLD_WALL_B! flux from bld to wall
REAL, DIMENSION(:), INTENT(INOUT) :: PDQS_WALL_B    ! heat storage inside the wall 
REAL, DIMENSION(:), INTENT(INOUT) :: PT_WIN2      ! indoor window temperature [K]
REAL, DIMENSION(:), INTENT(OUT)   :: PEMIT_LW_FAC ! LW flux emitted by the facade (W/m2 of facade)
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_WALL_A ! absorbed infrared rad.
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_WALL_B ! absorbed infrared rad.
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_WIN  ! window absorbed shortwave radiation [W m-2] 
REAL, DIMENSION(:), INTENT(OUT)   :: PH_WALL_A    ! Sensible heat flux from wall to air [W/m2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PH_WALL_B    ! Sensible heat flux from wall to air [W/m2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PIMB_WALL    ! wall residual energy imbalance 
REAL, DIMENSION(:),   INTENT(IN)  :: PF_WALL_FLOOR ! View factor wall-floor
REAL, DIMENSION(:),   INTENT(IN)  :: PF_WALL_MASS  ! View factor wall-mass
REAL, DIMENSION(:),   INTENT(IN)  :: PF_WALL_WIN   ! View factor wall-win
REAL, DIMENSION(:),   INTENT(IN)  :: PRADHT_IN     ! Indoor radiant heat transfer coefficient
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_ROOF_WALL ! rad. fluxes from roof to wall [W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_ROOF_WIN  ! rad. fluxes from roof to win [W m-2(roof)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WALL_FLOOR! rad. fluxes from wall to floor [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WALL_MASS ! rad. fluxes from wall to mass [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WIN_FLOOR ! rad. fluxes from window to floor [W m-2(win)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WIN_MASS  ! rad. fluxes from window to mass [W m-2(win)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_WALL_BLD ! conv. fluxes from wall to bld [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_WIN_BLD  ! conv. fluxes from window to bld [W m-2(wind.)]
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WIN        ! window aerodynamic conductance
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_IN_WALL  ! solar + inter. heat gains W/m2 [wall]
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_IN_WIN   ! solar + inter. heat gains W/m2 [win]
END SUBROUTINE FACADE_E_BUDGET
END INTERFACE
END MODULE MODI_FACADE_E_BUDGET
