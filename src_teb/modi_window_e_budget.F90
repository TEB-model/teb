!auto_modi:spll_window_e_budget.D
MODULE MODI_WINDOW_E_BUDGET
INTERFACE
SUBROUTINE WINDOW_E_BUDGET(PU_WIN, PEMIS_WIN, PF_WIN_FLOOR, PF_WIN_WALL, PF_WIN_MASS, &
                           PLW_W_TO_WIN, PLW_R_TO_WIN, PLW_G_TO_WIN,                  &
                           PLW_NR_TO_WIN, PLW_S_TO_WIN, PRAD_ROOF_WIN, PRAD_WALL_WIN, &
                           PABS_SW_WIN, PLW_RAD, PAC_WALL,                 &
                           PRADHT_IN, PTI_BLD, PTS_MASS, PTS_FLOOR, PRHOA, PDN_ROAD,  &
                           PT_CANYON, PTS_WALL, PTS_ROAD, PTSNOW_ROAD,                &
                           PTS_GARDEN, PT_WIN1, PT_WIN2, PRAD_WIN_FLOOR,              &
                           PRAD_WIN_MASS, PCONV_WIN_BLD, PEMIT_LW_WIN, PABS_LW_WIN, PLOAD_IN_WIN)
IMPLICIT NONE
REAL, DIMENSION(:), INTENT(IN)    :: PU_WIN       ! window U-factor [W m-2 K-1]
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_WIN    ! window emissivity
REAL, DIMENSION(:), INTENT(IN)    :: PF_WIN_FLOOR ! window floor view factor
REAL, DIMENSION(:), INTENT(IN)    :: PF_WIN_WALL  ! window wall view factor
REAL, DIMENSION(:), INTENT(IN)    :: PF_WIN_MASS  ! window mass view factor
REAL, DIMENSION(:), INTENT(IN)    :: PLW_W_TO_WIN ! Radiative heat trasfer coeff window-wall outdoor
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_WIN ! Radiative heat trasfer coeff window-road 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_WIN ! Radiative heat trasfer coeff window-garden
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_WIN! Radiative heat trasfer coeff window-snow
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_WIN ! Radiative heat trasfer coeff window-sky 
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_ROOF_WIN  ! rad. fluxes from roof to win [W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_WALL_WIN  ! rad. fluxes from wall to win [W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WIN  ! window absorbed shortwave radiation [W m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD      ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WALL     ! aerodynamical conductance between wall and canyon
REAL, DIMENSION(:),   INTENT(IN)  :: PRADHT_IN     ! Indoor radiant heat transfer coefficient
REAL, DIMENSION(:), INTENT(IN)    :: PTI_BLD      ! inside building temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_MASS     ! surface mass temperature  [K]
REAL, DIMENSION(:), INTENT(IN)    :: PTS_FLOOR    ! floor layers temperatures [K]
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA        ! rho
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROAD     ! snow-covered fraction on roads
REAL, DIMENSION(:), INTENT(IN)    :: PT_CANYON    ! air canyon temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_WALL     ! wall outdoor surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_ROAD     ! road surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTSNOW_ROAD  ! road snow temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GARDEN   ! green area surface temperature
REAL, DIMENSION(:), INTENT(OUT)   :: PT_WIN1      ! outdoor window temperature [K]
REAL, DIMENSION(:), INTENT(INOUT) :: PT_WIN2      ! indoor window temperature [K]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WIN_FLOOR ! rad. fluxes from window to floor [W m-2(window)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WIN_MASS  ! rad. fluxes from window to mass [W m-2(window)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_WIN_BLD  ! conv. fluxes from window to bld [W m-2(window)]
REAL, DIMENSION(:), INTENT(OUT)   :: PEMIT_LW_WIN  ! Longwave radiation emitted by the window [W m-2(window)]
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_WIN   ! Longwave radiation absorbed by the window [W m-2(window)]
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_IN_WIN  ! solar + internal heat gain 
END SUBROUTINE WINDOW_E_BUDGET
END INTERFACE
END MODULE MODI_WINDOW_E_BUDGET
