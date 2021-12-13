!auto_modi:spll_window_e_budget.D
MODULE MODI_WINDOW_E_BUDGET
INTERFACE
SUBROUTINE WINDOW_E_BUDGET(B, PEMIS_WIN, PLW_W_TO_WIN, PLW_R_TO_WIN, PLW_G_TO_WIN, &
                           PLW_NR_TO_WIN, PLW_S_TO_WIN, PRAD_RF_WIN, PRAD_WL_WIN,  &
                           PABS_SW_WIN, PLW_RAD, PAC_WL, PRADHT_IN, PTS_FL, PRHOA, &
                           PDN_RD, PT_CANYON, PTS_WL, PTS_RD, PTSN_RD, PTS_GD,     &
                           PRAD_WIN_FL, PRAD_WIN_MA, PCONV_WIN_BLD, PEMIT_LW_WIN,  &
                           PABS_LW_WIN, PLOAD_IN_WIN)
USE MODD_BEM_n, ONLY : BEM_t
IMPLICIT NONE
TYPE(BEM_t), INTENT(INOUT) :: B
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_WIN    ! window emissivity
REAL, DIMENSION(:), INTENT(IN)    :: PLW_W_TO_WIN ! Radiative heat trasfer coeff window-wall outdoor
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_WIN ! Radiative heat trasfer coeff window-road 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_WIN ! Radiative heat trasfer coeff window-garden
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_WIN! Radiative heat trasfer coeff window-snow
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_WIN ! Radiative heat trasfer coeff window-sky 
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_RF_WIN  ! rad. fluxes from roof to win [W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_WL_WIN  ! rad. fluxes from wall to win [W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WIN  ! window absorbed shortwave radiation [W m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD      ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WL     ! aerodynamical conductance between wall and canyon
REAL, DIMENSION(:),   INTENT(IN)  :: PRADHT_IN     ! Indoor radiant heat transfer coefficient
REAL, DIMENSION(:), INTENT(IN)    :: PTS_FL    ! floor layers temperatures [K]
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA        ! rho
REAL, DIMENSION(:), INTENT(IN)    :: PDN_RD     ! snow-covered fraction on roads
REAL, DIMENSION(:), INTENT(IN)    :: PT_CANYON    ! air canyon temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_WL     ! wall outdoor surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_RD     ! road surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTSN_RD  ! road snow temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GD   ! green area surface temperature
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WIN_FL ! rad. fluxes from window to floor [W m-2(window)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WIN_MA  ! rad. fluxes from window to mass [W m-2(window)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_WIN_BLD  ! conv. fluxes from window to bld [W m-2(window)]
REAL, DIMENSION(:), INTENT(OUT)   :: PEMIT_LW_WIN  ! Longwave radiation emitted by the window [W m-2(window)]
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_WIN   ! Longwave radiation absorbed by the window [W m-2(window)]
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_IN_WIN  ! solar + internal heat gain 
END SUBROUTINE WINDOW_E_BUDGET
END INTERFACE
END MODULE MODI_WINDOW_E_BUDGET
