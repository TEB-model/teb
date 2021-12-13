!auto_modi:spll_facade_e_budget.D
MODULE MODI_FACADE_E_BUDGET
INTERFACE
SUBROUTINE FACADE_E_BUDGET(TOP, T, B, DMT, PTSTEP, PDN_RD, PRHOA, PAC_WL, PAC_BLD, &
                           PLW_RAD, PPS, PEXNS, PT_CANYON, PTS_RD, PTSN_RD, PTS_GD, &
                           PTS_FL, PLW_WA_TO_WB, PLW_R_TO_WA, PLW_R_TO_WB,          &
                           PLW_G_TO_WA, PLW_G_TO_WB, PLW_S_TO_WA, PLW_S_TO_WB,      &
                           PLW_NR_TO_WA,PLW_NR_TO_WB, PLW_WIN_TO_WA, PLW_WIN_TO_WB, &
                           PLW_S_TO_WIN, PLW_WA_TO_WIN, PLW_WB_TO_WIN, PLW_R_TO_WIN,&
                           PLW_G_TO_WIN, PLW_NR_TO_WIN, PFLX_BLD_WL_A, PDQS_WL_A,   &
                           PFLX_BLD_WL_B, PDQS_WL_B, PEMIT_LW_FAC, PIMB_WL,         &
                           PRADHT_IN, PRAD_RF_WL, PRAD_RF_WIN, PRAD_WL_FL,          &
                           PRAD_WL_MA, PRAD_WIN_FL, PRAD_WIN_MA, PCONV_WL_BLD,      &
                           PCONV_WIN_BLD, PAC_WIN,  PLOAD_IN_WL, PLOAD_IN_WIN )
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_t
IMPLICIT NONE
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(BEM_t), INTENT(INOUT) :: B
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: DMT
REAL,               INTENT(IN)    :: PTSTEP       ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PDN_RD     ! snow-covered fraction on roads
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA        ! rho
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WL     ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(IN)    :: PAC_BLD      ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD      ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PPS          ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS        ! surface Exner function
REAL, DIMENSION(:), INTENT(IN)    :: PT_CANYON    ! air canyon temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_RD     ! road surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTSN_RD  ! road snow temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GD   ! green area surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_FL    ! floor layers temperatures [K]
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
REAL, DIMENSION(:), INTENT(OUT)   :: PFLX_BLD_WL_A! flux from bld to wall
REAL, DIMENSION(:), INTENT(INOUT) :: PDQS_WL_A    ! heat storage inside the wall 
REAL, DIMENSION(:), INTENT(OUT)   :: PFLX_BLD_WL_B! flux from bld to wall
REAL, DIMENSION(:), INTENT(INOUT) :: PDQS_WL_B    ! heat storage inside the wall 
REAL, DIMENSION(:), INTENT(OUT)   :: PEMIT_LW_FAC ! LW flux emitted by the facade (W/m2 of facade)
REAL, DIMENSION(:), INTENT(OUT)   :: PIMB_WL    ! wall residual energy imbalance 
REAL, DIMENSION(:),   INTENT(IN)  :: PRADHT_IN     ! Indoor radiant heat transfer coefficient
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_RF_WL ! rad. fluxes from roof to wall [W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_RF_WIN  ! rad. fluxes from roof to win [W m-2(roof)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WL_FL! rad. fluxes from wall to floor [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WL_MA ! rad. fluxes from wall to mass [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WIN_FL ! rad. fluxes from window to floor [W m-2(win)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WIN_MA  ! rad. fluxes from window to mass [W m-2(win)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_WL_BLD ! conv. fluxes from wall to bld [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_WIN_BLD  ! conv. fluxes from window to bld [W m-2(wind.)]
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WIN        ! window aerodynamic conductance
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_IN_WL  ! solar + inter. heat gains W/m2 [wall]
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_IN_WIN   ! solar + inter. heat gains W/m2 [win]
END SUBROUTINE FACADE_E_BUDGET
END INTERFACE
END MODULE MODI_FACADE_E_BUDGET
