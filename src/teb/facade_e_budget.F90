!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
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
!
!   ##########################################################################
!
!!****  *FACADE_E_BUDGET*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the evolution of wall and window temperature from
!     wall_layer_e_budget et window_e_budget
!         
!     
!!**  METHOD
!     ------
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    MODD_CST
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!      G. Pigeon           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/11/2011 
!!      G. Pigeon     /09/2012: new conv. coef for indoor/outdoor 
!!      G. Pigeon     /10/2012: new arg. : solar heat load for indoor
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_t
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_CSTS, ONLY : XCPD
!
USE MODI_WALL_LAYER_E_BUDGET
USE MODI_WINDOW_E_BUDGET
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(BEM_t), INTENT(INOUT) :: B
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: DMT
!
REAL,               INTENT(IN)    :: PTSTEP       ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PDN_RD     ! snow-covered fraction on roads
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA        ! rho
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WL     ! aerodynamical conductance
!                                                 ! between wall and canyon
REAL, DIMENSION(:), INTENT(IN)    :: PAC_BLD      ! aerodynamical conductance
                                                  ! inside the building itself
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD      ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PPS          ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS        ! surface Exner function
REAL, DIMENSION(:), INTENT(IN)    :: PT_CANYON    ! air canyon temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_RD     ! road surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTSN_RD  ! road snow temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GD   ! green area surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_FL    ! floor layers temperatures [K]
!
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
                                                  ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WB_TO_WIN! Radiative heat trasfer coeff wall-window 
                                                  ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_WIN ! Radiative heat trasfer coeff window-sky 
                                                  ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WIN_TO_WA! Radiative heat trasfer coeff window-wall
                                                  ! [W K-1 m-2] 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WIN_TO_WB! Radiative heat trasfer coeff window-wall
                                                  ! [W K-1 m-2] 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_WIN ! Radiative heat trasfer coeff window-road 
                                                  ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_WIN ! Radiative heat trasfer coeff window-garden
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_WIN! Radiative heat trasfer coeff window-snow
!
REAL, DIMENSION(:), INTENT(OUT)   :: PFLX_BLD_WL_A! flux from bld to wall
REAL, DIMENSION(:), INTENT(INOUT) :: PDQS_WL_A    ! heat storage inside the wall 
REAL, DIMENSION(:), INTENT(OUT)   :: PFLX_BLD_WL_B! flux from bld to wall
REAL, DIMENSION(:), INTENT(INOUT) :: PDQS_WL_B    ! heat storage inside the wall 
REAL, DIMENSION(:), INTENT(OUT)   :: PEMIT_LW_FAC ! LW flux emitted by the facade (W/m2 of facade)
REAL, DIMENSION(:), INTENT(OUT)   :: PIMB_WL    ! wall residual energy imbalance 
                                                  ! for verification [W m-2]
REAL, DIMENSION(:),   INTENT(IN)  :: PRADHT_IN     ! Indoor radiant heat transfer coefficient
                                                    ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_RF_WL ! rad. fluxes from roof to wall [W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_RF_WIN  ! rad. fluxes from roof to win [W m-2(roof)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WL_FL! rad. fluxes from wall to floor [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WL_MA ! rad. fluxes from wall to mass [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WIN_FL ! rad. fluxes from window to floor [W m-2(win)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WIN_MA  ! rad. fluxes from window to mass [W m-2(win)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_WL_BLD ! conv. fluxes from wall to bld [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_WIN_BLD  ! conv. fluxes from window to bld [W m-2(wind.)]
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WIN        ! window aerodynamic conductance
!
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_IN_WL  ! solar + inter. heat gains W/m2 [wall]
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_IN_WIN   ! solar + inter. heat gains W/m2 [win]
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZRAD_WL_WIN   ! rad. fluxes between averaged wall and win [W m-2(wall)]
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZRAD_WL_A_WIN ! rad. fluxes between one wall and win [W m-2(wall)]
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZEMIS_WIN     ! window emissivity
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZEMIT_LW_WIN  ! Longwave radiation emitted by the window [W m-2(window)]
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZEMIT_LW_WL_A ! Longwave radiation emitted by the wall [W m-2(wall)]
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZEMIT_LW_WL_B ! Longwave radiation emitted by the wall [W m-2(wall)]
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZTS_WL_A    ! surface temperature of wall A at previous time-step
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZTS_WL_B    ! surface temperature of wall B at previous time-step
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZTI_WL_A    ! internal temperature of wall A at previous time-step
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZTI_WL_B    ! internal temperature of wall B at previous time-step
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZTS_WL      ! surface temperature of averaged wall at new time-step
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZIMB_WL     ! wall residual energy imbalance for verification [W m-2]
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZRAD_WL_FL! rad. fluxes from wall to floor [W m-2(wall)]
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZRAD_WL_MA ! rad. fluxes from wall to mass [W m-2(wall)]
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZCONV_WL_BLD ! conv. fluxes from wall to bld [W m-2(wall)]
REAL, DIMENSION(SIZE(B%XTI_BLD)) :: ZLW_W_TO_WIN   ! Radiative heat trasfer coeff wall-window [W K-1 m-2]
INTEGER                        :: IWL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FACADE_E_BUDGET',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
PRAD_WL_FL= XUNDEF
PRAD_WL_MA = XUNDEF
PRAD_WIN_FL = XUNDEF
PRAD_WIN_MA  = XUNDEF
PCONV_WL_BLD = XUNDEF
PCONV_WIN_BLD  = XUNDEF
!
!* surface temperature of the opposite wall
!  ----------------------------------------
!
ZTS_WL_A(:) = T%XT_WALL_A(:,1)
ZTS_WL_B(:) = T%XT_WALL_B(:,1)
IWL = SIZE(T%XT_WALL_A,2)
ZTI_WL_A(:) = T%XT_WALL_A(:,IWL)
ZTI_WL_B(:) = T%XT_WALL_B(:,IWL)
!
! *Convection heat transfer coefficients [W m-2 K-1] from EP Engineering Reference
! --------------------------------------------------------------------------------
!
!
! *opaque balance of the facade for wall A
!  ---------------------------------------
!
 CALL WALL_LAYER_E_BUDGET(TOP, T, B, T%XT_WALL_A, ZTS_WL_B, ZTI_WL_B, PTSTEP, PDN_RD, &
                          PRHOA, PAC_WL, PAC_BLD, PLW_RAD, PPS, PEXNS, DMT%XABS_SW_WALL_A,&
                          PT_CANYON, PTS_RD, PTSN_RD, PTS_GD, PTS_FL, PLW_WA_TO_WB,   &
                          PLW_R_TO_WA, PLW_G_TO_WA, PLW_NR_TO_WA, PLW_WIN_TO_WA,      &
                          PLW_S_TO_WA, PFLX_BLD_WL_A, PDQS_WL_A, DMT%XABS_LW_WALL_A, &
                          ZEMIT_LW_WL_A, DMT%XH_WALL_A, ZIMB_WL, PRADHT_IN, PRAD_RF_WL, &
                          ZRAD_WL_A_WIN, ZRAD_WL_FL, ZRAD_WL_MA, ZCONV_WL_BLD,        &
                          PLOAD_IN_WL)
!
PIMB_WL      = 0.5 * ZIMB_WL
PRAD_WL_FL   = 0.5 * ZRAD_WL_FL
PRAD_WL_MA   = 0.5 * ZRAD_WL_MA 
PCONV_WL_BLD = 0.5 * ZCONV_WL_BLD 
ZRAD_WL_WIN  = 0.5 * ZRAD_WL_A_WIN
!
! *opaque balance of the facade for wall B
!  ---------------------------------------
!
IF (TOP%CWALL_OPT/='UNIF') THEN
 CALL WALL_LAYER_E_BUDGET(TOP, T, B, T%XT_WALL_B, ZTS_WL_A, ZTI_WL_A, PTSTEP, PDN_RD, &
                          PRHOA, PAC_WL, PAC_BLD, PLW_RAD, PPS, PEXNS, DMT%XABS_SW_WALL_B, &
                          PT_CANYON, PTS_RD, PTSN_RD, PTS_GD, PTS_FL, PLW_WA_TO_WB,   &
                          PLW_R_TO_WB, PLW_G_TO_WB, PLW_NR_TO_WB, PLW_WIN_TO_WB,      &
                          PLW_S_TO_WB, PFLX_BLD_WL_B, PDQS_WL_B, DMT%XABS_LW_WALL_B,   &
                          ZEMIT_LW_WL_B, DMT%XH_WALL_B, ZIMB_WL, PRADHT_IN, PRAD_RF_WL,&
                          ZRAD_WL_A_WIN, ZRAD_WL_FL, ZRAD_WL_MA, ZCONV_WL_BLD,        & 
                          PLOAD_IN_WL                                                )
ELSE
  T%XT_WALL_B    = T%XT_WALL_A
  DMT%XH_WALL_B        = DMT%XH_WALL_A
  DMT%XABS_LW_WALL_B   = DMT%XABS_LW_WALL_A
  PDQS_WL_B      = PDQS_WL_A
  PFLX_BLD_WL_B  = PFLX_BLD_WL_A
  ZEMIT_LW_WL_B  = ZEMIT_LW_WL_A
END IF
!
PIMB_WL      = PIMB_WL      + 0.5 * ZIMB_WL
PRAD_WL_FL   = PRAD_WL_FL   + 0.5 * ZRAD_WL_FL
PRAD_WL_MA   = PRAD_WL_MA   + 0.5 * ZRAD_WL_MA 
PCONV_WL_BLD = PCONV_WL_BLD + 0.5 * ZCONV_WL_BLD 
ZRAD_WL_WIN  = ZRAD_WL_WIN  + 0.5 * ZRAD_WL_A_WIN
!
!-------------------------------------------------------------------------------
!
! *Energy Balance for windows (averaged on both walls)
!  ---------------------------------------------------
!
IF (TOP%CBEM == 'DEF') THEN
  ZEMIT_LW_WIN(:) = 0.
END IF
!
IF (TOP%CBEM == 'BEM') THEN
   ZEMIS_WIN(:) = 0.84 ! From Engineering Reference documentation p219 
                       ! (see http://apps1.eere.energy.gov/buildings/energyplus/
   ZLW_W_TO_WIN(:) = PLW_WA_TO_WIN(:) + PLW_WB_TO_WIN(:)
   ZTS_WL(:) = XUNDEF
   WHERE (ZLW_W_TO_WIN(:)>0.) &
   ZTS_WL(:) = ( PLW_WA_TO_WIN(:)*T%XT_WALL_A(:,1)+PLW_WB_TO_WIN(:)*T%XT_WALL_B(:,1) ) / ZLW_W_TO_WIN(:)
   CALL WINDOW_E_BUDGET(B, ZEMIS_WIN, ZLW_W_TO_WIN, PLW_R_TO_WIN, PLW_G_TO_WIN, PLW_NR_TO_WIN, &
                        PLW_S_TO_WIN, PRAD_RF_WIN, ZRAD_WL_WIN, DMT%XABS_SW_WIN, PLW_RAD, PAC_WIN, &
                        PRADHT_IN, PTS_FL, PRHOA, PDN_RD,  PT_CANYON, ZTS_WL, PTS_RD, PTSN_RD, &
                        PTS_GD, PRAD_WIN_FL, PRAD_WIN_MA, PCONV_WIN_BLD, ZEMIT_LW_WIN,          &
                        DMT%XABS_LW_WIN, PLOAD_IN_WIN                                   )
ENDIF

!*        wall, and win emitted LW radiation on snow-free surfaces
!         ----------------------------------------------------------------
!
!
PEMIT_LW_FAC (:) = 0.5*(ZEMIT_LW_WL_A(:)+ZEMIT_LW_WL_B(:)) * (1-B%XGR(:)) +  ZEMIT_LW_WIN(:) * B%XGR(:)
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FACADE_E_BUDGET',1,ZHOOK_HANDLE)
END SUBROUTINE FACADE_E_BUDGET
