!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
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
 CHARACTER(LEN=4), INTENT(IN)      :: HWALL_OPT     ! TEB option for walls
                                                   ! 'UNIF' : uniform walls
                                                   ! 'TWO ' : two different opposite walls
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
!                                                 ! between wall and canyon
REAL, DIMENSION(:), INTENT(IN)    :: PAC_BLD      ! aerodynamical conductance
                                                  ! inside the building itself
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
                                                  ! wall = facade - glazing
REAL, DIMENSION(:), INTENT(OUT)   :: PH_WALL_B    ! Sensible heat flux from wall to air [W/m2(wall)]
                                                  ! wall = facade - glazing
REAL, DIMENSION(:), INTENT(OUT)   :: PIMB_WALL    ! wall residual energy imbalance 
                                                  ! for verification [W m-2]
REAL, DIMENSION(:),   INTENT(IN)  :: PF_WALL_FLOOR ! View factor wall-floor
REAL, DIMENSION(:),   INTENT(IN)  :: PF_WALL_MASS  ! View factor wall-mass
REAL, DIMENSION(:),   INTENT(IN)  :: PF_WALL_WIN   ! View factor wall-win
REAL, DIMENSION(:),   INTENT(IN)  :: PRADHT_IN     ! Indoor radiant heat transfer coefficient
                                                    ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_ROOF_WALL ! rad. fluxes from roof to wall [W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_ROOF_WIN  ! rad. fluxes from roof to win [W m-2(roof)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WALL_FLOOR! rad. fluxes from wall to floor [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WALL_MASS ! rad. fluxes from wall to mass [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WIN_FLOOR ! rad. fluxes from window to floor [W m-2(win)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WIN_MASS  ! rad. fluxes from window to mass [W m-2(win)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_WALL_BLD ! conv. fluxes from wall to bld [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_WIN_BLD  ! conv. fluxes from window to bld [W m-2(wind.)]
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WIN        ! window aerodynamic conductance
!
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_IN_WALL  ! solar + inter. heat gains W/m2 [wall]
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_IN_WIN   ! solar + inter. heat gains W/m2 [win]
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PTI_BLD)) :: ZRAD_WALL_WIN   ! rad. fluxes between averaged wall and win [W m-2(wall)]
REAL, DIMENSION(SIZE(PTI_BLD)) :: ZRAD_WALL_A_WIN ! rad. fluxes between one wall and win [W m-2(wall)]
REAL, DIMENSION(SIZE(PTI_BLD)) :: ZEMIS_WIN     ! window emissivity
REAL, DIMENSION(SIZE(PTI_BLD)) :: ZEMIT_LW_WIN  ! Longwave radiation emitted by the window [W m-2(window)]
REAL, DIMENSION(SIZE(PTI_BLD)) :: ZEMIT_LW_WALL_A ! Longwave radiation emitted by the wall [W m-2(wall)]
REAL, DIMENSION(SIZE(PTI_BLD)) :: ZEMIT_LW_WALL_B ! Longwave radiation emitted by the wall [W m-2(wall)]
REAL, DIMENSION(SIZE(PTI_BLD)) :: ZTS_WALL_A    ! surface temperature of wall A at previous time-step
REAL, DIMENSION(SIZE(PTI_BLD)) :: ZTS_WALL_B    ! surface temperature of wall B at previous time-step
REAL, DIMENSION(SIZE(PTI_BLD)) :: ZTI_WALL_A    ! internal temperature of wall A at previous time-step
REAL, DIMENSION(SIZE(PTI_BLD)) :: ZTI_WALL_B    ! internal temperature of wall B at previous time-step
REAL, DIMENSION(SIZE(PTI_BLD)) :: ZTS_WALL      ! surface temperature of averaged wall at new time-step
REAL, DIMENSION(SIZE(PTI_BLD)) :: ZIMB_WALL     ! wall residual energy imbalance for verification [W m-2]
REAL, DIMENSION(SIZE(PTI_BLD)) :: ZRAD_WALL_FLOOR! rad. fluxes from wall to floor [W m-2(wall)]
REAL, DIMENSION(SIZE(PTI_BLD)) :: ZRAD_WALL_MASS ! rad. fluxes from wall to mass [W m-2(wall)]
REAL, DIMENSION(SIZE(PTI_BLD)) :: ZCONV_WALL_BLD ! conv. fluxes from wall to bld [W m-2(wall)]
REAL, DIMENSION(SIZE(PTI_BLD)) :: ZLW_W_TO_WIN   ! Radiative heat trasfer coeff wall-window [W K-1 m-2]
INTEGER                        :: IWALL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FACADE_E_BUDGET',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
PRAD_WALL_FLOOR= XUNDEF
PRAD_WALL_MASS = XUNDEF
PRAD_WIN_FLOOR = XUNDEF
PRAD_WIN_MASS  = XUNDEF
PCONV_WALL_BLD = XUNDEF
PCONV_WIN_BLD  = XUNDEF
!
!* surface temperature of the opposite wall
!  ----------------------------------------
!
ZTS_WALL_A(:) = PT_WALL_A(:,1)
ZTS_WALL_B(:) = PT_WALL_B(:,1)
IWALL = SIZE(PT_WALL_A,2)
ZTI_WALL_A(:) = PT_WALL_A(:,IWALL)
ZTI_WALL_B(:) = PT_WALL_B(:,IWALL)
!
! *Convection heat transfer coefficients [W m-2 K-1] from EP Engineering Reference
! --------------------------------------------------------------------------------
!
!
! *opaque balance of the facade for wall A
!  ---------------------------------------
!
 CALL WALL_LAYER_E_BUDGET(HBEM, PT_WALL_A, ZTS_WALL_B, ZTI_WALL_B, PTSTEP,                     &
                         PHC_WALL, PTC_WALL, PD_WALL, PDN_ROAD, PRHOA, PAC_WALL,  &
                         PAC_BLD, PTI_BLD, PLW_RAD, PPS, PEXNS, PABS_SW_WALL_A, PT_CANYON,    &
                         PTS_ROAD, PTSNOW_ROAD, PTS_GARDEN, PTS_MASS, PTS_FLOOR, PEMIS_WALL,  &
                         PLW_WA_TO_WB, PLW_R_TO_WA, PLW_G_TO_WA, PLW_NR_TO_WA,                &
                         PLW_WIN_TO_WA, PLW_S_TO_WA,                                          &
                         PT_WIN1, PFLX_BLD_WALL_A, PDQS_WALL_A, PT_WIN2, PABS_LW_WALL_A,      &
                         ZEMIT_LW_WALL_A, PH_WALL_A, ZIMB_WALL, PF_WALL_MASS,                 &
                         PF_WALL_FLOOR, PF_WALL_WIN, PRADHT_IN, PRAD_ROOF_WALL,               &
                         ZRAD_WALL_A_WIN, ZRAD_WALL_FLOOR, ZRAD_WALL_MASS, ZCONV_WALL_BLD,    &
                         PLOAD_IN_WALL)
!
PIMB_WALL       = 0.5 * ZIMB_WALL
PRAD_WALL_FLOOR = 0.5 * ZRAD_WALL_FLOOR
PRAD_WALL_MASS  = 0.5 * ZRAD_WALL_MASS 
PCONV_WALL_BLD  = 0.5 * ZCONV_WALL_BLD 
ZRAD_WALL_WIN   = 0.5 * ZRAD_WALL_A_WIN
!
! *opaque balance of the facade for wall B
!  ---------------------------------------
!
IF (HWALL_OPT/='UNIF') THEN
 CALL WALL_LAYER_E_BUDGET(HBEM, PT_WALL_B, ZTS_WALL_A, ZTI_WALL_A, PTSTEP,                     &
                         PHC_WALL, PTC_WALL, PD_WALL, PDN_ROAD, PRHOA, PAC_WALL,              &
                         PAC_BLD, PTI_BLD, PLW_RAD, PPS, PEXNS, PABS_SW_WALL_B, PT_CANYON,    &
                         PTS_ROAD, PTSNOW_ROAD, PTS_GARDEN, PTS_MASS, PTS_FLOOR, PEMIS_WALL,  &
                         PLW_WA_TO_WB, PLW_R_TO_WB, PLW_G_TO_WB, PLW_NR_TO_WB,                &
                         PLW_WIN_TO_WB, PLW_S_TO_WB,                                          &
                         PT_WIN1, PFLX_BLD_WALL_B, PDQS_WALL_B, PT_WIN2, PABS_LW_WALL_B,      &
                         ZEMIT_LW_WALL_B, PH_WALL_B, ZIMB_WALL, PF_WALL_MASS,                 &
                         PF_WALL_FLOOR, PF_WALL_WIN, PRADHT_IN, PRAD_ROOF_WALL,               &
                         ZRAD_WALL_A_WIN, ZRAD_WALL_FLOOR, ZRAD_WALL_MASS, ZCONV_WALL_BLD,    & 
                         PLOAD_IN_WALL                                                        )
ELSE
  PT_WALL_B        = PT_WALL_A
  PH_WALL_B        = PH_WALL_A
  PABS_LW_WALL_B   = PABS_LW_WALL_A
  PDQS_WALL_B      = PDQS_WALL_A
  PFLX_BLD_WALL_B  = PFLX_BLD_WALL_A
  ZEMIT_LW_WALL_B  = ZEMIT_LW_WALL_A
END IF
!
PIMB_WALL       = PIMB_WALL       + 0.5 * ZIMB_WALL
PRAD_WALL_FLOOR = PRAD_WALL_FLOOR + 0.5 * ZRAD_WALL_FLOOR
PRAD_WALL_MASS  = PRAD_WALL_MASS  + 0.5 * ZRAD_WALL_MASS 
PCONV_WALL_BLD  = PCONV_WALL_BLD  + 0.5 * ZCONV_WALL_BLD 
ZRAD_WALL_WIN   = ZRAD_WALL_WIN   + 0.5 * ZRAD_WALL_A_WIN
!
!-------------------------------------------------------------------------------
!
! *Energy Balance for windows (averaged on both walls)
!  ---------------------------------------------------
!
IF (HBEM == 'DEF') THEN
  ZEMIT_LW_WIN(:) = 0.
  PABS_LW_WIN (:) = XUNDEF
END IF
!
IF (HBEM == 'BEM') THEN
        ZEMIS_WIN(:) = 0.84 ! From Engineering Reference documentation p219 
                            ! (see http://apps1.eere.energy.gov/buildings/energyplus/
   ZLW_W_TO_WIN(:) = PLW_WA_TO_WIN(:) + PLW_WB_TO_WIN(:)
   ZTS_WALL(:) = XUNDEF
   WHERE (ZLW_W_TO_WIN(:)>0.) &
   ZTS_WALL(:) = ( PLW_WA_TO_WIN(:)*PT_WALL_A(:,1)+PLW_WB_TO_WIN(:)*PT_WALL_B(:,1) ) / ZLW_W_TO_WIN(:)
   CALL WINDOW_E_BUDGET(PU_WIN, ZEMIS_WIN, PF_WIN_FLOOR, PF_WIN_WALL, PF_WIN_MASS,                &
                        ZLW_W_TO_WIN, PLW_R_TO_WIN, PLW_G_TO_WIN, PLW_NR_TO_WIN, PLW_S_TO_WIN,    &
                        PRAD_ROOF_WIN, ZRAD_WALL_WIN, PABS_SW_WIN, PLW_RAD, PAC_WIN,              &
                        PRADHT_IN, PTI_BLD, PTS_MASS, PTS_FLOOR, PRHOA, PDN_ROAD,                 &
                        PT_CANYON, ZTS_WALL, PTS_ROAD, PTSNOW_ROAD, PTS_GARDEN,                   &
                        PT_WIN1, PT_WIN2, PRAD_WIN_FLOOR, PRAD_WIN_MASS, PCONV_WIN_BLD,           &
                        ZEMIT_LW_WIN, PABS_LW_WIN, PLOAD_IN_WIN                                   )
ENDIF

!*        wall, and win emitted LW radiation on snow-free surfaces
!         ----------------------------------------------------------------
!
!
PEMIT_LW_FAC (:) = 0.5*(ZEMIT_LW_WALL_A(:)+ZEMIT_LW_WALL_B(:)) * (1-PGR(:)) +  ZEMIT_LW_WIN(:) * PGR(:)
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FACADE_E_BUDGET',1,ZHOOK_HANDLE)
END SUBROUTINE FACADE_E_BUDGET
