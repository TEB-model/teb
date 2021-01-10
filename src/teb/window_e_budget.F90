!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
SUBROUTINE WINDOW_E_BUDGET(B, PEMIS_WIN, PLW_W_TO_WIN, PLW_R_TO_WIN, PLW_G_TO_WIN, &
                           PLW_NR_TO_WIN, PLW_S_TO_WIN, PRAD_RF_WIN, PRAD_WL_WIN,  &
                           PABS_SW_WIN, PLW_RAD, PAC_WL, PRADHT_IN, PTS_FL, PRHOA, &
                           PDN_RD, PT_CANYON, PTS_WL, PTS_RD, PTSN_RD, PTS_GD,     &
                           PRAD_WIN_FL, PRAD_WIN_MA, PCONV_WIN_BLD, PEMIT_LW_WIN,  &
                           PABS_LW_WIN, PLOAD_IN_WIN)
!###################################################################################################################################
!
!!****  *WINDOW_E_BUDGET*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the evoultion of window temperature
!         
!     
!!**  METHOD
!     ------
!
! window is supposed double pane with no thermal capacity
! a steady state energy balance is applied to both layers that exchanges U_WIN(T1 - T2)
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
!!      B. Bueno           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2011 
!!      G. Pigeon   10/11/2011 exteranlized window balance from wall_layer_e_budget to window_e_budget and use fluxes from layers for
!                              which the balance has been computed before (roof and wall) and returns the fluxes to the compartments
!                              computed afterwards (floor, mass, bld)
!!      G. Pigeon      09/2012 new indoor conv. coef
!!      G. Pigeon      10/2012 separate abs from outdoor on both side of window
!!                            + add loads for indoor face of the window
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_BEM_n, ONLY : BEM_t
!
USE MODD_CSTS,ONLY : XCPD, XSTEFAN
USE MODE_CONV_DOE
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(BEM_t), INTENT(INOUT) :: B
!
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_WIN    ! window emissivity
REAL, DIMENSION(:), INTENT(IN)    :: PLW_W_TO_WIN ! Radiative heat trasfer coeff window-wall outdoor
                                                  ! [W K-1 m-2] 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_WIN ! Radiative heat trasfer coeff window-road 
                                                  ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_WIN ! Radiative heat trasfer coeff window-garden
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_WIN! Radiative heat trasfer coeff window-snow
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_WIN ! Radiative heat trasfer coeff window-sky 
                                                  ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_RF_WIN  ! rad. fluxes from roof to win [W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_WL_WIN  ! rad. fluxes from wall to win [W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WIN  ! window absorbed shortwave radiation [W m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD      ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WL     ! aerodynamical conductance between wall and canyon
REAL, DIMENSION(:),   INTENT(IN)  :: PRADHT_IN     ! Indoor radiant heat transfer coefficient
                                                    ! [W K-1 m-2]
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
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PT_CANYON)) :: ZT_SKY         ! sky temperature [K]
REAL, DIMENSION(SIZE(PT_CANYON)) :: ZA12           ! 
REAL, DIMENSION(SIZE(PT_CANYON)) :: ZB2            ! auxiliar variables for window temperature
REAL, DIMENSION(SIZE(PT_CANYON)) :: ZB1            ! calculations
REAL, DIMENSION(SIZE(PT_CANYON)) :: ZA21           !
REAL, DIMENSION(SIZE(PT_CANYON)) :: ZCHTC_IN_WIN   ! indoor convective heat transfer coef. window [W m-2 K-1]
INTEGER :: JJ

!*      Preliminaries
!       -------------
ZCHTC_IN_WIN = CHTC_VERT_DOE(B%XT_WIN2, B%XTI_BLD)
DO JJ=1,SIZE(ZCHTC_IN_WIN)
   ZCHTC_IN_WIN(JJ) = MAX(1., ZCHTC_IN_WIN(JJ))
ENDDO

!*      1.    Sky temperature
!             ---------------
ZT_SKY(:) = (PLW_RAD(:)/XSTEFAN)**0.25
!  
ZA12(:) = ZCHTC_IN_WIN(:)   &
         + PRADHT_IN(:) * B%XF_WIN_MASS (:)  &
         + PRADHT_IN(:) * B%XF_WIN_FLOOR(:) &
         + B%XUGG_WIN(:)
!
ZB1(:)  = ZCHTC_IN_WIN   (:) * B%XTI_BLD (:)             &
         + PRADHT_IN     (:) * B%XF_WIN_MASS(:) * B%XT_MASS(:,1)             &
         + B%XF_WIN_WALL (:) * PRAD_WL_WIN(:) &
         + B%XF_WIN_FLOOR(:) * PRAD_RF_WIN(:) &
         + PRADHT_IN     (:) * B%XF_WIN_FLOOR(:)* PTS_FL(:) &
         + PLOAD_IN_WIN(:) + PABS_SW_WIN(:) /2.
!
ZA21(:) = - PAC_WL       (:) * PRHOA(:)*XCPD   &
          - PLW_S_TO_WIN (:)                  &
          - PLW_W_TO_WIN (:)                  &
          - PLW_R_TO_WIN (:) *(1.-PDN_RD(:))&
          - PLW_G_TO_WIN (:)                  &
          - PLW_NR_TO_WIN(:) * PDN_RD(:)    &
          - B%XUGG_WIN(:)
!
ZB2(:)  = - PAC_WL       (:) * PRHOA(:)*XCPD*PT_CANYON(:)   &
          - PLW_S_TO_WIN (:) * ZT_SKY(:)                     &
          - PLW_W_TO_WIN (:) * PTS_WL(:)                   &
          - PLW_R_TO_WIN (:) *(1.-PDN_RD(:)) * PTS_RD(:) &
          - PLW_NR_TO_WIN(:) * PDN_RD(:) * PTSN_RD(:) &
          - PABS_SW_WIN  (:) /2.
!
IF (SIZE(PTS_GD)>0) THEN
  ZB2(:) = ZB2(:) - PLW_G_TO_WIN (:) * PTS_GD(:)
ENDIF
!
! compute outdoor temperature
B%XT_WIN1(:) = ( ZB2(:) - ZB1(:)*B%XUGG_WIN(:)/ZA12(:) ) / &
                  ( ZA21(:) + B%XUGG_WIN(:)**2/ZA12(:) )
!
! compute indoor temperature
B%XT_WIN2(:) = (ZB1(:) + B%XUGG_WIN(:)*B%XT_WIN1(:) ) / ZA12(:)
!
! outdoor infrared radiation absorded by the window
PABS_LW_WIN(:)    = PLW_S_TO_WIN (:) * (ZT_SKY (:) - B%XT_WIN1(:)) + &
  (1.-PDN_RD(:)) *  PLW_R_TO_WIN (:) * (PTS_RD (:) - B%XT_WIN1(:)) + &
                    PLW_W_TO_WIN (:) * (PTS_WL (:) - B%XT_WIN1(:)) + &
        PDN_RD(:) * PLW_NR_TO_WIN(:) * (PTSN_RD(:) - B%XT_WIN1(:))
!
IF (SIZE(PTS_GD)>0) THEN
  PABS_LW_WIN(:) = PABS_LW_WIN(:) +  PLW_G_TO_WIN (:) * (PTS_GD (:) - B%XT_WIN1(:))      
ENDIF
!
! outdoor infrared radiation emited by the window
PEMIT_LW_WIN(:) = XSTEFAN * B%XT_WIN1(:)**4 + (1 - PEMIS_WIN(:))/PEMIS_WIN(:) * PABS_LW_WIN(:)
!
PRAD_WIN_FL  (:)  = PRADHT_IN(:)    * (B%XT_WIN2(:) - PTS_FL(:))
PRAD_WIN_MA  (:)  = PRADHT_IN(:)    * (B%XT_WIN2(:) - B%XT_MASS(:,1))
PCONV_WIN_BLD(:)  = ZCHTC_IN_WIN(:) * (B%XT_WIN2(:) - B%XTI_BLD(:))
!
END SUBROUTINE WINDOW_E_BUDGET
