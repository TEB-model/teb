!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
SUBROUTINE WINDOW_E_BUDGET(PU_WIN, PEMIS_WIN, PF_WIN_FLOOR, PF_WIN_WALL, PF_WIN_MASS, &
                           PLW_W_TO_WIN, PLW_R_TO_WIN, PLW_G_TO_WIN,                  &
                           PLW_NR_TO_WIN, PLW_S_TO_WIN, PRAD_ROOF_WIN, PRAD_WALL_WIN, &
                           PABS_SW_WIN, PLW_RAD, PAC_WALL,                 &
                           PRADHT_IN, PTI_BLD, PTS_MASS, PTS_FLOOR, PRHOA, PDN_ROAD,  &
                           PT_CANYON, PTS_WALL, PTS_ROAD, PTSNOW_ROAD,                &
                           PTS_GARDEN, PT_WIN1, PT_WIN2, PRAD_WIN_FLOOR,              &
                           PRAD_WIN_MASS, PCONV_WIN_BLD, PEMIT_LW_WIN, PABS_LW_WIN, PLOAD_IN_WIN)
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
USE MODD_CSTS,ONLY : XCPD, XSTEFAN
USE MODE_CONV_DOE
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN)    :: PU_WIN       ! window U-factor [W m-2 K-1]
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_WIN    ! window emissivity
REAL, DIMENSION(:), INTENT(IN)    :: PF_WIN_FLOOR ! window floor view factor
REAL, DIMENSION(:), INTENT(IN)    :: PF_WIN_WALL  ! window wall view factor
REAL, DIMENSION(:), INTENT(IN)    :: PF_WIN_MASS  ! window mass view factor
REAL, DIMENSION(:), INTENT(IN)    :: PLW_W_TO_WIN ! Radiative heat trasfer coeff window-wall outdoor
                                                  ! [W K-1 m-2] 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_WIN ! Radiative heat trasfer coeff window-road 
                                                  ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_WIN ! Radiative heat trasfer coeff window-garden
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_WIN! Radiative heat trasfer coeff window-snow
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_WIN ! Radiative heat trasfer coeff window-sky 
                                                  ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_ROOF_WIN  ! rad. fluxes from roof to win [W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_WALL_WIN  ! rad. fluxes from wall to win [W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WIN  ! window absorbed shortwave radiation [W m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD      ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WALL     ! aerodynamical conductance between wall and canyon
REAL, DIMENSION(:),   INTENT(IN)  :: PRADHT_IN     ! Indoor radiant heat transfer coefficient
                                                    ! [W K-1 m-2]
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
ZCHTC_IN_WIN = CHTC_VERT_DOE(PT_WIN2, PTI_BLD)
DO JJ=1,SIZE(ZCHTC_IN_WIN)
   ZCHTC_IN_WIN(JJ) = MAX(1., ZCHTC_IN_WIN(JJ))
ENDDO

!*      1.    Sky temperature
!             ---------------
ZT_SKY(:) = (PLW_RAD(:)/XSTEFAN)**0.25
!  
ZA12(:) = ZCHTC_IN_WIN(:)   &
         + PRADHT_IN(:) * PF_WIN_MASS(:)  &
         + PRADHT_IN(:) * PF_WIN_FLOOR(:) &
         + PU_WIN(:)
!
ZB1(:)  = ZCHTC_IN_WIN   (:) * PTI_BLD (:)             &
         + PRADHT_IN(:) * PF_WIN_MASS(:) * PTS_MASS(:)             &
         + PF_WIN_WALL(:) * PRAD_WALL_WIN(:) &
         + PF_WIN_FLOOR(:) * PRAD_ROOF_WIN(:) &
         + PRADHT_IN(:) * PF_WIN_FLOOR(:)* PTS_FLOOR(:) &
         + PLOAD_IN_WIN(:) + PABS_SW_WIN(:) /2.
!
ZA21(:) = - PAC_WALL     (:) * PRHOA(:)*XCPD   &
           - PLW_S_TO_WIN (:)                  &
           - PLW_W_TO_WIN (:)                  &
           - PLW_R_TO_WIN (:) *(1.-PDN_ROAD(:))&
           - PLW_G_TO_WIN (:)                  &
           - PLW_NR_TO_WIN(:) * PDN_ROAD(:)    &
           - PU_WIN(:)
!
ZB2(:)  = - PAC_WALL     (:) * PRHOA(:)*XCPD*PT_CANYON(:)   &
           - PLW_S_TO_WIN (:) * ZT_SKY(:)                     &
           - PLW_W_TO_WIN (:) * PTS_WALL(:)                   &
           - PLW_R_TO_WIN (:) *(1.-PDN_ROAD(:)) * PTS_ROAD(:) &
           - PLW_G_TO_WIN (:) * PTS_GARDEN(:)                 &
           - PLW_NR_TO_WIN(:) * PDN_ROAD(:) * PTSNOW_ROAD(:) &
           - PABS_SW_WIN  (:) /2.
!
! compute outdoor temperature
PT_WIN1(:) = ( ZB2(:) - ZB1(:)*PU_WIN(:)/ZA12(:) ) / &
              ( ZA21(:) + PU_WIN(:)**2/ZA12(:) )
!
! compute indoor temperature
PT_WIN2(:) = (ZB1(:) + PU_WIN(:)*PT_WIN1(:) ) / ZA12(:)
!
! outdoor infrared radiation absorded by the window
PABS_LW_WIN(:)    = PLW_S_TO_WIN (:) * (ZT_SKY     (:) - PT_WIN1(:)) + &
  (1.-PDN_ROAD(:))*  PLW_R_TO_WIN (:) * (PTS_ROAD   (:) - PT_WIN1(:)) + &
                     PLW_G_TO_WIN (:) * (PTS_GARDEN (:) - PT_WIN1(:)) + &
                     PLW_W_TO_WIN (:) * (PTS_WALL  (:)  - PT_WIN1(:)) + &
      PDN_ROAD(:) * PLW_NR_TO_WIN(:) * (PTSNOW_ROAD(:) - PT_WIN1(:))
!
! outdoor infrared radiation emited by the window
PEMIT_LW_WIN(:) = XSTEFAN * PT_WIN1(:)**4 + (1 - PEMIS_WIN(:))/PEMIS_WIN(:) * PABS_LW_WIN(:)
!
PRAD_WIN_FLOOR(:) = PRADHT_IN(:)   * (PT_WIN2(:) - PTS_FLOOR(:))
PRAD_WIN_MASS(:)  = PRADHT_IN(:)   * (PT_WIN2(:) - PTS_MASS(:))
PCONV_WIN_BLD(:)  = ZCHTC_IN_WIN(:) * (PT_WIN2(:) - PTI_BLD(:))
!
END SUBROUTINE WINDOW_E_BUDGET
