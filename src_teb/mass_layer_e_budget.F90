!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ##########################################################################
    SUBROUTINE MASS_LAYER_E_BUDGET(PT_MASS, PTSTEP, PHC_MASS, PTC_MASS, PD_MASS, &
                                    PFLX_BLD_MASS, PDQS_MASS, PIMB_MASS,           &
                                    PF_MASS_WALL, PF_MASS_WIN,       &
                                    PF_MASS_FLOOR, PRADHT_IN,           &
                                    PRAD_WALL_MASS, PRAD_ROOF_MASS,       &
                                    PRAD_WIN_MASS, PLOAD_MASS, PTI_BLD,             &
                                    PRAD_FLOOR_MASS, PCONV_MASS_BLD                  )
!   ##########################################################################
!
!!****  *MASS_LAYER_E_BUDGET*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the evoultion of building floor temperatures
!         
!     
!!**  METHOD
!     ------
!
!    6 : equations for evolution of Ts_floor 
!        *************************************************************
!
!     dTf_k(t) / dt = 1/(df_k*Cf_k) * (- 2*Kf_k-1*(Tf_k-Tf_k-1)/(df_k-1 +df_k) 
!                                      - 2*Kf_k  *(Tf_k-Tf_k+1)/(df_k+1 +df_k) )
!
!     dTf_1(t) / dt = 1/(df_1*Cf_1) * (- 2*Kw_1*(Tw_1-Tw_2)/(dw_1 +dw_2))
!
!       with
!
!   K*_k  = (d*_k+ d*_k+1)/(d*_k/k*_k+ d*_k+1/k*_k+1)
!
!
! The system is implicited (or semi-implicited).
!
! ZIMPL=1    ---> implicit system
! ZIMPL=0.5  ---> semi-implicit system
! ZIMPL=0    ---> explicit system
!
!
!
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
!!	G. Pigeon           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    11/11
!!      G. Pigeon   09/12 modif internal convection coef
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODI_LAYER_E_BUDGET_GET_COEF
USE MODI_LAYER_E_BUDGET
USE MODE_CONV_DOE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_MASS     ! floor layers temperatures
REAL,                 INTENT(IN)    :: PTSTEP       ! time step
REAL, DIMENSION(:,:), INTENT(IN)    :: PHC_MASS    ! heat capacity for road layers
REAL, DIMENSION(:,:), INTENT(IN)    :: PTC_MASS    ! thermal conductivity for 
                                                    !road layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_MASS       ! depth of road layers
REAL, DIMENSION(:),   INTENT(OUT)  :: PFLX_BLD_MASS !flux from building to floor
REAL, DIMENSION(:),   INTENT(OUT) :: PDQS_MASS !heat storage inside the floor
REAL, DIMENSION(:),   INTENT(OUT) :: PIMB_MASS !floor energy residual imbalance for verification
REAL, DIMENSION(:), INTENT(IN)    :: PF_MASS_WALL  ! View factor mass-wall
REAL, DIMENSION(:), INTENT(IN)    :: PF_MASS_WIN   ! View factor mass-window
REAL, DIMENSION(:), INTENT(IN)    :: PF_MASS_FLOOR ! View factor mass-floor
REAL, DIMENSION(:), INTENT(IN)    :: PRADHT_IN      ! Indoor radiant heat transfer coefficient
                                                    ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_ROOF_MASS ! rad. fluxes from roof to floor[W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_WALL_MASS ! rad. fluxes from wall to floor[W m-2(wall)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_WIN_MASS  ! rad. fluxes from win to floor[W m-2(win)]
REAL, DIMENSION(:), INTENT(IN)    :: PTI_BLD   ! indoor air temp.
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_MASS ! solar and internal load to the floor
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_FLOOR_MASS  ! rad. fluxes from floor to mass [W m-2(floor)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_MASS_BLD  ! conv. fluxes from floor to bld [W m-2(floor)]
!
!*      0.2    declarations of local variables
!
!
REAL :: ZIMPL=1.0      ! implicit coefficient
REAL :: ZEXPL=0.0      ! explicit coefficient
!
REAL, DIMENSION(SIZE(PT_MASS,1),SIZE(PT_MASS,2)) :: ZA,& ! lower diag.
                                                      ZB,& ! main  diag.
                                                      ZC,& ! upper diag.
                                                      ZY   ! r.h.s.
!
REAL, DIMENSION(SIZE(PT_MASS,1)) :: ZTS_MASS  ! surf. mass temp.
                                              ! used during calculation
REAL, DIMENSION(SIZE(PT_MASS,1)) :: ZTS_MASS_CONV  ! surf. mass temp. used for conv flux
REAL, DIMENSION(SIZE(PT_MASS,1)) :: ZCHTC_IN_MASS ! Indoor floor convec heat transfer coefficient
                                                ! [W K-1 m-2(bld)]
INTEGER :: JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MASS_LAYER_E_BUDGET',0,ZHOOK_HANDLE)
!
! *Convection heat transfer coefficients [W m-2 K-1] from EP Engineering Reference
!
ZCHTC_IN_MASS(:) = CHTC_VERT_DOE(PT_MASS(:,1), PTI_BLD(:)) 
DO JJ=1,SIZE(ZCHTC_IN_MASS)
   ZCHTC_IN_MASS(JJ) = MAX(1., ZCHTC_IN_MASS(JJ))
ENDDO
!
 CALL LAYER_E_BUDGET_GET_COEF( PT_MASS, PTSTEP, ZIMPL, PHC_MASS, PTC_MASS, PD_MASS, &
                              ZA, ZB, ZC, ZY )
!
ZTS_MASS(:) = PT_MASS(:,1) 

ZB(:,1) = ZB(:,1) + ZIMPL * 4./3. * ZCHTC_IN_MASS(:)

ZY(:,1) = ZY(:,1)  &
   + ZCHTC_IN_MASS(:) * (PTI_BLD(:) - 1./3. * PT_MASS(:, 1) * (4 * ZEXPL -1.))  &
   + PF_MASS_WIN  (:) * PRAD_WIN_MASS(:)                    &
   + PF_MASS_WALL (:) * PRAD_WALL_MASS(:)                   &
   + PF_MASS_FLOOR (:) * (PRAD_ROOF_MASS(:) +PRAD_FLOOR_MASS(:)) &
   + PLOAD_MASS(:)
!
 CALL LAYER_E_BUDGET( PT_MASS, PTSTEP, ZIMPL, PHC_MASS, PTC_MASS, PD_MASS, &
                     ZA, ZB, ZC, ZY, PDQS_MASS )
!
!*      calculation of temperature used in energy balance calculation
!       -------------------------------------------------------------
!
ZTS_MASS_CONV(:) = ZIMPL * 4./3. * PT_MASS(:,1) +1./3 * ZTS_MASS(:) * (4 * ZEXPL -1.)
ZTS_MASS(:) = ZEXPL * ZTS_MASS(:) + ZIMPL * PT_MASS(:,1)
!
!*      calculation of convection flux between mass and building air
!       ------------------------------------------------------------
!
PCONV_MASS_BLD(:) = ZCHTC_IN_MASS(:) * (ZTS_MASS_CONV(:) - PTI_BLD(:))
!
!*      For diagnostics calculation of flux exchanged between the mass and the
!       indoor
!       ------------------------------------------------
!
PFLX_BLD_MASS(:) = - PCONV_MASS_BLD(:) &
       + PF_MASS_WIN  (:) * PRAD_WIN_MASS(:)                    &
       + PF_MASS_WALL (:) * PRAD_WALL_MASS(:)                   &
       + PF_MASS_FLOOR (:) * (PRAD_ROOF_MASS(:) + PRAD_FLOOR_MASS(:)) &
       + PLOAD_MASS(:)
!
!*      Floor residual energy imbalance for verification
!       ------------------------------------------------
!
PIMB_MASS(:) = PFLX_BLD_MASS(:) - PDQS_MASS(:)
!
IF (LHOOK) CALL DR_HOOK('MASS_LAYER_E_BUDGET',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE MASS_LAYER_E_BUDGET
