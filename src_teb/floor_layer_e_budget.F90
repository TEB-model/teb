!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ##########################################################################
    SUBROUTINE FLOOR_LAYER_E_BUDGET(PT_FLOOR, PTSTEP, PHC_FLOOR, PTC_FLOOR, PD_FLOOR, &
                                    PFLX_BLD_FLOOR, PDQS_FLOOR, PIMB_FLOOR,           &
                                    PF_FLOOR_MASS, PF_FLOOR_WALL, PF_FLOOR_WIN,       &
                                    PF_FLOOR_ROOF, PRADHT_IN,           &
                                    PTS_MASS, PRAD_WALL_FLOOR, PRAD_ROOF_FLOOR,       &
                                    PRAD_WIN_FLOOR, PLOAD_FLOOR, PTI_BLD,             &
                                    PRAD_FLOOR_MASS, PCONV_FLOOR_BLD                  )
!   ##########################################################################
!
!!****  *FLOOR_LAYER_E_BUDGET*  
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
!!      Original    15/04/09 
!!      G. Pigeon      08/10 computation of residual of energy balance 
!!                           modification of the limit condition for
!!                           the deep temp. from the deep road temp.
!!                           to zero flux condition. idem for sfce T
!!      G. Pigeon      11/11 split in 2 for floor and mass to be flux conservative
!!      G. Pigeon      09/12 modif of conv. coef + implicitation of the exchange
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
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_FLOOR     ! floor layers temperatures
REAL,                 INTENT(IN)    :: PTSTEP       ! time step
REAL, DIMENSION(:,:), INTENT(IN)    :: PHC_FLOOR    ! heat capacity for road layers
REAL, DIMENSION(:,:), INTENT(IN)    :: PTC_FLOOR    ! thermal conductivity for 
                                                    !road layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_FLOOR       ! depth of road layers
REAL, DIMENSION(:),   INTENT(OUT)  :: PFLX_BLD_FLOOR !flux from building to floor
REAL, DIMENSION(:),   INTENT(OUT) :: PDQS_FLOOR !heat storage inside the floor
REAL, DIMENSION(:),   INTENT(OUT) :: PIMB_FLOOR !floor energy residual imbalance for verification
REAL, DIMENSION(:), INTENT(IN)    :: PF_FLOOR_MASS  ! View factor floor-mass
REAL, DIMENSION(:), INTENT(IN)    :: PF_FLOOR_WALL  ! View factor floor-wall
REAL, DIMENSION(:), INTENT(IN)    :: PF_FLOOR_WIN   ! View factor floor-window
REAL, DIMENSION(:), INTENT(IN)    :: PF_FLOOR_ROOF  ! View factor floor-roof
REAL, DIMENSION(:), INTENT(IN)    :: PRADHT_IN      ! Indoor radiant heat transfer coefficient
                                                    ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PTS_MASS  ! surf. mass temp. (contact with bld air)
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_ROOF_FLOOR ! rad. fluxes from roof to floor[W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_WALL_FLOOR ! rad. fluxes from wall to floor[W m-2(wall)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_WIN_FLOOR  ! rad. fluxes from win to floor[W m-2(win)]
REAL, DIMENSION(:), INTENT(IN)    :: PTI_BLD   ! indoor air temp.
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_FLOOR ! solar and internal load to the floor
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_FLOOR_MASS  ! rad. fluxes from floor to mass [W m-2(floor)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_FLOOR_BLD  ! conv. fluxes from floor to bld [W m-2(floor)]
!
!*      0.2    declarations of local variables
!
!
REAL :: ZIMPL=1.0      ! implicit coefficient
REAL :: ZEXPL=0.0      ! explicit coefficient
!
REAL, DIMENSION(SIZE(PT_FLOOR,1),SIZE(PT_FLOOR,2)) :: ZA,& ! lower diag.
                                                      ZB,& ! main  diag.
                                                      ZC,& ! upper diag.
                                                      ZY   ! r.h.s.
!
REAL, DIMENSION(SIZE(PT_FLOOR,1)) :: ZTS_FLOOR ! surf. floor temp.  used for rad. exchanges
REAL, DIMENSION(SIZE(PT_FLOOR,1)) :: ZTS_FLOOR_CONV ! surf. floor temp. used for conv exchanges
                                               ! used during calculation
REAL, DIMENSION(SIZE(PT_FLOOR,1)) :: ZCHTC_IN_FLOOR   ! Indoor floor convec heat transfer coefficient
                                                    ! [W K-1 m-2(bld)]

REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER :: JJ
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FLOOR_LAYER_E_BUDGET',0,ZHOOK_HANDLE)
!
! *Convection heat transfer coefficients [W m-2 K-1]
!  From EP Engineering Reference
ZCHTC_IN_FLOOR(:) = CHTC_UP_DOE(PT_FLOOR(:,1), PTI_BLD(:))
DO JJ=1,SIZE(ZCHTC_IN_FLOOR)
   ZCHTC_IN_FLOOR(JJ) = MAX(1., ZCHTC_IN_FLOOR(JJ))
ENDDO


!
 CALL LAYER_E_BUDGET_GET_COEF( PT_FLOOR, PTSTEP, ZIMPL, PHC_FLOOR, PTC_FLOOR, PD_FLOOR, &
                              ZA, ZB, ZC, ZY )
!
ZTS_FLOOR(:) = PT_FLOOR(:,1) 

ZB(:,1) = ZB(:,1) + ZIMPL * &
   (ZCHTC_IN_FLOOR(:)*4./3. + PRADHT_IN(:) * PF_FLOOR_MASS(:))


ZY(:,1) = ZY(:,1)  &
   + ZCHTC_IN_FLOOR(:) * (PTI_BLD(:) - 1./3. * PT_FLOOR(:, 1) * (4* ZEXPL -1))  &
   + PF_FLOOR_WIN  (:) * PRAD_WIN_FLOOR(:)                    &
   + PF_FLOOR_WALL (:) * PRAD_WALL_FLOOR(:)                   &
   + PF_FLOOR_ROOF (:) * PRAD_ROOF_FLOOR(:)                   &
   + PRADHT_IN(:) * PF_FLOOR_MASS (:) * (PTS_MASS(:) - ZEXPL * PT_FLOOR(:,1))  &
   + PLOAD_FLOOR(:)
!
 CALL LAYER_E_BUDGET( PT_FLOOR, PTSTEP, ZIMPL, PHC_FLOOR, PTC_FLOOR, PD_FLOOR, &
                     ZA, ZB, ZC, ZY, PDQS_FLOOR )
!
!*      floor surface temperature used in the implicit formulation
!       ----------------------------------------------------------
ZTS_FLOOR_CONV(:) = 4./3. * ZIMPL * PT_FLOOR(:,1) + 1./3. * ZTS_FLOOR(:) * (4 * ZEXPL - 1.)
ZTS_FLOOR(:) = ZEXPL * ZTS_FLOOR(:) + ZIMPL * PT_FLOOR(:,1)
!
!*      fluxes with mass and indoor air
!       ----------------------------------------------------------
PRAD_FLOOR_MASS(:) = PRADHT_IN(:)    * (ZTS_FLOOR(:) - PTS_MASS(:))
PCONV_FLOOR_BLD(:) = ZCHTC_IN_FLOOR(:) * (ZTS_FLOOR_CONV(:) - PTI_BLD (:))
!
!*     Flux between floor and indoor surfaces and air
!       ------------------------------------------------
PFLX_BLD_FLOOR(:) = - PCONV_FLOOR_BLD(:)  &
       + PF_FLOOR_WIN  (:) * PRAD_WIN_FLOOR(:)  &
       + PF_FLOOR_WALL (:) * PRAD_WALL_FLOOR(:) &
       + PF_FLOOR_ROOF (:) * PRAD_ROOF_FLOOR(:) &
       + PRADHT_IN(:) * PF_FLOOR_MASS(:) * (PTS_MASS (:)  - ZTS_FLOOR(:)) &
       + PLOAD_FLOOR(:)
!
!*     Floor residual energy imbalance for verification
!       ------------------------------------------------
PIMB_FLOOR(:) = PFLX_BLD_FLOOR(:) - PDQS_FLOOR(:)
!
IF (LHOOK) CALL DR_HOOK('FLOOR_LAYER_E_BUDGET',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE FLOOR_LAYER_E_BUDGET
