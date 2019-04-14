!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!   ##########################################################################
    SUBROUTINE FLOOR_LAYER_E_BUDGET(B, PTSTEP, PFLX_BLD_FL, PDQS_FL, PIMB_FL, PRADHT_IN,  &
                                    PRAD_WL_FL, PRAD_RF_FL, PRAD_WIN_FL, PLOAD_FL,        &
                                    PRAD_FL_MA, PCONV_FL_BLD                  )
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
!!      G. Pigeon           * Meteo-France *
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
USE MODD_BEM_n, ONLY : BEM_t
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
TYPE(BEM_t), INTENT(INOUT) :: B
!
REAL,                 INTENT(IN)    :: PTSTEP       ! time step
REAL, DIMENSION(:),   INTENT(OUT)  :: PFLX_BLD_FL !flux from building to floor
REAL, DIMENSION(:),   INTENT(OUT) :: PDQS_FL !heat storage inside the floor
REAL, DIMENSION(:),   INTENT(OUT) :: PIMB_FL !floor energy residual imbalance for verification
REAL, DIMENSION(:), INTENT(IN)    :: PRADHT_IN      ! Indoor radiant heat transfer coefficient
                                                    ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_RF_FL ! rad. fluxes from roof to floor[W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_WL_FL ! rad. fluxes from wall to floor[W m-2(wall)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_WIN_FL  ! rad. fluxes from win to floor[W m-2(win)]
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_FL ! solar and internal load to the floor
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_FL_MA  ! rad. fluxes from floor to mass [W m-2(floor)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_FL_BLD  ! conv. fluxes from floor to bld [W m-2(floor)]
!
!*      0.2    declarations of local variables
!
!
REAL :: ZIMPL=1.0      ! implicit coefficient
REAL :: ZEXPL=0.0      ! explicit coefficient
!
REAL, DIMENSION(SIZE(B%XT_FLOOR,1),SIZE(B%XT_FLOOR,2)) :: ZA,& ! lower diag.
                                                      ZB,& ! main  diag.
                                                      ZC,& ! upper diag.
                                                      ZY   ! r.h.s.
!
REAL, DIMENSION(SIZE(B%XT_FLOOR,1)) :: ZTS_FL ! surf. floor temp.  used for rad. exchanges
REAL, DIMENSION(SIZE(B%XT_FLOOR,1)) :: ZTS_FL_CONV ! surf. floor temp. used for conv exchanges
                                               ! used during calculation
REAL, DIMENSION(SIZE(B%XT_FLOOR,1)) :: ZCHTC_IN_FL   ! Indoor floor convec heat transfer coefficient
                                                    ! [W K-1 m-2(bld)]

REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER :: JJ
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FLOOR_LAYER_E_BUDGET',0,ZHOOK_HANDLE)
!
! *Convection heat transfer coefficients [W m-2 K-1]
!  From EP Engineering Reference
ZCHTC_IN_FL(:) = CHTC_UP_DOE(B%XT_FLOOR(:,1), B%XTI_BLD(:))
DO JJ=1,SIZE(ZCHTC_IN_FL)
   ZCHTC_IN_FL(JJ) = MAX(1., ZCHTC_IN_FL(JJ))
ENDDO


!
 CALL LAYER_E_BUDGET_GET_COEF( B%XT_FLOOR, PTSTEP, ZIMPL, B%XHC_FLOOR, B%XTC_FLOOR, B%XD_FLOOR, &
                              ZA, ZB, ZC, ZY )
!
ZTS_FL(:) = B%XT_FLOOR(:,1) 

ZB(:,1) = ZB(:,1) + ZIMPL * (ZCHTC_IN_FL(:)*4./3. + PRADHT_IN(:) * B%XF_FLOOR_MASS(:))


ZY(:,1) = ZY(:,1)  &
   + ZCHTC_IN_FL(:) * (B%XTI_BLD(:) - 1./3. * B%XT_FLOOR(:, 1) * (4* ZEXPL -1))  &
   + B%XF_FLOOR_WIN  (:) * PRAD_WIN_FL(:) + B%XF_FLOOR_WALL (:) * PRAD_WL_FL(:)  &
   + B%XF_FLOOR_ROOF (:) * PRAD_RF_FL(:)                   &
   + PRADHT_IN(:) * B%XF_FLOOR_MASS (:) * (B%XT_MASS(:,1) - ZEXPL * B%XT_FLOOR(:,1))  &
   + PLOAD_FL(:)
!
 CALL LAYER_E_BUDGET( B%XT_FLOOR, PTSTEP, ZIMPL, B%XHC_FLOOR, B%XTC_FLOOR, B%XD_FLOOR, &
                     ZA, ZB, ZC, ZY, PDQS_FL )
!
!*      floor surface temperature used in the implicit formulation
!       ----------------------------------------------------------
ZTS_FL_CONV(:) = 4./3. * ZIMPL * B%XT_FLOOR(:,1) + 1./3. * ZTS_FL(:) * (4 * ZEXPL - 1.)
ZTS_FL(:) = ZEXPL * ZTS_FL(:) + ZIMPL * B%XT_FLOOR(:,1)
!
!*      fluxes with mass and indoor air
!       ----------------------------------------------------------
PRAD_FL_MA(:) = PRADHT_IN(:)     * (ZTS_FL(:) - B%XT_MASS(:,1))
PCONV_FL_BLD(:) = ZCHTC_IN_FL(:) * (ZTS_FL_CONV(:) - B%XTI_BLD (:))
!
!*     Flux between floor and indoor surfaces and air
!       ------------------------------------------------
PFLX_BLD_FL(:) = - PCONV_FL_BLD(:) + B%XF_FLOOR_WIN  (:) * PRAD_WIN_FL(:)  &
       + B%XF_FLOOR_WALL (:) * PRAD_WL_FL(:) + B%XF_FLOOR_ROOF (:) * PRAD_RF_FL(:) &
       + PRADHT_IN(:) * B%XF_FLOOR_MASS(:) * (B%XT_MASS(:,1)  - ZTS_FL(:)) &
       + PLOAD_FL(:)
!
!*     Floor residual energy imbalance for verification
!       ------------------------------------------------
PIMB_FL(:) = PFLX_BLD_FL(:) - PDQS_FL(:)
!
IF (LHOOK) CALL DR_HOOK('FLOOR_LAYER_E_BUDGET',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE FLOOR_LAYER_E_BUDGET
