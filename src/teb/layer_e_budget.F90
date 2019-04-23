!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!   ##########################################################################
    SUBROUTINE LAYER_E_BUDGET( PT, PTSTEP, PIMPL, PHC, PTC, PD, PA, PB, PC, PY, PDQS )
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
!!                     08/10 (G. Pigeon) computation of residual of energy balance 
!!                                       modification of the limit condition for
!!                                       the deep temp. from the deep road temp.
!!                                       to zero flux condition. idem for sfce T
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODI_TRIDIAG_GROUND
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT      ! floor layers temperatures
REAL,                 INTENT(IN)    :: PTSTEP  ! time step
REAL, DIMENSION(:,:), INTENT(IN)    :: PHC     ! heat capacity for road layers
REAL, DIMENSION(:,:), INTENT(IN)    :: PTC     ! thermal conductivity for 
                                               !road layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD      ! depth of road layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PA 
REAL, DIMENSION(:,:), INTENT(IN)  :: PB 
REAL, DIMENSION(:,:), INTENT(IN)  :: PC
REAL, DIMENSION(:,:), INTENT(IN)  :: PY
REAL, DIMENSION(:), INTENT(OUT) :: PDQS
REAL,               INTENT(IN)  :: PIMPL ! implicitation coefficient
!
!*      0.2    declarations of local variables
!
!
REAL :: ZIMPL          ! implicit coefficient
REAL :: ZEXPL          ! explicit coefficient
!
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2)) :: ZX   ! solution
!
REAL, DIMENSION(SIZE(PT,1)) :: ZEI  ! internal energy of floor at t
REAL, DIMENSION(SIZE(PT,1)) :: ZPEI ! internal energy of floor at time t+
!
INTEGER :: ILAYER          ! number of floor layers
INTEGER :: JLAYER          ! loop counter
INTEGER :: JJ              ! loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LAYER_E_BUDGET',0,ZHOOK_HANDLE)
!
ILAYER = SIZE(PT,2)
!
ZIMPL = PIMPL
ZEXPL = 1.-PIMPL
!
!-------------------------------------------------------------------------------
!
!*      1.    Preliminaries : internal energy of floor at the current time step
!             -----------------------------------------------------------------
!
ZEI(:) = 0.
DO JLAYER=1,ILAYER
  DO JJ=1,SIZE(PT,1)
     ZEI(JJ)=ZEI(JJ) + ( PHC(JJ,JLAYER)*PD(JJ,JLAYER)*PT(JJ,JLAYER) )
  ENDDO
END DO
!
!-------------------------------------------------------------------------------
!
!*     2.     Tri-diagonal system resolution
!              ------------------------------
!
 CALL TRIDIAG_GROUND(PA,PB,PC,PY,ZX)
!
DO JLAYER=1,ILAYER
  PT(:,JLAYER) = ZX(:,JLAYER)
END DO
!
!*      3.     heat storage inside floor and flux toward the floor
!              ---------------------------------------------------
!
!       3.1    internal energy of the floor at the next time step
!              --------------------------------------------------
!
ZPEI(:) = 0.0
DO JLAYER=1,ILAYER
  DO JJ=1,SIZE(PT,1)
     ZPEI(JJ) = ZPEI(JJ)+ ( PHC(JJ,JLAYER)*PD(JJ,JLAYER)*PT(JJ,JLAYER) )
  ENDDO
END DO
!
!        3.2   heat storage flux inside floor 
!              ------------------------------
!
PDQS(:)=(ZPEI(:)-ZEI(:))/PTSTEP
!
IF (LHOOK) CALL DR_HOOK('LAYER_E_BUDGET',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE LAYER_E_BUDGET
