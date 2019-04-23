!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!   ##########################################################################
    SUBROUTINE LAYER_E_BUDGET_GET_COEF(PT, PTSTEP, PIMPL, PHC, PTC, PD, PA, PB, PC, PY)
!   ##########################################################################
!
!!****  *LAYER_E_BUDGET_GET_COEF*  
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
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(IN) :: PT     ! floor layers temperatures
REAL,                 INTENT(IN) :: PTSTEP ! time step
REAL, DIMENSION(:,:), INTENT(IN) :: PHC    ! heat capacity for road layers
REAL, DIMENSION(:,:), INTENT(IN) :: PTC    ! thermal conductivity for road layers
REAL, DIMENSION(:,:), INTENT(IN) :: PD     ! depth of road layers
REAL, DIMENSION(:,:), INTENT(OUT) :: PA 
REAL, DIMENSION(:,:), INTENT(OUT) :: PB 
REAL, DIMENSION(:,:), INTENT(OUT) :: PC
REAL, DIMENSION(:,:), INTENT(OUT)  :: PY
REAL,                 INTENT(IN)   :: PIMPL  ! implicitation coefficient
!
!*      0.2    declarations of local variables
!
REAL :: ZIMPL          ! implicit coefficient
REAL :: ZEXPL          ! explicit coefficient
!
! mean thermal conductivity over distance between 2 layers
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2)) :: ZMTC_O_D
! thermal capacity times layer depth
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2)) :: ZHC_D
!
INTEGER :: ILAYER          ! number of floor layers
INTEGER :: JLAYER          ! loop counter
INTEGER :: JJ              ! loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LAYER_E_BUDGET_GET_COEF',0,ZHOOK_HANDLE)
!
!*      1.     layer thermal properties
!              ------------------------
!
ZIMPL = PIMPL
ZEXPL = 1. - PIMPL
!
ILAYER = SIZE(PT,2)
ZMTC_O_D(:,:) = 0.
!
DO JLAYER=1,ILAYER-1
  DO JJ=1,SIZE(PT,1)
    ZMTC_O_D(JJ,JLAYER) = 2./( PD(JJ,JLAYER)/PTC(JJ,JLAYER) + PD(JJ,JLAYER+1)/PTC(JJ,JLAYER+1) )
    ZHC_D   (JJ,JLAYER) = PHC(JJ,JLAYER) * PD(JJ,JLAYER)
  ENDDO
END DO
!
DO JJ=1,SIZE(PT,1)
  ZHC_D(JJ,ILAYER) = PHC(JJ,ILAYER) * PD(JJ,ILAYER)
ENDDO
!
!-------------------------------------------------------------------------------
!
!*      3.    Surface layer coefficients
!             ------------------------------
!
! no implication for Ts and use of flux from BEM or BLD_E_BUDGET
DO JLAYER = 1, ILAYER
  DO JJ = 1, SIZE(PT,1)
    !
    IF ( JLAYER == 1 ) THEN 
      PA(JJ,JLAYER) = 0.
    ELSE
      PA(JJ,JLAYER) = - ZIMPL * ZMTC_O_D(JJ,JLAYER-1)
    ENDIF
    !
    IF ( JLAYER ==  ILAYER ) THEN
      PC(JJ,JLAYER) = 0. 
    ELSE
      PC(JJ,JLAYER) = - ZIMPL * ZMTC_O_D(JJ,JLAYER)
    ENDIF
    !
    PB(JJ,JLAYER) = ZHC_D(JJ,JLAYER)/PTSTEP
    PY(JJ,JLAYER) = ZHC_D(JJ,JLAYER)/PTSTEP * PT(JJ,JLAYER)
    !
    IF ( JLAYER .GT. 1.) THEN
      PB(JJ,JLAYER) = PB(JJ,JLAYER) + ZIMPL * ZMTC_O_D(JJ,JLAYER-1)
      PY(JJ,JLAYER) = PY(JJ,JLAYER) + ZEXPL * ZMTC_O_D(JJ,JLAYER-1) * ( PT(JJ,JLAYER-1) - PT(JJ,JLAYER) )
    ENDIF
    !
    IF ( JLAYER .LT. ILAYER) THEN
      PB(JJ,JLAYER) = PB(JJ,JLAYER) + ZIMPL * ZMTC_O_D(JJ,JLAYER)
      PY(JJ,JLAYER) = PY(JJ,JLAYER) + ZEXPL * ZMTC_O_D(JJ,JLAYER) * ( PT(JJ,JLAYER+1) - PT(JJ,JLAYER) )
    ENDIF
    !
  ENDDO
  !
ENDDO
!
IF (LHOOK) CALL DR_HOOK('LAYER_E_BUDGET_GET_COEF',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE LAYER_E_BUDGET_GET_COEF
