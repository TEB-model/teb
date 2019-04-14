!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ###############################################
SUBROUTINE ROOF_IMPL_COEF(PTSTEP, KROOF_LAYER, PD_ROOF, PTC_ROOF, PHC_ROOF, PT_ROOF, PTDEEP_A,PTDEEP_B)
!   ###############################################
!
!!
!!    PURPOSE
!!    -------
!
!     Computes the corfficients for implicitation of upper
!     roof layer with what is above
!         
!     
!!**  METHOD
!     ------
!
!  One computes a guess assuming a zero flux condition at the base
!  of the roof. One solves the half part of the tridiagonal matrix
!  fromm bottom to top.
!
!!     The classical tridiagonal algorithm is used to invert the 
!!     implicit operator (from bottom to top only). Its matrix is given by:
!!
!!     (  b(1)      c(1)      0        0        0         0        0        0  )
!!     (  a(2)      b(2)     c(2)      0  ...    0        0        0        0  ) 
!!     (   0        a(3)     b(3)     c(3)       0        0        0        0  ) 
!!      .......................................................................
!!     (   0   ...   0      a(k)      b(k)     c(k)       0   ...  0        0  ) 
!!      .......................................................................
!!     (   0         0        0        0        0 ...  a(n-1)   b(n-1)   c(n-1))
!!     (   0         0        0        0        0 ...     0      a(n)     b(n) )
!!
!!
!!    AUTHOR
!!    ------
!!
!!	V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2013
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_LAYER_E_BUDGET_GET_COEF
!
IMPLICIT NONE
!
!
!*      0.1    Declarations of arguments
!
REAL              ,   INTENT(IN)  :: PTSTEP      ! time step
INTEGER           ,   INTENT(IN)  :: KROOF_LAYER ! number of roof layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_ROOF   ! thickness            of each layer
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_ROOF  ! thermal conductivity of each layer
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_ROOF  ! heat capacity        of each layer
REAL, DIMENSION(:,:), INTENT(IN)  :: PT_ROOF   ! temperature          of each layer
REAL, DIMENSION(:),   INTENT(OUT) :: PTDEEP_A, PTDEEP_B     
                                              ! Deep soil temperature (prescribed)
!                                      PTDEEP_A = Deep soil temperature
!                                                 coefficient depending on flux
!                                      PTDEEP_B = Deep soil temperature (prescribed)
!                                               which models heating/cooling from
!                                               below the diurnal wave penetration
!                                               (surface temperature) depth. If it
!                                               is FLAGGED as undefined, then the zero
!                                               flux lower BC is applied.
!                                      Tdeep = PTDEEP_B + PTDEEP_A * PDEEP_FLUX
!                                              (with PDEEP_FLUX in W/m2)
!
!*      0.2    Local variables
!
REAL :: ZIMPL = 1.0        ! implicit coefficient
INTEGER :: JK              ! loop counter
!
REAL, DIMENSION(SIZE(PTDEEP_A),KROOF_LAYER) :: ZA,& ! lower diag.
                                               ZB,& ! main  diag.
                                               ZC,& ! upper diag.
                                               ZY   ! r.h.s.

REAL, DIMENSION(SIZE(PTDEEP_A))             :: ZDET ! work array
REAL, DIMENSION(SIZE(PTDEEP_A),KROOF_LAYER) :: ZW   ! work array
REAL, DIMENSION(SIZE(PTDEEP_A),KROOF_LAYER) :: ZT   ! guess of T
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ROOF_IMPL_COEF',0,ZHOOK_HANDLE)
!
!*      1.0    Coefficients of the tridioagonal matrix for heat conduction eq.
!              ---------------------------------------------------------------
!
 CALL LAYER_E_BUDGET_GET_COEF( PT_ROOF, PTSTEP, ZIMPL, PHC_ROOF, PTC_ROOF, PD_ROOF, &
                              ZA, ZB, ZC, ZY )
!
!-------------------------------------------------------------------------------
!
!*      2.0    Solving of the equation from bottom to top
!              ------------------------------------------
!
! layer at bottom of roof
!
ZDET(:)   = ZB(:,KROOF_LAYER)
!
ZT  (:,KROOF_LAYER) = ZY(:,KROOF_LAYER) / ZDET(:)
!
! internal layers & top layer (but without the external heat flux term)
!
DO JK=KROOF_LAYER-1,1,-1
  ZW  (:,JK)  = ZA(:,JK+1)/ZDET(:)
  ZDET(:)     = ZB(:,JK  ) - ZC(:,JK)*ZW(:,JK)
  ZT  (:,JK)  = ( ZY(:,JK) - ZC(:,JK)*ZT(:,JK+1) ) / ZDET(:) ! + FLUX / ZDET
  !                                                          ! for layer 1
  !                                                          ! because the external
  !                                                          ! flux would be
  !                                                          ! included in the Y
  !                                                          ! term
END DO
!
! Implicit coefficients for the heat flux
!
PTDEEP_B = ZT  (:,1)
PTDEEP_A = 1. / ZDET(:) 
!
!* The following lines are here if you want to test the explicit coupling
!PTDEEP_B = PT_ROOF(:,1)
!PTDEEP_A = 0.
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ROOF_IMPL_COEF',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE ROOF_IMPL_COEF
