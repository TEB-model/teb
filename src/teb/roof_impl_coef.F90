!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!   ###############################################
SUBROUTINE ROOF_IMPL_COEF(T, PTSTEP, PTDEEP_A, PTDEEP_B)
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
!!      V. Masson           * Meteo-France *
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
USE MODD_TEB_n, ONLY : TEB_t
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
TYPE(TEB_t), INTENT(INOUT) :: T
!
REAL              ,   INTENT(IN)  :: PTSTEP      ! time step
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
INTEGER :: ILAYER
!
REAL, DIMENSION(SIZE(PTDEEP_A),SIZE(T%XT_ROOF,2)) :: ZA,& ! lower diag.
                                               ZB,& ! main  diag.
                                               ZC,& ! upper diag.
                                               ZY   ! r.h.s.

REAL, DIMENSION(SIZE(PTDEEP_A))             :: ZDET ! work array
REAL, DIMENSION(SIZE(PTDEEP_A),SIZE(T%XT_ROOF,2)) :: ZW   ! work array
REAL, DIMENSION(SIZE(PTDEEP_A),SIZE(T%XT_ROOF,2)) :: ZT   ! guess of T
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ROOF_IMPL_COEF',0,ZHOOK_HANDLE)
!
!*      1.0    Coefficients of the tridioagonal matrix for heat conduction eq.
!              ---------------------------------------------------------------
!
ILAYER = SIZE(T%XT_ROOF,2)
!
 CALL LAYER_E_BUDGET_GET_COEF(T%XT_ROOF, PTSTEP, ZIMPL, T%XHC_ROOF, T%XTC_ROOF, T%XD_ROOF, &
                              ZA, ZB, ZC, ZY )
!
!-------------------------------------------------------------------------------
!
!*      2.0    Solving of the equation from bottom to top
!              ------------------------------------------
!
! layer at bottom of roof
!
ZDET(:)   = ZB(:,ILAYER)
!
ZT  (:,ILAYER) = ZY(:,ILAYER) / ZDET(:)
!
! internal layers & top layer (but without the external heat flux term)
!
DO JK=ILAYER-1,1,-1
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
!PTDEEP_B = T%XT_ROOF(:,1)
!PTDEEP_A = 0.
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ROOF_IMPL_COEF',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE ROOF_IMPL_COEF
