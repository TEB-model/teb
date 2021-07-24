!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL licence
!SFX_LIC version 2.1. See Licence_CeCILL_V2.1-en.txt and Licence_CeCILL_V2.1-fr.txt  
!SFX_LIC for details.
!     #############################################################
SUBROUTINE WINDOW_DATA_STRUCT(KI,PSHGC, PU_WIN, PALB_WIN, PABS_WIN, PUGG_WIN, PTRAN_WIN)
!     #############################################################
!
!
!!****  *WINDOW_DATA*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the radiative properties of the window
!         
!     
!!**  METHOD
!!     ------
!! uses U_WIN (window conductivity) and SHGC (window solar heat gain coef)
!! to compute the normal incident transmittance/reflectance
!! a factor to normal incident tran/ref is applied to take into account all
!! possible incidence between 18deg and 72deg
!!
!!
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
!! Simple Window Model from EnergyPlus Engineering Reference V7, 2012, p. 217.
!!      
!!    AUTHOR
!!    ------
!!
!!      B. Bueno           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    11/10
!!      G. Pigeon & J. Le Bras  09/12 apply new factor to normal incident tran/ref
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!

!
!
USE MODI_WINDOW_DATA
!
USE MODD_BEM_n, ONLY : BEM_t
!
IMPLICIT NONE
!
INTEGER,             INTENT(IN)  :: KI       ! number of points
REAL, DIMENSION(KI), INTENT(IN)  :: PSHGC    ! solar heat gain coef. of windows
REAL, DIMENSION(KI), INTENT(IN)  :: PU_WIN   ! window U-factor [K m W-2]
REAL, DIMENSION(KI), INTENT(OUT) :: PALB_WIN ! window albedo
REAL, DIMENSION(KI), INTENT(OUT) :: PABS_WIN ! window absortance
REAL, DIMENSION(KI), INTENT(OUT) :: PUGG_WIN ! window glass-to-glass U-factor [W m-2 K-1]
REAL, DIMENSION(KI), INTENT(OUT) :: PTRAN_WIN! window transmittance (-)
!
!
!*      0.2    Declarations of local variables 
!
TYPE(BEM_t) :: B
!
!
!------------------------------------------------------------------------------
ALLOCATE(B%XSHGC    (1))
ALLOCATE(B%XU_WIN   (1))
ALLOCATE(B%XALB_WIN (1))
ALLOCATE(B%XABS_WIN (1))
ALLOCATE(B%XUGG_WIN (1))
ALLOCATE(B%XTRAN_WIN(1))
!------------------------------------------------------------------------------
!
! put input variables in structures
B%XSHGC  = PSHGC    ! solar heat gain coef. of windows
B%XU_WIN = PU_WIN   ! window U-factor [K m W-2]
!------------------------------------------------------------------------------
!
CALL WINDOW_DATA(KI,B)
!
!------------------------------------------------------------------------------
!
! put computed variables in output variables
PALB_WIN = B%XALB_WIN ! window albedo
PABS_WIN = B%XABS_WIN ! window absortance
PUGG_WIN = B%XUGG_WIN ! window glass-to-glass U-factor [W m-2 K-1]
PTRAN_WIN= B%XTRAN_WIN! window transmittance (-)!
!------------------------------------------------------------------------------
DEALLOCATE(B%XSHGC    )
DEALLOCATE(B%XU_WIN   )
DEALLOCATE(B%XALB_WIN )
DEALLOCATE(B%XABS_WIN )
DEALLOCATE(B%XUGG_WIN )
DEALLOCATE(B%XTRAN_WIN)
!------------------------------------------------------------------------------
!
END SUBROUTINE WINDOW_DATA_STRUCT
