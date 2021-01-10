!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
SUBROUTINE INI_CTURBS
!!
!!****     *INI_CTURBS*  - routine to initialize the turbulence scheme 
!!                        constants.
!!
!!      PURPOSE
!!      -------
!         The purpose of this routine is to initialize the turbulence 
!       scheme constants that are stored in module MODD_CTURB
!
!!      METHOD
!!      ------
!!        The constants are set to their numerical values
!!
!!      AUTHOR
!!      ------
!!        V. Masson         * Meteo-France *
!!
!!      MODIFICATIONS
!!      -------------
!! --------------------------------------------------------------------------
!
!*        0. DECLARATIONS
!            ------------
!
USE MODD_CSTS,        ONLY : XKARMAN
USE MODD_CANOPY_TURB, ONLY : XTKEMIN, XCEP, XCED, XALPSBL, XA0, XCMFS, XASBL,&
                               XCSHF, XCTP  
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!  ---------------------------------------------------------------------------
!
!         1. SETTING THE NUMERICAL VALUES
!            ----------------------------
!
!         1.1 Constant for dissipation of Tke
!

REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('INI_CTURBS',0,ZHOOK_HANDLE)
XCED  = 0.845
!       Redelsperger-Sommeria (1981) = 0.70
!       Schmidt-Schumann      (1989) = 0.845
!       Cheng-Canuto-Howard   (2002) = 0.845
!
!
!         1.2 Constant for wind pressure-correlations
!
XCEP  = 2.11
!       Redelsperger-Sommeria (1981) = 4.
!       Schmidt-Schumann      (1989) = 3.5
!       Cheng-Canuto-Howard   (2002) = 2.11
!
!
!         1.3 Constant a0 for wind pressure-correlations
!
XA0   = 0.6
!       Redelsperger-Sommeria (1981) = 0.6
!       Schmidt-Schumann      (1989) = 0.55
!       Cheng-Canuto-Howard   (2002) = 0.6
!
!
!         1.4 Constant for temperature and vapor pressure-correlations
!
XCTP  = 4.65
!       Redelsperger-Sommeria (1981) = 4.
!       Schmidt-Schumann      (1989) = 3.25
!       Cheng-Canuto-Howard   (2002) = 4.65
!
!
!         1.5 Value related to the TKE universal function within SBL
!
XALPSBL = 4.63
!       Redelsperger et al 2001     = 4.63
!       Wyngaard et al. 1974        = 3.75
!       Stull 1988                  = 4.75
!
!
!
!         2. Derivated constants
!            -------------------
!
!         2.1 Constant in fluxes equations
!
XCMFS= 2./3./XCEP*(1.-XA0)   !Constant for the momentum flux due to shear (RS)
!
! Redelsperger-Sommeria (1981) ......... 0.066
! Schmidt-Schumann      (1989) ......... 0.086
! Cheng-Canuto-Howard   (2002) ......... 0.126
!
XCSHF= 2./3./XCTP            !Constant for the sensible heat flux(RS)
!
! Redelsperger-Sommeria (1981) ......... 0.167
! Schmidt-Schumann      (1989) ......... 0.204
! Cheng-Canuto-Howard   (2002) ......... 0.143
!
!
!         2.2 Value related to the TKE universal function within SBL
!
!
XASBL   = 0.5*( XALPSBL**(3./2.)*XKARMAN*XCED + XKARMAN/SQRT(XALPSBL)/XCMFS )
!       Redelsperger et al 2001
!
!
!
!         3. MINIMUM VALUES 
!            --------------
!
XTKEMIN=1.E-4
IF (LHOOK) CALL DR_HOOK('INI_CTURBS',1,ZHOOK_HANDLE)
!

END SUBROUTINE INI_CTURBS
