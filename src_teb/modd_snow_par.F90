!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ######################
      MODULE MODD_SNOW_PAR
!     ######################
!
!!****  *MODD_SNOW_PAR* - declaration of parameters related
!!                          to the snow parameterization
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     parameters related to the surface parameterization of snow.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004                    
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!--------------------------------------------------------------------------------
! Snow on the ground: Given in ini_surf_csts and/or in NAM_SURF_CSTS
!--------------------------------------------------------------------------------
!
! Snow emissivity:
!
REAL, SAVE       :: XEMISSN = 1.0
!
! Minimum and maximum values of the albedo of snow:
!
REAL, SAVE       :: XANSMIN = 0.50
REAL, SAVE       :: XANSMAX = 0.85
!
! Minimum and maximum values of the albedo of permanet snow/ice:
!
REAL, SAVE       :: XAGLAMIN = 0.8 
REAL, SAVE       :: XAGLAMAX = 0.85
! 
! Prescribed ice albedo in 3 spectral bands for glacier simulation with CROCUS scheme.
REAL, SAVE       :: XALBICE1,XALBICE2,XALBICE3
!

! Density threshold for ice detection in CROCUS scheme.
REAL, SAVE       :: XRHOTHRESHOLD_ICE

!for ageing effects
REAL, SAVE      :: XVAGING_NOGLACIER, XVAGING_GLACIER

! Height (m) of aged snow in glacier case (allows Pn=1)
!
REAL, SAVE       :: XHGLA
! 
! Coefficient for calculation of snow fraction over vegetation
!
REAL, SAVE       :: XWSNV  = 5.0
!
! Roughness length of pure snow surface (m)
!
REAL, SAVE       :: XZ0SN   = 0.001
!
! Roughness length for heat of pure snow surface (m)
!
REAL, SAVE       :: XZ0HSN  = 0.0001
!
! Roughness length ratio between ice and snow
REAL, SAVE       :: XZ0ICEZ0SNOW
!
! Logicals for snowdrift and sublimation
LOGICAL,SAVE    :: LSNOWDRIFT, LSNOWDRIFT_SUBLIM

! Snow Melt timescale with D95 (s): needed to prevent time step 
! dependence of melt when snow fraction < unity.
!
REAL, SAVE       :: XTAU_SMELT = 300.
!
!--------------------------------------------------------------------------------
! Snow on the ground: PARAMETER
!--------------------------------------------------------------------------------
!
! Critical value of the equivalent water content
! of the snow reservoir for snow fractional coverage and albedo computations
!
REAL, PARAMETER       :: XWCRN = 10.0   ! (kg m-2) Veg (default value)
REAL, PARAMETER       :: XWCRN_EXPL =  1.0   ! (kg m-2) Veg explicit
REAL, PARAMETER       :: XWCRN_ROOF =  1.0   ! (kg m-2)  Roofs 
REAL, PARAMETER       :: XWCRN_ROAD =  1.0   ! (kg m-2)  Roads
REAL, PARAMETER       :: XWCRN_VEG  =  1.0   ! (kg m-2)  Urban veg
!
!
! Critical value of snow emissivity
!
REAL, PARAMETER       :: XEMCRIN = 0.98
!
! Minimum and maximum values of the albedo of snow:
!
REAL, PARAMETER       :: XANSMIN_ROOF = 0.30 ! (-)   Roofs
REAL, PARAMETER       :: XANSMIN_ROAD = 0.15 ! (-)   Roads
!
REAL, PARAMETER       :: XANSMAX_ROOF = 0.85 ! (-)   Roofs
REAL, PARAMETER       :: XANSMAX_ROAD = 0.85 ! (-)   Roads
!
! Snow aging coefficients (albedo and Force-Restore density):
!
REAL, PARAMETER       :: XANS_TODRY    = 0.008     ! (-) Veg (default value)
REAL, PARAMETER       :: XANS_TODRY_ROOF = 0.008   ! (-)  Roofs
REAL, PARAMETER       :: XANS_TODRY_ROAD = 0.008   ! (-)  Roads
!
REAL, PARAMETER       :: XANS_T        = 0.240     ! (-) Veg (default value)
REAL, PARAMETER       :: XANS_T_ROOF     = 0.174   ! (-)  Roofs
REAL, PARAMETER       :: XANS_T_ROAD     = 0.174   ! (-)  Roads (alley simul)
!
! Minimum and maximum values of the density of snow 
! for Force-Restore snow option
!
REAL, PARAMETER       :: XRHOSMIN = 100.       ! (kg m-3)   Veg (Default value)
REAL, PARAMETER       :: XRHOSMIN_ROOF = 100.  ! (kg m-3)   Roofs
REAL, PARAMETER       :: XRHOSMIN_ROAD = 100.  ! (kg m-3)   Roads
!
REAL, PARAMETER       :: XRHOSMAX = 300.       ! (kg m-3)   Veg (Default value)
REAL, PARAMETER       :: XRHOSMAX_ROOF = 300.  ! (kg m-3)   Roofs
REAL, PARAMETER       :: XRHOSMAX_ROAD = 350.  ! (kg m-3)   Roads
!
! Minimum and maximum values of the density of snow 
! for ISBA-ES snow option
!
REAL, PARAMETER       :: XRHOSMIN_ES =  50.  ! (kg m-3)
REAL, PARAMETER       :: XRHOSMAX_ES = 750.  ! (kg m-3)
!
! ISBA-ES Critical snow depth at which snow grid thicknesses constant
!
REAL, PARAMETER                      :: XSNOWCRITD = 0.03  ! (m)
!                                       
! ISBA-ES Minimum total snow depth for model 
!
! REAL, PARAMETER                      :: XSNOWDMIN = 0.000001  ! (m)
REAL, PARAMETER                      :: XSNOWDMIN = 4E-5  ! (m)
!                                       
! Maximum Richardson number limit for very stable conditions using the ISBA-ES 'RIL' option
!
REAL, PARAMETER                      :: X_RI_MAX = 0.20
!                                       
! ISBA-ES Maximum snow liquid water holding capacity (fraction by mass) parameters:
!
REAL, PARAMETER                      :: XWSNOWHOLDMAX2   = 0.10  ! (-) 
REAL, PARAMETER                      :: XWSNOWHOLDMAX1   = 0.03  ! (-)
REAL, PARAMETER                      :: XSNOWRHOHOLD     = 200.0 ! (kg/m3)
!

!--------------------------------------------------------------------------------
! Calibration coefficients for CROCUS albedo computation
!--------------------------------------------------------------------------------
! for grains effects:
! REAL, PARAMETER :: XD1=1., XD2=3., XD3=4., XX=99.,XVALB2=.96, XVALB3=1.58,&
!                     XVALB4=.94,XVALB5=.95,XVALB6=15.4,XVALB7=346.3, XVALB8=32.31,  &
!                     XVALB9=.88, XVALB10=.175,XVALB11=.7,XVDIOP1=2.3E-3, &
!                     XVRPRE1=.5,XVRPRE2=1. 
! ! for ageing effects:
! REAL, PARAMETER :: XVAGING_SOIL=90. , XVAGING_GLACIER=900. , XVPRES1=87000.

! modifs SM 20110805 tests SIberie - albedo
REAL, PARAMETER :: XD1=1., XD2=3., XD3=4., XX=99.,XVALB2=.96, XVALB3=1.58,&
                   XVALB4=.92,XVALB5=.90,XVALB6=15.4,XVALB7=346.3, XVALB8=32.31,&  
                   XVALB9=.88,XVALB10=.200,XVALB11=.6,XVDIOP1=2.3E-3, &
                   XVRPRE1=.5,XVRPRE2=1.5
! for ageing effects:
REAL, PARAMETER :: XVPRES1=87000.
! REAL, PARAMETER :: XZ_SWE_ALB=5.
! end modifs SM 20110805 tests SIberie - albedo


! for spectral distribution and thickness effects
REAL, PARAMETER :: XVSPEC1=.71, XVSPEC2=.21 , XVSPEC3=.08
! modif SM 20110519
!REAL, PARAMETER :: XVSPEC1=.68, XVSPEC2=.25 , XVSPEC3=.07
! end modif SM 20110519
! for thickness effects
REAL, PARAMETER :: XVW1=.80, XVW2=.20 , XVD1=.02, XVD2=.01
!--------------------------------------------------------------------------------

END MODULE MODD_SNOW_PAR












