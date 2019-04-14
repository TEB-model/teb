!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
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
!! P. Samuelsson  10/2014   MEB complements
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
REAL, SAVE       :: XANSMAX = 0.85
!
! Coefficient for calculation of snow fraction over vegetation
!
REAL, SAVE       :: XWSNV
!
! Roughness length of pure snow surface (m)
!
REAL, SAVE       :: XZ0SN   = 0.001
!
! Roughness length for heat of pure snow surface (m)
!
REAL, SAVE       :: XZ0HSN  = 0.0001
!
!--------------------------------------------------------------------------------
! Snow on the ground: PARAMETER
!--------------------------------------------------------------------------------
!
! Critical value of the equivalent water content
! of the snow reservoir for snow fractional coverage and albedo computations
!
REAL, PARAMETER       :: XWCRN      = 10.0   ! (kg m-2) Veg (default value)
REAL, PARAMETER       :: XWCRN_ROOF =  1.0   ! (kg m-2)  Roofs 
REAL, PARAMETER       :: XWCRN_ROAD =  1.0   ! (kg m-2)  Roads
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
REAL, PARAMETER       :: XANS_TODRY_ROOF = 0.008   ! (-)  Roofs
REAL, PARAMETER       :: XANS_TODRY_ROAD = 0.008   ! (-)  Roads
!
REAL, PARAMETER       :: XANS_T_ROOF     = 0.174   ! (-)  Roofs
REAL, PARAMETER       :: XANS_T_ROAD     = 0.174   ! (-)  Roads (alley simul)
!
! Minimum and maximum values of the density of snow 
! for Force-Restore snow option
!
REAL, PARAMETER       :: XRHOSMIN_ROOF = 100.  ! (kg m-3)   Roofs
REAL, PARAMETER       :: XRHOSMIN_ROAD = 100.  ! (kg m-3)   Roads
!
REAL, PARAMETER       :: XRHOSMAX = 300.       ! (kg m-3)   Veg (Default value)
REAL, PARAMETER       :: XRHOSMAX_ROOF = 300.  ! (kg m-3)   Roofs
REAL, PARAMETER       :: XRHOSMAX_ROAD = 350.  ! (kg m-3)   Roads
!
!--------------------------------------------------------------------------------
END MODULE MODD_SNOW_PAR
