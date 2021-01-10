!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ###############
      MODULE MODD_OCEAN_CSTS      
!     ###############
!
!!****  *MODD_OCEAN_CSTS* - declaration of Physic constants for ocean
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare  the 
!     Physics constants for ocean.    
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!          
!!    AUTHOR
!!    ------
!!      C. Lebeaupin Brossier   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2008
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE 
!
REAL,SAVE :: XRHOSW,XRHOSWREF   ! densities of seawater
REAL,SAVE :: XCPSW              ! Cp (for seawater)
!
REAL,SAVE :: XGAMA              ! fraction of Coriolis parameter in time 
!
REAL,SAVE :: XZCE,XCKL          ! turbulence constants
!
REAL,SAVE :: XR                 ! from Jerlov's (1976)
REAL,SAVE :: XD1,XD2            ! optical water C type I
REAL,SAVE :: XALBEDOSW          ! Albedo of seawater
!
!
END MODULE MODD_OCEAN_CSTS
