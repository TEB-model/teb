!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ###############
      MODULE MODD_BEM_CST      
!     ###############
!
!!****  *MODD_BEM_CST* - declaration of Technical constants for Building Energy Model
!!
!!    PURPOSE
!!    -------
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
!!      G. Pigeon    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2012
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE 
REAL, PARAMETER :: XWIN_SW_MAX = 150. ! Maximum shortwave radiation received by
!                                     ! windows before activation of shading 
!                                     ! (if available and pertinent).
REAL, PARAMETER :: XCOMFORT_TEMP_FOR_SHADING_USE = 297.16
!                                     ! temperature above which one considers
!                                     ! the use of shading protections
!
END MODULE MODD_BEM_CST

