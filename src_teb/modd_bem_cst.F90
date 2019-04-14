!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

