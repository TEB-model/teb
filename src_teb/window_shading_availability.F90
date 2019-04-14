!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ###########################################################################################################
      SUBROUTINE WINDOW_SHADING_AVAILABILITY(OSHADE, PTI_BLD, PTCOOL_TARGET,OSHADE_POSSIBLE)
!     ###########################################################################################################
!
!!****  *WINDOW_SHADING_AVAILABILITY* 
!!
!!    PURPOSE
!!    -------
!!* checks if one encouters the conditions for closing the windows at day
!  (needs to be equiped of solar protections and not to be confortable inside)
!  
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	G. Pigeon   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2011
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_BEM_CST, ONLY : XCOMFORT_TEMP_FOR_SHADING_USE
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
LOGICAL, DIMENSION(:), INTENT(IN)  :: OSHADE          ! TRUE if solar protections exist
REAL,    DIMENSION(:), INTENT(IN)  :: PTI_BLD         ! indoor air temperature
REAL,    DIMENSION(:), INTENT(IN)  :: PTCOOL_TARGET   ! Cooling setpoint of HVAC system
LOGICAL, DIMENSION(:), INTENT(OUT) :: OSHADE_POSSIBLE ! TRUE if solar protections 
!                                                     ! are likely to be used today if
!                                                     ! solar irradiance is enough
!
!*       0.2   Declarations of local variables
!
!
!-------------------------------------------------------------------------------
!
!
!
OSHADE_POSSIBLE(:) = OSHADE(:)                                   .AND. &
          (       PTI_BLD(:) .GE. XCOMFORT_TEMP_FOR_SHADING_USE        &
             .OR. PTI_BLD(:) .EQ. PTCOOL_TARGET(:)              )
!
END SUBROUTINE WINDOW_SHADING_AVAILABILITY
