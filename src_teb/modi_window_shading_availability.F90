!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_WINDOW_SHADING_AVAILABILITY
INTERFACE
      SUBROUTINE WINDOW_SHADING_AVAILABILITY(OSHADE, PTI_BLD, PTCOOL_TARGET,OSHADE_POSSIBLE)
LOGICAL, DIMENSION(:), INTENT(IN)  :: OSHADE          ! TRUE if solar protections exist
REAL,    DIMENSION(:), INTENT(IN)  :: PTI_BLD         ! indoor air temperature
REAL,    DIMENSION(:), INTENT(IN)  :: PTCOOL_TARGET   ! Cooling setpoint of HVAC system
LOGICAL, DIMENSION(:), INTENT(OUT) :: OSHADE_POSSIBLE ! TRUE if solar protections 
END SUBROUTINE WINDOW_SHADING_AVAILABILITY
END INTERFACE
END MODULE MODI_WINDOW_SHADING_AVAILABILITY
