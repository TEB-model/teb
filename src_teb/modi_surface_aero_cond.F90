!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_SURFACE_AERO_COND
INTERFACE
    SUBROUTINE SURFACE_AERO_COND(PRI, PZREF, PUREF, PVMOD, PZ0,&
                                     PZ0H, PAC, PRA, PCH           ) 
REAL, DIMENSION(:), INTENT(IN)    :: PRI      ! Richardson number
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD    ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PZREF    ! reference height of the first
REAL, DIMENSION(:), INTENT(IN)    :: PUREF    ! reference height of the wind
REAL, DIMENSION(:), INTENT(IN)    :: PZ0      ! roughness length for momentum
REAL, DIMENSION(:), INTENT(IN)    :: PZ0H     ! roughness length for heat
REAL, DIMENSION(:), INTENT(OUT)   :: PAC      ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PRA      ! aerodynamical resistance
REAL, DIMENSION(:), INTENT(OUT)   :: PCH      ! drag coefficient for heat
END SUBROUTINE SURFACE_AERO_COND
END INTERFACE
END MODULE MODI_SURFACE_AERO_COND
