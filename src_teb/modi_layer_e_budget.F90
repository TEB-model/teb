!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_LAYER_E_BUDGET
INTERFACE
    SUBROUTINE LAYER_E_BUDGET( PT, PTSTEP, PIMPL, PHC, PTC, PD, PA, PB, PC, PY, PDQS )
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT      ! floor layers temperatures
REAL,                 INTENT(IN)    :: PTSTEP  ! time step
REAL, DIMENSION(:,:), INTENT(IN)    :: PHC     ! heat capacity for road layers
REAL, DIMENSION(:,:), INTENT(IN)    :: PTC     ! thermal conductivity for 
REAL, DIMENSION(:,:), INTENT(IN)  :: PD      ! depth of road layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PA 
REAL, DIMENSION(:,:), INTENT(IN)  :: PB 
REAL, DIMENSION(:,:), INTENT(IN)  :: PC
REAL, DIMENSION(:,:), INTENT(IN)  :: PY
REAL, DIMENSION(:), INTENT(OUT) :: PDQS
REAL,               INTENT(IN)  :: PIMPL ! implicitation coefficient
END SUBROUTINE LAYER_E_BUDGET
END INTERFACE
END MODULE MODI_LAYER_E_BUDGET
