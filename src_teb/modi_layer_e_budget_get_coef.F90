!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_LAYER_E_BUDGET_GET_COEF
INTERFACE
    SUBROUTINE LAYER_E_BUDGET_GET_COEF(PT, PTSTEP, PIMPL, PHC, PTC, PD, PA, PB, PC, PY)
REAL, DIMENSION(:,:), INTENT(IN) :: PT     ! floor layers temperatures
REAL,                 INTENT(IN) :: PTSTEP ! time step
REAL, DIMENSION(:,:), INTENT(IN) :: PHC    ! heat capacity for road layers
REAL, DIMENSION(:,:), INTENT(IN) :: PTC    ! thermal conductivity for road layers
REAL, DIMENSION(:,:), INTENT(IN) :: PD     ! depth of road layers
REAL, DIMENSION(:,:), INTENT(OUT) :: PA 
REAL, DIMENSION(:,:), INTENT(OUT) :: PB 
REAL, DIMENSION(:,:), INTENT(OUT) :: PC
REAL, DIMENSION(:,:), INTENT(OUT)  :: PY
REAL,                 INTENT(IN)   :: PIMPL  ! implicitation coefficient
END SUBROUTINE LAYER_E_BUDGET_GET_COEF
END INTERFACE
END MODULE MODI_LAYER_E_BUDGET_GET_COEF
