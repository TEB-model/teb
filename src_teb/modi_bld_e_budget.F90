!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_BLD_E_BUDGET
INTERFACE
    SUBROUTINE BLD_E_BUDGET( OTI_EVOL, PTSTEP, PBLD, PWALL_O_HOR,      &
                             PRHOA, PT_ROOF, PT_WALL, PTI_BLD, PTS_FLOOR )  
LOGICAL,              INTENT(IN)   :: OTI_EVOL      ! true --> internal temp. of
REAL,                 INTENT(IN)   :: PTSTEP        ! time step
REAL, DIMENSION(:),   INTENT(IN)   :: PBLD          ! building fraction
REAL, DIMENSION(:),   INTENT(IN)   :: PWALL_O_HOR   ! wall surf. / hor. surf.
REAL, DIMENSION(:),   INTENT(IN)   :: PRHOA         ! air density
REAL, DIMENSION(:,:), INTENT(IN)   :: PT_ROOF       ! roof layers temperatures
REAL, DIMENSION(:,:), INTENT(IN)   :: PT_WALL       ! wall layers temperatures
REAL, DIMENSION(:),   INTENT(INOUT):: PTI_BLD       ! building air temperature
REAL, DIMENSION(:),   INTENT(IN)  :: PTS_FLOOR     ! floor surface temperature
END SUBROUTINE BLD_E_BUDGET
END INTERFACE
END MODULE MODI_BLD_E_BUDGET
