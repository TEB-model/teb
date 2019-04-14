!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_WIND_THRESHOLD
INTERFACE
    FUNCTION WIND_THRESHOLD(PWIND,PUREF) RESULT(PWIND_NEW)
REAL, DIMENSION(:), INTENT(IN)   :: PWIND      ! wind
REAL, DIMENSION(:), INTENT(IN)   :: PUREF      ! forcing level
REAL, DIMENSION(SIZE(PWIND))     :: PWIND_NEW  ! modified wind
END FUNCTION WIND_THRESHOLD
END INTERFACE
END MODULE MODI_WIND_THRESHOLD
