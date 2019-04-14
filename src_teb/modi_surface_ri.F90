!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_SURFACE_RI
INTERFACE
    SUBROUTINE SURFACE_RI(PTG, PQS, PEXNS, PEXNA, PTA, PQA,   &
                               PZREF, PUREF, PDIRCOSZW, PVMOD, PRI )  
REAL, DIMENSION(:), INTENT(IN)    :: PTG      ! surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PQS      ! surface specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS    ! surface exner function
REAL, DIMENSION(:), INTENT(IN)    :: PTA      ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA      ! specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PEXNA    ! exner function
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD    ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PZREF    ! reference height of the first
REAL, DIMENSION(:), INTENT(IN)    :: PUREF    ! reference height of the wind
REAL, DIMENSION(:), INTENT(IN)    :: PDIRCOSZW! Cosine of the angle between
REAL, DIMENSION(:), INTENT(OUT)   :: PRI      ! Richardson number
END SUBROUTINE SURFACE_RI
END INTERFACE
END MODULE MODI_SURFACE_RI
