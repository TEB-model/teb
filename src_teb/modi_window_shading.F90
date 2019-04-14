!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_WINDOW_SHADING
INTERFACE
SUBROUTINE WINDOW_SHADING(PSHGC, PSHGC_SH, O_SHADE, PALB_WALL,      &
                          PABS_WIN, PABS_WINSH, PALB_WIN, PTRAN_WIN )
REAL,    DIMENSION(:), INTENT(IN) :: PSHGC     !Window solar heat gain coefficient
REAL,    DIMENSION(:), INTENT(IN) :: PSHGC_SH  !Window + shading solar heat gain coefficient
LOGICAL, DIMENSION(:), INTENT(IN) :: O_SHADE   !use of shadings TRUE -> shadings ;
REAL, DIMENSION(:), INTENT(IN)  :: PALB_WALL !albedo of the wall     
REAL, DIMENSION(:), INTENT(IN)  :: PABS_WIN  !Window absorptivity
REAL, DIMENSION(:), INTENT(OUT) :: PABS_WINSH!Window absorptivity after shading
REAL, DIMENSION(:), INTENT(OUT) :: PALB_WIN  !Albedo of the ensemble window + shading
REAL, DIMENSION(:), INTENT(INOUT) :: PTRAN_WIN !Window transmitivity
END SUBROUTINE WINDOW_SHADING
END INTERFACE
END MODULE MODI_WINDOW_SHADING
