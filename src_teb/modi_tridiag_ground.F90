!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_TRIDIAG_GROUND
INTERFACE
       SUBROUTINE TRIDIAG_GROUND(PA,PB,PC,PY,PX)
REAL,    DIMENSION(:,:), INTENT(IN)  :: PA  ! lower diag. elements of A matrix
REAL,    DIMENSION(:,:), INTENT(IN)  :: PB  ! main  diag. elements of A matrix
REAL,    DIMENSION(:,:), INTENT(IN)  :: PC  ! upper diag. elements of A matrix
REAL,    DIMENSION(:,:), INTENT(IN)  :: PY  ! r.h.s. term   
REAL,    DIMENSION(:,:), INTENT(OUT) :: PX  ! solution of A.X = Y 
END SUBROUTINE TRIDIAG_GROUND
END INTERFACE
END MODULE MODI_TRIDIAG_GROUND
