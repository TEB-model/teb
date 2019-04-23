!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!* The modules below are just proxi of a debugging code internal
! to Meteo-France and IFS softwares. These proxi have no use at all
! except the one to allow the compiling of the routines without having
! to change them from the ones in the SURFEX platform (containing TEB 
! and other surface schemes).
!
MODULE PARKIND1
INTEGER, PARAMETER :: JPRB = 4
END MODULE PARKIND1
!
MODULE YOMHOOK
LOGICAL :: LHOOK = .FALSE.
CONTAINS
SUBROUTINE DR_HOOK(HSTRING,KINTEGER,PREAL)
USE PARKIND1
  CHARACTER(LEN=*) :: HSTRING
  INTEGER          :: KINTEGER
  REAL(KIND=JPRB)  :: PREAL
print*, HSTRING,KINTEGER
PREAL=0.
END SUBROUTINE DR_HOOK
END MODULE YOMHOOK
