!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL licence version 2.1.
! See the following links for details:
! https://cecill.info/licences/Licence_CeCILL_V2.1-en.txt
! https://cecill.info/licences/Licence_CeCILL_V2.1-fr.txt
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
