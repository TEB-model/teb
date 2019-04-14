!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_DAY_OF_WEEK
INTERFACE
SUBROUTINE DAY_OF_WEEK(PYEAR, PMONTH, PDAY, PDOW)
INTEGER,INTENT(IN) :: PYEAR !current year (UTC)
INTEGER,INTENT(IN) :: PMONTH!current month (UTC)
INTEGER,INTENT(IN) :: PDAY  !current day (UTC)
INTEGER,INTENT(OUT):: PDOW  !current day of the week
END SUBROUTINE DAY_OF_WEEK
END INTERFACE
END MODULE MODI_DAY_OF_WEEK
