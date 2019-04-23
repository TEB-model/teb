!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!################################################
SUBROUTINE DAY_OF_WEEK(PYEAR, PMONTH, PDAY, PDOW)
!################################################
!! DAY_OF_WEEK
!!
!!  PURPOSE
!!  -------
!! DAY_OF_WEEK COMPUTES THE DAY OF THE WEEK BASED UPON THE GIVEN DATE,
!! MONTH AND YEAR.  IT USES THE ZELLER CONGRUENCE ALGORITHIM.
!! PDAY IS THE DAY OF THE MONTH, 1 - 31
!! PMONTH IS THE MONTH OF THE YEAR, 1 - 12
!! PYEAR IS THE YEAR, E.G., 1977
!! IT RETURNS 1 FOR SUNDAY, 2 FOR MONDAY, ETC.
!!
!! AUTHOR
!! ------
!! G. Pigeon *Météo-France*
!! 
!! MODIFICATIONS
!! -------------
!! Original  02/2010

USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
!! DECLARATION
IMPLICIT NONE
!! 1. declaraction of arguments
INTEGER,INTENT(IN) :: PYEAR !current year (UTC)
INTEGER,INTENT(IN) :: PMONTH!current month (UTC)
INTEGER,INTENT(IN) :: PDAY  !current day (UTC)
INTEGER,INTENT(OUT):: PDOW  !current day of the week
!!
!! 2. declaration of local variables
INTEGER :: DAY, YR, MN, N1, N2
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DAY_OF_WEEK',0,ZHOOK_HANDLE)
!
YR = PYEAR
MN = PMONTH
   
!   IF JANUARY OR FEBRUARY, ADJUST MONTH AND YEAR
   
IF (MN .LE. 2) THEN
   MN = MN + 12
   YR = YR - 1
END IF
N1 = (26 * (MN + 1)) / 10
N2 = (125 * YR) / 100
DAY = PDAY + N1 + N2 - (YR / 100) + (YR / 400) - 1
PDOW = MOD(DAY, 7) + 1
!
IF (LHOOK) CALL DR_HOOK('DAY_OF_WEEK',1,ZHOOK_HANDLE)
!
END SUBROUTINE DAY_OF_WEEK
