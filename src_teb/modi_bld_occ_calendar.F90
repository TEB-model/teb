!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_BLD_OCC_CALENDAR
INTERFACE
SUBROUTINE BLD_OCC_CALENDAR(KYEAR, KMONTH, KDAY, PTSUN, PRESIDENTIAL, PTCOOL_TARGET_IN, PTHEAT_TARGET_IN, &
                            PQIN_IN, PDT_RES, PDT_OFF, PQIN_FRAC, PTCOOL_TARGET, PTHEAT_TARGET, PQIN)
INTEGER,INTENT(IN)                :: KYEAR         ! current year  (UTC)
INTEGER,INTENT(IN)                :: KMONTH        ! current month (UTC)
INTEGER,INTENT(IN)                :: KDAY          ! current day   (UTC)
REAL,   DIMENSION(:) , INTENT(IN) :: PTSUN         ! current solar time  (s, UTC)
REAL, DIMENSION(:) , INTENT(IN)  :: PRESIDENTIAL     ! Residential use fraction
REAL, DIMENSION(:) , INTENT(IN)  :: PTCOOL_TARGET_IN ! Cooling setpoint of HVAC system [K] input
REAL, DIMENSION(:) , INTENT(IN)  :: PTHEAT_TARGET_IN ! Heating setpoint of HVAC system [K] input
REAL, DIMENSION(:) , INTENT(IN)  :: PQIN_IN          ! Internal heat gains [W m-2(floor)] input
REAL,                INTENT(IN)  :: PDT_RES          ! Target temperature change when unoccupied (K) (residential buildings)
REAL,                INTENT(IN)  :: PDT_OFF          ! Target temperature change when unoccupied (K) (offices and commercial buildings)
REAL,                INTENT(IN)  :: PQIN_FRAC        ! Fraction of internal gains when unoccupied (-)
REAL, DIMENSION(:) , INTENT(OUT)  :: PTCOOL_TARGET ! Cooling setpoint of HVAC system [K]
REAL, DIMENSION(:) , INTENT(OUT)  :: PTHEAT_TARGET ! Heating setpoint of HVAC system [K]
REAL, DIMENSION(:) , INTENT(OUT)  :: PQIN          ! Internal heat gains [W m-2(floor)]
END SUBROUTINE BLD_OCC_CALENDAR
END INTERFACE
END MODULE MODI_BLD_OCC_CALENDAR
