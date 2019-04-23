!auto_modi:spll_teb_irrig.D
MODULE MODI_TEB_IRRIG
INTERFACE
    SUBROUTINE TEB_IRRIG(OIRRIG, PTSTEP, KMONTH, PSOLAR_TIME,   &
               PSTART_MONTH, PEND_MONTH, PSTART_HOUR, PEND_HOUR,&
               P24H_IRRIG, PIRRIG                               )
IMPLICIT NONE
LOGICAL,            INTENT(IN)  :: OIRRIG   ! Flag for irrigation or not
REAL,               INTENT(IN)  :: PTSTEP   ! time step
INTEGER,            INTENT(IN)  :: KMONTH   ! Present Month
REAL, DIMENSION(:), INTENT(IN)  :: PSOLAR_TIME ! solar time (s since midnight)
REAL, DIMENSION(:), INTENT(IN)  :: PSTART_MONTH! First month for irrigation (included)
REAL, DIMENSION(:), INTENT(IN)  :: PEND_MONTH  ! Last  month for irrigation (included)
REAL, DIMENSION(:), INTENT(IN)  :: PSTART_HOUR ! First month for irrigation (included)
REAL, DIMENSION(:), INTENT(IN)  :: PEND_HOUR   ! Last  month for irrigation (ecluded)
REAL, DIMENSION(:), INTENT(IN)  :: P24H_IRRIG  ! 24h total irrigation
REAL, DIMENSION(:), INTENT(OUT) :: PIRRIG      ! present irrigation
END SUBROUTINE TEB_IRRIG
END INTERFACE
END MODULE MODI_TEB_IRRIG
