!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE TEB_IRRIG(OIRRIG, PTSTEP, KMONTH, PSOLAR_TIME,   &
               PSTART_MONTH, PEND_MONTH, PSTART_HOUR, PEND_HOUR,&
               P24H_IRRIG, PIRRIG                               )
!   ##########################################################################
!
!!****  *TEB_IRRIG*  
!!
!!    PURPOSE
!!    -------
!
! Computes the amount of water for irrigation
!         
!     
!!**  METHOD
!     ------
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!      C. DeMunck & V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/2013
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
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
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*      0.2    declarations of local variables
!
LOGICAL, DIMENSION(SIZE(PIRRIG)):: GMONTH!if irrigation occurs during present month
LOGICAL, DIMENSION(SIZE(PIRRIG)):: GHOUR !if irrigation occurs at     present time
REAL,    DIMENSION(SIZE(PIRRIG)):: ZDURATION ! duration of irrigation per day (s)
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TEB_IRRIG',0,ZHOOK_HANDLE)
!
!*      1.     Default
!              -------
!
PIRRIG = 0.
!
IF (.NOT. OIRRIG) THEN
  IF (LHOOK) CALL DR_HOOK('TEB_IRRIG',1,ZHOOK_HANDLE)
  RETURN
END IF
!
!-------------------------------------------------------------------------------
!
!*      2.     Is present month irrigated ?
!              ----------------------------
!
WHERE (PSTART_MONTH <= PEND_MONTH)   ! for summer in North hemisphere
  GMONTH = KMONTH>=PSTART_MONTH .AND. KMONTH<= PEND_MONTH
ELSEWHERE                            ! for summer in South hemisphere
                                     ! change of year during irrigation period
  GMONTH = KMONTH>=PSTART_MONTH .OR.  KMONTH<= PEND_MONTH
END WHERE
!
!-------------------------------------------------------------------------------
!
!*      3.     Is present time irrigated (for each location) ?
!              -----------------------------------------------
!
GHOUR = .FALSE.
ZDURATION = 0.
!
WHERE (PSTART_HOUR  <= PEND_HOUR)      ! typically irrigation during day
  GHOUR  = PSOLAR_TIME+PTSTEP>=PSTART_HOUR*3600. .AND. PSOLAR_TIME<=PEND_HOUR*3600.
  ZDURATION = 3600.* (PEND_HOUR-PSTART_HOUR)
ELSE WHERE                             ! typically irrigation during night
                                       ! midnight is during irrigation period
  GHOUR  = PSOLAR_TIME+PTSTEP>=PSTART_HOUR*3600. .OR.  PSOLAR_TIME< PEND_HOUR*3600.
  ZDURATION = 3600.* (PEND_HOUR + 24.-PSTART_HOUR)
END WHERE
!
!-------------------------------------------------------------------------------
!
!*      4.     Computes irrigation
!              -------------------
!
!* if duration is very short, all water is irrigated during the timestep
ZDURATION = MAX(ZDURATION, PTSTEP)
!
!* Concentration of the total 24h irrigation during the period of irrigation
WHERE(GMONTH(:) .AND. GHOUR(:)) PIRRIG(:) = P24H_IRRIG(:) / ZDURATION(:)
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TEB_IRRIG',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE TEB_IRRIG
