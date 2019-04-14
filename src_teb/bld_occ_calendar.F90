!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!#####################################################################################
SUBROUTINE BLD_OCC_CALENDAR(KYEAR, KMONTH, KDAY, PTSUN, PRESIDENTIAL, PTCOOL_TARGET_IN, PTHEAT_TARGET_IN, &
                            PQIN_IN, PDT_RES, PDT_OFF, PQIN_FRAC, PTCOOL_TARGET, PTHEAT_TARGET, PQIN)
!#####################################################################################
!! **** BLD_OCC_CALENDAR *
!!
!!  PURPOSE
!!  -------
!! BLD_OCC_CALENDAR defines the parameters relevant for BEM (TCOOL_TARGET, THEAT_TARGET and QIN)
!! depending on the building use type (hence the DOW and local French TOD) * based on MUSCADE scenarii *
!!
!! AUTHOR
!! ------
!! C. de Munck     *Météo-France*
!! 
!! MODIFICATIONS
!! -------------
!! Original  02/2013
!
USE MODD_CSTS,         ONLY : XTT
!
USE MODI_DAY_OF_WEEK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!! 1. declaration of arguments
!
INTEGER,INTENT(IN)                :: KYEAR         ! current year  (UTC)
INTEGER,INTENT(IN)                :: KMONTH        ! current month (UTC)
INTEGER,INTENT(IN)                :: KDAY          ! current day   (UTC)
REAL,   DIMENSION(:) , INTENT(IN) :: PTSUN         ! current solar time  (s, UTC)
!
REAL, DIMENSION(:) , INTENT(IN)  :: PRESIDENTIAL     ! Residential use fraction
REAL, DIMENSION(:) , INTENT(IN)  :: PTCOOL_TARGET_IN ! Cooling setpoint of HVAC system [K] input
REAL, DIMENSION(:) , INTENT(IN)  :: PTHEAT_TARGET_IN ! Heating setpoint of HVAC system [K] input
REAL, DIMENSION(:) , INTENT(IN)  :: PQIN_IN          ! Internal heat gains [W m-2(floor)] input
REAL,                INTENT(IN)  :: PDT_RES          ! Target temperature change when unoccupied (K) (residential buildings)
REAL,                INTENT(IN)  :: PDT_OFF          ! Target temperature change when unoccupied (K) (offices and commercial buildings)
REAL,                INTENT(IN)  :: PQIN_FRAC        ! Fraction of internal gains when unoccupied (-)
!
REAL, DIMENSION(:) , INTENT(OUT)  :: PTCOOL_TARGET ! Cooling setpoint of HVAC system [K]
REAL, DIMENSION(:) , INTENT(OUT)  :: PTHEAT_TARGET ! Heating setpoint of HVAC system [K]
REAL, DIMENSION(:) , INTENT(OUT)  :: PQIN          ! Internal heat gains [W m-2(floor)]
!
!! 2. declaration of local variables
!
 CHARACTER(LEN = 6)                :: CTIME         ! Local time scheme (winter or summer time)
INTEGER                           :: JDOW          ! day of week
INTEGER                           :: JJ
!
REAL, DIMENSION(SIZE(PQIN))       :: ZTOD_BEG      ! first Time Of Day that building unoccupied (UTC, s)
REAL, DIMENSION(SIZE(PQIN))       :: ZTOD_END      ! last Time Of Day that building unoccupied (UTC, s)
REAL, DIMENSION(SIZE(PQIN))       :: ZDT           ! Target temperature change when unoccupied (K)

REAL(KIND=JPRB) :: ZTODOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('BLD_OCC_CALENDAR',0,ZTODOOK_HANDLE)
!
!--------------------------------------------------------------------------------------
!  3. determine the day of the week and the local time scheme in France
!--------------------------------------------------------------------------------------
!
 CALL DAY_OF_WEEK(KYEAR, KMONTH, KDAY, JDOW)
!
 CTIME = 'WINTER'
IF (KMONTH >= 4 .AND. KMONTH <= 10) CTIME = 'SUMMER'
!
!--------------------------------------------------------------------------------------
!  4. initialisation of parameters
!--------------------------------------------------------------------------------------
!
! Parameters assigned to the occupied values - read in namelist via BATI.csv :
!
PTHEAT_TARGET(:) = PTHEAT_TARGET_IN(:)
!  
PTCOOL_TARGET(:) = PTCOOL_TARGET_IN(:)
! 
PQIN(:)          = PQIN_IN(:)
!
ZTOD_BEG(:) = 0.
ZTOD_END(:) = 0.
!
!--------------------------------------------------------------------------------------
!  5. computes beginning and end of unoccupied calendar based on building USE TYPE
!--------------------------------------------------------------------------------------
!
DO JJ =1,SIZE(PTSUN)
!
 IF (PRESIDENTIAL(JJ) > 0.5) THEN ! RESIDENTIAL
   !
   IF (JDOW >= 2 .AND. JDOW <=6) THEN ! week days
     ZTOD_BEG(JJ) =  9. * 3600.       !  9 UTC - WINTER time
     ZTOD_END(JJ) = 17. * 3600.       ! 17 UTC - WINTER time
   END IF
   ZDT(JJ) = PDT_RES
   !
 ELSE
   !     
   IF (JDOW >= 2 .AND. JDOW <=7) THEN ! week days
     ZTOD_BEG(JJ) = 17. * 3600.       ! 17 UTC
     ZTOD_END(JJ) =  7. * 3600.       !  7 UTC
   ELSE                               ! week-end
     ZTOD_BEG(JJ) =  0. * 3600.       !   0 UTC
     ZTOD_END(JJ) = 24. * 3600.       !  24 UTC
   END IF
   ZDT(JJ) = PDT_OFF
   !
 END IF
! adjustment of unoccupied TOD based on time scheme
 IF (CTIME == 'SUMMER') THEN
  ZTOD_BEG(JJ) = ZTOD_BEG(JJ) - 3600.   
  ZTOD_END(JJ) = ZTOD_END(JJ) - 3600.   
 END IF
!
ENDDO
!
!--------------------------------------------------------------------------------------
!  6. modulate BEM input values for unoccupied building calendar
!--------------------------------------------------------------------------------------
!
DO JJ =1,SIZE(PTSUN)
!
     IF (( (ZTOD_BEG(JJ) < ZTOD_END(JJ)) .AND. (PTSUN(JJ) > ZTOD_BEG(JJ) .AND. PTSUN(JJ) < ZTOD_END(JJ))    )      &
        .OR.                                                                                         &
         ( (ZTOD_BEG(JJ) > ZTOD_END(JJ)) .AND. ((PTSUN(JJ) > 0 .AND. PTSUN(JJ) < ZTOD_END(JJ)) .OR.            &
                                        (PTSUN(JJ) > ZTOD_BEG(JJ) .AND. PTSUN(JJ) < 24 * 3600.)))) THEN             
         !
            PTHEAT_TARGET(JJ) = PTHEAT_TARGET_IN(JJ) - ZDT(JJ)
            PTCOOL_TARGET(JJ) = PTCOOL_TARGET_IN(JJ) + ZDT(JJ)
            PQIN         (JJ) = PQIN_FRAC * PQIN(JJ)
      ENDIF
!
ENDDO
!
!--------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('BLD_OCC_CALENDAR',1,ZTODOOK_HANDLE)
!
END SUBROUTINE BLD_OCC_CALENDAR
