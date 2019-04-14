!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ##################
      MODULE MODI_SUNPOS
!     ##################
!
INTERFACE
!
      SUBROUTINE SUNPOS (KYEAR, KMONTH, KDAY, PTIME, PLON, PLAT, PTSUN, PZENITH, PAZIMSOL)
!
INTEGER,                      INTENT(IN)   :: KYEAR      ! current year                        
INTEGER,                      INTENT(IN)   :: KMONTH     ! current month                        
INTEGER,                      INTENT(IN)   :: KDAY       ! current day                        
REAL,                         INTENT(IN)   :: PTIME      ! current time                        
REAL, DIMENSION(:),           INTENT(IN)   :: PLON       ! longitude
REAL, DIMENSION(:),           INTENT(IN)   :: PLAT       ! latitude

REAL, DIMENSION(:),           INTENT(OUT)  :: PZENITH    ! Solar zenithal angle
REAL, DIMENSION(:),           INTENT(OUT)  :: PAZIMSOL   ! Solar azimuthal angle
REAL, DIMENSION(:),           INTENT(OUT)  :: PTSUN      ! Solar time
!
END SUBROUTINE SUNPOS
!
END INTERFACE
!
END MODULE MODI_SUNPOS
!
!
!     ####################################################################################
      SUBROUTINE SUNPOS (KYEAR, KMONTH, KDAY, PTIME, PLON, PLAT, PTSUN, PZENITH, PAZIMSOL)
!     ####################################################################################
!
!!****  *SUNPOS * - routine to compute the position of the sun
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the cosine and sinus of the 
!!    solar zenithal angle (angle defined by the local vertical at the position
!!    XLAT, XLON and the direction of the sun) and the azimuthal solar
!!    angle (angle between an horizontal direction (south or north according
!!    to the terrestrial hemisphere) and the horizontal projection of the
!!    direction of the sun.
!!
!!**  METHOD
!!    ------
!!      The cosine and sinus of the zenithal solar angle  and the azimuthal 
!!    solar angle are computed from the true universal time, valid for the (XLAT,
!!    XLON) location, and from the solar declination angle of the day. There
!!    is a special convention to define the azimuthal solar angle.
!!     
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!      "Radiative Processes in Meteorology and Climatology"  
!!                          (1976)   Paltridge and Platt 
!!
!!    AUTHOR
!!    ------
!!	J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             16/10/94 
!!      Revised              12/09/95
!!      (J.Stein)            01:04/96  bug correction for ZZEANG     
!!      (K. Suhre)           14/02/97  bug correction for ZLON0     
!!      (V. Masson)          01/03/03  add zenithal angle output
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CSTS,          ONLY : XPI, XDAY
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                      INTENT(IN)   :: KYEAR      ! current year                        
INTEGER,                      INTENT(IN)   :: KMONTH     ! current month                        
INTEGER,                      INTENT(IN)   :: KDAY       ! current day                        
REAL,                         INTENT(IN)   :: PTIME      ! current time                        
REAL, DIMENSION(:),           INTENT(IN)   :: PLON       ! longitude
REAL, DIMENSION(:),           INTENT(IN)   :: PLAT       ! latutude
!
REAL, DIMENSION(:),           INTENT(OUT)  :: PZENITH    ! Solar zenithal angle
REAL, DIMENSION(:),           INTENT(OUT)  :: PAZIMSOL   ! Solar azimuthal angle
REAL, DIMENSION(:),           INTENT(OUT)  :: PTSUN      ! Solar time
!
!*       0.2   declarations of local variables
!
!
REAL                                       :: ZUT        ! Universal time
!
REAL, DIMENSION(SIZE(PLON))                :: ZTUT    ,&! True (absolute) Universal Time
                                              ZSOLANG ,&! Hourly solar angle
                                              ZSINAZI ,&! Sine of the solar azimuthal angle
                                              ZCOSAZI ,&! Cosine of the solar azimuthal angle
                                              ZLAT,    &
                                              ZLON,    &! Array of latitudes and longitudes
                                              ZSINZEN, &!Sine of zenithal angle
                                              ZCOSZEN   !Cosine of zenithal angle
INTEGER, DIMENSION(0:11)                   :: IBIS, INOBIS ! Cumulative number of days per month
                                                           ! for bissextile and regular years
REAL                                       :: ZDATE         ! Julian day of the year
REAL                                       :: ZAD           ! Angular Julian day of the year
REAL                                       :: ZDECSOL       ! Daily solar declination angle 
REAL                                       :: ZA1, ZA2      ! Ancillary variables
REAL                                       :: ZTSIDER, &
                                              ZSINDEL, &!azimuthal angle
                                              ZCOSDEL !azimuthal angle
!                                            
INTEGER                                    :: JI, JJ
!
!-------------------------------------------------------------------------------
!
!*       1.    TO COMPUTE THE TRUE SOLAR TIME
!              -------------------------------
!
ZUT  = MOD( 24.0+MOD(PTIME/3600.,24.0),24.0 )

INOBIS(:) = (/0,31,59,90,120,151,181,212,243,273,304,334/)
IBIS(0:1) = INOBIS(0:1)
DO JI=2,11
  IBIS(JI) = INOBIS(JI)+1
END DO
IF( MOD(KYEAR,4).EQ.0 ) THEN
  ZDATE = FLOAT(KDAY +   IBIS(KMONTH-1)) - 1
  ZAD = 2.0*XPI*ZDATE/366.0
ELSE
  ZDATE = FLOAT(KDAY + INOBIS(KMONTH-1)) - 1
  ZAD = 2.0*XPI*ZDATE/365.0
END IF

ZA1 = (1.00554*ZDATE- 6.28306)*(XPI/180.0)
ZA2 = (1.93946*ZDATE+23.35089)*(XPI/180.0)
ZTSIDER = (7.67825*SIN(ZA1)+10.09176*SIN(ZA2)) / 60.0
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTE THE SOLAR DECLINATION ANGLE
!	        -----------------------------------
!
ZDECSOL = 0.006918-0.399912*COS(ZAD)   +0.070257*SIN(ZAD)    &
         -0.006758*COS(2.*ZAD)+0.000907*SIN(2.*ZAD) &
         -0.002697*COS(3.*ZAD)+0.00148 *SIN(3.*ZAD)
ZSINDEL = SIN(ZDECSOL)
ZCOSDEL = COS(ZDECSOL)
!-------------------------------------------------------------------------------
!
DO JJ=1,SIZE(PLON)
!
!*       3.    LOADS THE ZLAT, ZLON ARRAYS
!              ---------------------------
!
  ZLAT(JJ) = PLAT(JJ)*(XPI/180.)
  ZLON(JJ) = PLON(JJ)*(XPI/180.)
!

!-------------------------------------------------------------------------------
!
!
!*       4.    COMPUTE THE TRUE SOLAR TIME
!              ----------------------------
!

  ZTUT(JJ) = ZUT - ZTSIDER + ZLON(JJ)*((180./XPI)/15.0)
!
  PTSUN(JJ) = MOD(PTIME -ZTSIDER*3600. +PLON(JJ)*240., XDAY)
!
!-------------------------------------------------------------------------------
!*       3.    COMPUTES THE COSINE AND SINUS OF THE ZENITHAL SOLAR ANGLE
!              ---------------------------------------------------------
!
  ZSOLANG(JJ) = (ZTUT(JJ)-12.0)*15.0*(XPI/180.)          ! hour angle in radians
!
  ZCOSZEN(JJ) = SIN(ZLAT(JJ))*ZSINDEL +                 &! Cosine of the zenithal
               COS(ZLAT(JJ))*ZCOSDEL*COS(ZSOLANG(JJ))  !       solar angle
!
  ZSINZEN(JJ)  = SQRT( 1. - ZCOSZEN(JJ)*ZCOSZEN(JJ) )
!
!-------------------------------------------------------------------------------
!
!*       5.    ZENITHAL SOLAR ANGLE
!              --------------------
!
  PZENITH(JJ) = ACOS(ZCOSZEN(JJ))
!
!-------------------------------------------------------------------------------
!
!*       6.    COMPUTE THE AZINUTHAL SOLAR ANGLE (PAZIMSOL)
!              --------------------------------------------
!
  IF (ZSINZEN(JJ)/=0.) THEN
    !Azimuth is measured counter-opposite-clockwise from due south
    ZSINAZI(JJ)  = - ZCOSDEL * SIN(ZSOLANG(JJ)) / ZSINZEN(JJ)
    ZCOSAZI(JJ)  = (-SIN(ZLAT(JJ))*ZCOSDEL*COS(ZSOLANG(JJ))     & 
                     +COS(ZLAT(JJ))*ZSINDEL                       &
                    ) / ZSINZEN(JJ)
    PAZIMSOL(JJ) = ATAN2(ZSINAZI(JJ),ZCOSAZI(JJ))
  ELSE
    PAZIMSOL(JJ) = XPI
  ENDIF
!
ENDDO
!-------------------------------------------------------------------------------
!
END SUBROUTINE SUNPOS
