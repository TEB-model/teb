!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL licence version 2.1.
! See the following links for details:
! https://cecill.info/licences/Licence_CeCILL_V2.1-en.txt
! https://cecill.info/licences/Licence_CeCILL_V2.1-fr.txt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ################################################################
      SUBROUTINE OPEN_CLOSE_BIN_ASC_FORC(HACTION,HFORCING,KNI,HACTION2)
!     ################################################################
!
!!****  *OPEN_CLOSE_BIN_ASC_FORC* - routine to open and close atmospheric forcing files
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	A. Lemonsu  *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2008
!!      Modified by P. Le Moigne 07/2008: HACTION2 added
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=5),  INTENT(IN)  :: HACTION  ! action to do
CHARACTER(LEN=6),  INTENT(IN)  :: HFORCING ! forcing file type
INTEGER,           INTENT(IN)  :: KNI      ! number of points
CHARACTER(LEN=1),  INTENT(IN)  :: HACTION2 ! 'R': read, 'W': write
CHARACTER(LEN=7)               :: YSTATUS  ! file status (OLD/NEW)

!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
IF (HACTION2=='R') THEN
   YSTATUS='UNKNOWN'
ELSE IF (HACTION2=='W') THEN
   YSTATUS='NEW    '
ELSE
   CALL ABOR1_SFX('OPEN_CLOSE_BIN_ASC_FORC: UNKNOWN FILE STATUS, '//YSTATUS)      
ENDIF
!      
IF (HACTION=='CONF ') THEN
  OPEN(UNIT=261,FILE='Params_config.txt',FORM='FORMATTED',STATUS=YSTATUS)
  RETURN
END IF
!
IF (HACTION=='OPEN ') THEN
  IF (HFORCING=='ASCII ') THEN
    OPEN(UNIT=250,FILE='input/Forc_TA.txt      ',FORM='FORMATTED',STATUS=YSTATUS)
    OPEN(UNIT=251,FILE='input/Forc_QA.txt      ',FORM='FORMATTED',STATUS=YSTATUS)
    OPEN(UNIT=252,FILE='input/Forc_WIND.txt    ',FORM='FORMATTED',STATUS=YSTATUS)
    OPEN(UNIT=253,FILE='input/Forc_LW.txt      ',FORM='FORMATTED',STATUS=YSTATUS)
    OPEN(UNIT=254,FILE='input/Forc_DIR_SW.txt  ',FORM='FORMATTED',STATUS=YSTATUS)
    OPEN(UNIT=255,FILE='input/Forc_SCA_SW.txt  ',FORM='FORMATTED',STATUS=YSTATUS)
    OPEN(UNIT=256,FILE='input/Forc_RAIN.txt    ',FORM='FORMATTED',STATUS=YSTATUS)
    OPEN(UNIT=257,FILE='input/Forc_SNOW.txt    ',FORM='FORMATTED',STATUS=YSTATUS)
    OPEN(UNIT=258,FILE='input/Forc_PS.txt      ',FORM='FORMATTED',STATUS=YSTATUS)
    OPEN(UNIT=259,FILE='input/Forc_DIR.txt     ',FORM='FORMATTED',STATUS=YSTATUS)
    OPEN(UNIT=260,FILE='input/Forc_CO2.txt     ',FORM='FORMATTED',STATUS=YSTATUS)
!
  ELSE IF (HFORCING=='BINARY') THEN
    OPEN(UNIT=250,FILE='input/Forc_TA.bin      ',FORM='UNFORMATTED',STATUS=YSTATUS,ACCESS='DIRECT',RECL=KNI*4)
    OPEN(UNIT=251,FILE='input/Forc_QA.bin      ',FORM='UNFORMATTED',STATUS=YSTATUS,ACCESS='DIRECT',RECL=KNI*4)
    OPEN(UNIT=252,FILE='input/Forc_WIND.bin    ',FORM='UNFORMATTED',STATUS=YSTATUS,ACCESS='DIRECT',RECL=KNI*4)
    OPEN(UNIT=253,FILE='input/Forc_LW.bin      ',FORM='UNFORMATTED',STATUS=YSTATUS,ACCESS='DIRECT',RECL=KNI*4)
    OPEN(UNIT=254,FILE='input/Forc_DIR_SW.bin  ',FORM='UNFORMATTED',STATUS=YSTATUS,ACCESS='DIRECT',RECL=KNI*4)
    OPEN(UNIT=255,FILE='input/Forc_SCA_SW.bin  ',FORM='UNFORMATTED',STATUS=YSTATUS,ACCESS='DIRECT',RECL=KNI*4)
    OPEN(UNIT=256,FILE='input/Forc_RAIN.bin    ',FORM='UNFORMATTED',STATUS=YSTATUS,ACCESS='DIRECT',RECL=KNI*4)
    OPEN(UNIT=257,FILE='input/Forc_SNOW.bin    ',FORM='UNFORMATTED',STATUS=YSTATUS,ACCESS='DIRECT',RECL=KNI*4)
    OPEN(UNIT=258,FILE='input/Forc_PS.bin      ',FORM='UNFORMATTED',STATUS=YSTATUS,ACCESS='DIRECT',RECL=KNI*4)
    OPEN(UNIT=259,FILE='input/Forc_DIR.bin     ',FORM='UNFORMATTED',STATUS=YSTATUS,ACCESS='DIRECT',RECL=KNI*4)
    OPEN(UNIT=260,FILE='input/Forc_CO2.bin     ',FORM='UNFORMATTED',STATUS=YSTATUS,ACCESS='DIRECT',RECL=KNI*4)
  ENDIF
ENDIF

IF (HACTION=='CLOSE') THEN
  CLOSE(250)
  CLOSE(251)
  CLOSE(252)
  CLOSE(253)
  CLOSE(254)
  CLOSE(255)
  CLOSE(256)
  CLOSE(257)
  CLOSE(258)
  CLOSE(259)
  CLOSE(260)
  CLOSE(261)
!*** Modified by the coupled scheme ***
  CLOSE(262)
  CLOSE(263)
  CLOSE(264)
!
END IF
!-------------------------------------------------------------------------------
!
END SUBROUTINE OPEN_CLOSE_BIN_ASC_FORC
