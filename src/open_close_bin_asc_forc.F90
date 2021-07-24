!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL version 2.1
!SFX_LIC version 1. See LICENSE, Licence_CeCILL_V2.1-en.txt and Licence_CeCILL_V2.1-fr.txt  
!SFX_LIC for details. version 1.
!     ################################################################
      SUBROUTINE OPEN_CLOSE_BIN_ASC_FORC(HACTION,HFORCING,HACTION2)
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
!!      A. Lemonsu  *Meteo France*
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
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=5),  INTENT(IN)  :: HACTION  ! action to do
 CHARACTER(LEN=6),  INTENT(IN)  :: HFORCING ! forcing file type
 CHARACTER(LEN=1),  INTENT(IN)  :: HACTION2 ! 'R': read, 'W': write
 CHARACTER(LEN=7)               :: YSTATUS  ! file status (OLD/NEW)
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: INI
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('OPEN_CLOSE_BIN_ASC_FORC',0,ZHOOK_HANDLE)
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
  IF (LHOOK) CALL DR_HOOK('OPEN_CLOSE_BIN_ASC_FORC',1,ZHOOK_HANDLE)
  RETURN
END IF
!
IF (HACTION=='OPEN ') THEN
  IF (HFORCING=='ASCII ') THEN
      OPEN(UNIT=122,FILE='input/Forc_TA.txt      ',FORM='FORMATTED',STATUS=YSTATUS)
      OPEN(UNIT=123,FILE='input/Forc_QA.txt      ',FORM='FORMATTED',STATUS=YSTATUS)
      OPEN(UNIT=124,FILE='input/Forc_WIND.txt    ',FORM='FORMATTED',STATUS=YSTATUS)
      OPEN(UNIT=125,FILE='input/Forc_LW.txt      ',FORM='FORMATTED',STATUS=YSTATUS)
      OPEN(UNIT=126,FILE='input/Forc_DIR_SW.txt  ',FORM='FORMATTED',STATUS=YSTATUS)
      OPEN(UNIT=127,FILE='input/Forc_SCA_SW.txt  ',FORM='FORMATTED',STATUS=YSTATUS)
      OPEN(UNIT=128,FILE='input/Forc_RAIN.txt    ',FORM='FORMATTED',STATUS=YSTATUS)
      OPEN(UNIT=129,FILE='input/Forc_SNOW.txt    ',FORM='FORMATTED',STATUS=YSTATUS)
      OPEN(UNIT=130,FILE='input/Forc_PS.txt      ',FORM='FORMATTED',STATUS=YSTATUS)
      OPEN(UNIT=131,FILE='input/Forc_DIR.txt     ',FORM='FORMATTED',STATUS=YSTATUS)
      OPEN(UNIT=132,FILE='input/Forc_CO2.txt     ',FORM='FORMATTED',STATUS=YSTATUS)
  ENDIF
ENDIF
!
IF (HACTION=='CLOSE') THEN
    CLOSE(121)
    CLOSE(122)
    CLOSE(123)
    CLOSE(124)
    CLOSE(125)
    CLOSE(126)
    CLOSE(127)
    CLOSE(128)
    CLOSE(129)
    CLOSE(130)
    CLOSE(131)
    CLOSE(132)
END IF
IF (LHOOK) CALL DR_HOOK('OPEN_CLOSE_BIN_ASC_FORC',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE OPEN_CLOSE_BIN_ASC_FORC
