!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!!    ###############
      MODULE MODE_POS_SURF
!!    ###############
!!
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
INTERFACE POS
!!
MODULE PROCEDURE POSNAM
MODULE PROCEDURE POSKEY
!!
END INTERFACE
!!
!!
 CONTAINS
!!
!!    ##############################################
      SUBROUTINE POSNAM(KULNAM,HDNAML,OFOUND,KLUOUT)
!!    ##############################################
!!
!!*** *POSNAM*
!!
!!    PURPOSE
!!    -------
!     To position namelist file at correct place for reading
!     namelist CDNAML.
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENT
!!    -----------------
!!
!!    REFERENCE
!!    ----------
!!       ECMWF Research Department documentation of the IFS (Hamrud)
!!
!!    AUTHOR
!!    -------
!!    I. Mallet  15/10/01 
!!
!!    MODIFICATIONS
!!    --------------
!!       I. Mallet  15/10/01     adaptation to MesoNH (F90 norm)
!------------------------------------------------------------------------------
!
IMPLICIT NONE
!
!*       0.    DECLARATIONS
!              ------------
!
!*       0.1   Declarations of arguments
!
INTEGER,          INTENT(IN) :: KULNAM
 CHARACTER(LEN=*), INTENT(IN) :: HDNAML
LOGICAL,          INTENT(OUT):: OFOUND
INTEGER, OPTIONAL,INTENT(IN) :: KLUOUT
!
!*       0.2   Declarations of local variables
!
 CHARACTER(LEN=120) :: YLINE
 CHARACTER(LEN=1)   :: YLTEST
INTEGER            :: ILEN,ILEY,INDL,IND1,IRET
INTEGER            :: J,JA, JFILE
LOGICAL            :: LLOPENED
!
 CHARACTER(LEN=1),DIMENSION(26) :: YLO=(/'a','b','c','d','e','f','g','h', &
     'i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'/)
 CHARACTER(LEN=1),DIMENSION(26) :: YUP=(/'A','B','C','D','E','F','G','H', &
     'I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       1.    POSITION FILE
!              -------------
!
IF (LHOOK) CALL DR_HOOK('MODE_POS_SURF:POSNAM',0,ZHOOK_HANDLE)
OFOUND=.FALSE.
ILEN=LEN(HDNAML)
!
!      CONTINUE READING THE FILE, THEN REWIND IF NOT FOUND
DO JFILE=1,2
  search_nam : DO
    YLINE=' '
    READ(UNIT=KULNAM,FMT='(A)',IOSTAT=IRET,END=100) YLINE

!   If file does not exist, most compilers would just create it and jump 
!   to the END label ; but a few of them would report an error:         
    IF (IRET /=0 ) THEN                 
      INQUIRE(KULNAM,OPENED=LLOPENED)
      IF (LLOPENED) THEN
        IF (PRESENT(KLUOUT)) THEN
          WRITE(KLUOUT,FMT=*) 'MODE_POS_SURF : error reading from unit ',&
                KULNAM,' file ',HDNAML,' line ',YLINE
          CALL FLUSH(KLUOUT)
        ENDIF        
        CALL ABOR1_SFX('MODE_POS_SURF: read error in namelist file') 
      ELSE
        EXIT search_nam
      END IF
    ELSE   
!     FIRST SEARCH for "&" IN THE LINE, THEN CORRECT LINE AND TEST :
      INDL=INDEX(YLINE,'&')
      IF (INDL .NE. 0 ) THEN
        ILEY=LEN(YLINE)
        DO J=1,ILEY
          DO JA=1,26
            IF (YLINE(J:J)==YLO(JA)) YLINE(J:J)=YUP(JA) 
          END DO
        END DO
        IND1=INDEX(YLINE,'&'//HDNAML)
        IF(IND1.NE.0) THEN
          YLTEST=YLINE(IND1+ILEN+1:IND1+ILEN+1)
          IF(YLTEST == ' ') THEN
!           NAMELIST FOUND : RETURN
            BACKSPACE(KULNAM)
            OFOUND=.TRUE.
            IF (PRESENT(KLUOUT)) WRITE(KLUOUT,FMT=*) '-- namelist ',HDNAML,' read'
            IF (LHOOK) CALL DR_HOOK('MODE_POS_SURF:POSNAM',1,ZHOOK_HANDLE)
            RETURN
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDDO search_nam
  100  CONTINUE
  IF(JFILE == 1) REWIND(KULNAM)
ENDDO

BACKSPACE(KULNAM)
! end of file: namelist name not found
IF (PRESENT(KLUOUT)) &
WRITE(KLUOUT,FMT=*)  &
'-- namelist ',HDNAML,' not found: default values used if required'
IF (LHOOK) CALL DR_HOOK('MODE_POS_SURF:POSNAM',1,ZHOOK_HANDLE)
!------------------------------------------------------------------
END SUBROUTINE POSNAM
!!
!!
!!    ################################################
      SUBROUTINE POSKEY(KULNAM,KLUOUT,HKEYWD1,HKEYWD2)
!!    ################################################
!!
!!*** *POSKEY*
!!
!!    PURPOSE
!!    -------
!     To position namelist file at correct place after reading
!     keyword HKEYWD
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENT
!!    -----------------
!!
!!    REFERENCE
!!    ----------
!!
!!    AUTHOR
!!    -------
!!       I. Mallet *Meteo-France*
!!
!!    MODIFICATIONS
!!    --------------
!!       Original : 15/10/01
!------------------------------------------------------------------------------
!
IMPLICIT NONE
!
!*       0.    DECLARATIONS
!              ------------
!
!*       0.1   Declarations of arguments
!
INTEGER,                    INTENT(IN) :: KULNAM
INTEGER,                    INTENT(IN) :: KLUOUT
 CHARACTER(LEN=*),           INTENT(IN) :: HKEYWD1
 CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: HKEYWD2
!
!*       0.2   Declarations of local variables
!
 CHARACTER(LEN=120) :: YLINE
INTEGER            :: ILEN1,ILEN2,IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       1.    POSITION FILE
!              -------------
!
IF (LHOOK) CALL DR_HOOK('MODE_POS_SURF:POSKEY',0,ZHOOK_HANDLE)
REWIND(KULNAM)
ILEN1=LEN(HKEYWD1)
IF (PRESENT(HKEYWD2)) ILEN2=LEN(HKEYWD2)
!
search_key : DO
      YLINE=' '
      READ(UNIT=KULNAM,FMT='(A)',IOSTAT=IRET,END=100) YLINE
      IF (IRET /=0 ) THEN
         WRITE(KLUOUT,FMT=*) '-> error when reading line from unit ',KULNAM
      ELSE
        YLINE=ADJUSTL(YLINE)
        IF (YLINE(1:ILEN1) .EQ. HKEYWD1(1:ILEN1)) EXIT search_key
      ENDIF
ENDDO search_key
!
WRITE(KLUOUT,FMT=*) '-- keyword ',HKEYWD1,' found'
!
IF (LHOOK) CALL DR_HOOK('MODE_POS_SURF:POSKEY',1,ZHOOK_HANDLE)
RETURN
!
! end of file: keyword not found
100  CONTINUE
IF (.NOT.PRESENT(HKEYWD2)) THEN
  CALL ABOR1_SFX('MODE_POS_SURF: KEYWORD NOT FOUND: '//HKEYWD1)
ELSE
!
!*       2.    SECOND KEYWORD: POSITION FILE
!              -----------------------------
!
  REWIND(KULNAM)
  search_key2 : DO
      YLINE=' '
      READ(UNIT=KULNAM,FMT='(A)',IOSTAT=IRET,END=101) YLINE
      IF (IRET /=0 ) THEN
        WRITE(KLUOUT,FMT=*) '-> error when reading line from unit ',KULNAM
      ELSE
        YLINE=ADJUSTL(YLINE)
        IF (YLINE(1:ILEN2) .EQ. HKEYWD2(1:ILEN2)) EXIT search_key2
      ENDIF
  ENDDO search_key2
  WRITE(KLUOUT,FMT=*) '-- keyword ',HKEYWD2,' found'
  IF (LHOOK) CALL DR_HOOK('MODE_POS_SURF:POSKEY',1,ZHOOK_HANDLE)
  RETURN
END IF
! end of file: scd keyword not found
101  CONTINUE
 CALL ABOR1_SFX('MODE_POS_SURF: KEYWORD NOT FOUND: '//HKEYWD2)
IF (LHOOK) CALL DR_HOOK('MODE_POS_SURF:POSKEY',1,ZHOOK_HANDLE)
!------------------------------------------------------------------
END SUBROUTINE POSKEY
!
END MODULE MODE_POS_SURF
