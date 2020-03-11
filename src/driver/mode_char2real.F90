!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL licence version 2.1.
! See the following links for details:
! https://cecill.info/licences/Licence_CeCILL_V2.1-en.txt
! https://cecill.info/licences/Licence_CeCILL_V2.1-fr.txt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODE_CHAR2REAL

IMPLICIT NONE

INTERFACE ASSIGNMENT (=)
  MODULE PROCEDURE CHAR_TO_REAL
END INTERFACE


CONTAINS 

SUBROUTINE CHAR_TO_REAL(KOUT,HSTR) 
USE MODD_ARCH, ONLY : LITTLE_ENDIAN_ARCH
CHARACTER(LEN=*),DIMENSION(:),INTENT(IN)  :: HSTR
REAL,DIMENSION(SIZE(HSTR)),   INTENT(OUT) :: KOUT
INTEGER :: JI
INTEGER :: ISIZE
CHARACTER(LEN=8),DIMENSION(:),ALLOCATABLE :: YTEMP
CHARACTER(LEN=4),DIMENSION(:),ALLOCATABLE :: YTEMP2

SELECT CASE(LEN(HSTR(1)))
CASE (8)
  ISIZE = SIZE(HSTR)
  IF (LITTLE_ENDIAN_ARCH) THEN! must swap 8 bytes
    ALLOCATE(YTEMP(ISIZE))
    DO JI=1,ISIZE
      YTEMP(JI)(1:1) = HSTR(JI)(8:8)
      YTEMP(JI)(2:2) = HSTR(JI)(7:7)
      YTEMP(JI)(3:3) = HSTR(JI)(6:6)
      YTEMP(JI)(4:4) = HSTR(JI)(5:5)
      YTEMP(JI)(5:5) = HSTR(JI)(4:4)
      YTEMP(JI)(6:6) = HSTR(JI)(3:3)
      YTEMP(JI)(7:7) = HSTR(JI)(2:2)
      YTEMP(JI)(8:8) = HSTR(JI)(1:1)
    END DO
    KOUT = TRANSFER(YTEMP,1.0_8,ISIZE)
    DEALLOCATE(YTEMP)
  ELSE
    KOUT = TRANSFER(HSTR,1.0_8,ISIZE)
  END IF
CASE (4) 
  ! EMULATE a 32 bits REAL
  ISIZE = SIZE(HSTR)
  IF (LITTLE_ENDIAN_ARCH) THEN
    ALLOCATE(YTEMP2(ISIZE))
    DO JI=1,ISIZE
      YTEMP2(JI)(1:1) = HSTR(JI)(4:4)
      YTEMP2(JI)(2:2) = HSTR(JI)(3:3)
      YTEMP2(JI)(3:3) = HSTR(JI)(2:2)
      YTEMP2(JI)(4:4) = HSTR(JI)(1:1)
    END DO
    KOUT = TRANSFER(YTEMP2,1.0_4,ISIZE)
    DEALLOCATE(YTEMP2)
  ELSE
    KOUT = TRANSFER(HSTR,1.0_4,ISIZE)
  END IF 
CASE (2)
  ! EMULATE a 16 bits signed INTEGER
  IF (LITTLE_ENDIAN_ARCH) THEN! must swap 2 bytes
    KOUT = ICHAR(HSTR(:)(2:2))+256*ICHAR(HSTR(:)(1:1))
  ELSE
    KOUT = ICHAR(HSTR(:)(1:1))+256*ICHAR(HSTR(:)(2:2))
 END IF
   WHERE (KOUT > 32767) 
      KOUT = KOUT - 65536.
   END WHERE
  
CASE(1)
  ! EMULATE an 8 bits signed INTEGER
  KOUT(:) = ICHAR(HSTR(:))
  WHERE (KOUT > 127) 
    KOUT = KOUT - 256.
  END WHERE

CASE default
  CALL ABOR1_SFX('MODE_CHAR2REAL: CONVERSION ERROR IN READ_DIRECT SUBROUTINE')
END SELECT
END SUBROUTINE CHAR_TO_REAL
END MODULE MODE_CHAR2REAL
