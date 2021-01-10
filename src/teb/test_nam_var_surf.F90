!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!############################
MODULE MODI_TEST_NAM_VAR_SURF
!############################
!
INTERFACE TEST_NAM_VAR_SURF
!
      SUBROUTINE TEST_NAM_VARC0_SURF(KLUOUT,HNAME,HVAR,       &
                                       HVALUE1,HVALUE2,HVALUE3, &
                                       HVALUE4,HVALUE5,HVALUE6, &
                                       HVALUE7,HVALUE8,HVALUE9, &
                                       HVALUE10  )  
!
INTEGER,          INTENT(IN)           ::KLUOUT   ! output listing logical unit
 CHARACTER(LEN=*) ,INTENT(IN)           ::HNAME    ! name of the variable to test
 CHARACTER(LEN=*) ,INTENT(IN)           ::HVAR     ! variable to test

 CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE1  ! first possible value
 CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE2  ! second possible value
 CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE3  ! third possible value
 CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE4  ! fourth possible value
 CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE5  ! fiveth possible value
 CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE6  ! sixth possible value
 CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE7  ! seventh possible value
 CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE8  ! eightth possible value
 CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE9  ! nineth possible value
 CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE10 ! tenth possible value
!
END SUBROUTINE TEST_NAM_VARC0_SURF
!
      SUBROUTINE TEST_NAM_VARL0_SURF(KLUOUT,HNAME,OVAR,OVALUE)
!
INTEGER,          INTENT(IN)           ::KLUOUT   ! output listing logical unit
 CHARACTER(LEN=*) ,INTENT(IN)           ::HNAME    ! name of the variable to test
LOGICAL          ,INTENT(IN)           ::OVAR     ! variable to test

LOGICAL          ,INTENT(IN), OPTIONAL ::OVALUE   ! possible value
!
END SUBROUTINE TEST_NAM_VARL0_SURF
!
      SUBROUTINE TEST_NAM_VARN0_SURF(KLUOUT,HNAME,KVAR,       &
                                       KVALUE1,KVALUE2,KVALUE3, &
                                       KVALUE4,KVALUE5,KVALUE6, &
                                       KVALUE7,KVALUE8,KVALUE9, &
                                       KVALUE10  )  
!
INTEGER,          INTENT(IN)           ::KLUOUT   ! output listing logical unit
 CHARACTER(LEN=*) ,INTENT(IN)           ::HNAME    ! name of the variable to test
INTEGER          ,INTENT(IN)           ::KVAR     ! variable to test

INTEGER          ,INTENT(IN), OPTIONAL ::KVALUE1  ! first possible value
INTEGER          ,INTENT(IN), OPTIONAL ::KVALUE2  ! second possible value
INTEGER          ,INTENT(IN), OPTIONAL ::KVALUE3  ! third possible value
INTEGER          ,INTENT(IN), OPTIONAL ::KVALUE4  ! fourth possible value
INTEGER          ,INTENT(IN), OPTIONAL ::KVALUE5  ! fiveth possible value
INTEGER          ,INTENT(IN), OPTIONAL ::KVALUE6  ! sixth possible value
INTEGER          ,INTENT(IN), OPTIONAL ::KVALUE7  ! seventh possible value
INTEGER          ,INTENT(IN), OPTIONAL ::KVALUE8  ! eightth possible value
INTEGER          ,INTENT(IN), OPTIONAL ::KVALUE9  ! nineth possible value
INTEGER          ,INTENT(IN), OPTIONAL ::KVALUE10 ! 10th possible value
!
END SUBROUTINE TEST_NAM_VARN0_SURF
!
      SUBROUTINE TEST_NAM_VARX0_SURF(KLUOUT,HNAME,PVAR,       &
                                     PVALUE1,PVALUE2,PVALUE3, &
                                     PVALUE4,PVALUE5,PVALUE6, &
                                     PVALUE7,PVALUE8,PVALUE9, &
                                     PVALUE10  )  
!
INTEGER,          INTENT(IN)          ::KLUOUT   ! output listing logical unit
 CHARACTER(LEN=*) ,INTENT(IN)          ::HNAME    ! name of the variable to test
REAL             ,INTENT(IN)          ::PVAR     ! variable to test

REAL             ,INTENT(IN), OPTIONAL ::PVALUE1  ! first possible value
REAL             ,INTENT(IN), OPTIONAL ::PVALUE2  ! second possible value
REAL             ,INTENT(IN), OPTIONAL ::PVALUE3  ! third possible value
REAL             ,INTENT(IN), OPTIONAL ::PVALUE4  ! fourth possible value
REAL             ,INTENT(IN), OPTIONAL ::PVALUE5  ! fiveth possible value
REAL             ,INTENT(IN), OPTIONAL ::PVALUE6  ! sixth possible value
REAL             ,INTENT(IN), OPTIONAL ::PVALUE7  ! seventh possible value
REAL             ,INTENT(IN), OPTIONAL ::PVALUE8  ! eightth possible value
REAL             ,INTENT(IN), OPTIONAL ::PVALUE9  ! nineth possible value
REAL             ,INTENT(IN), OPTIONAL ::PVALUE10 ! 10th possible value
!
END SUBROUTINE TEST_NAM_VARX0_SURF
!
END INTERFACE
!
END MODULE MODI_TEST_NAM_VAR_SURF
!
!
!     #########################################################
      SUBROUTINE TEST_NAM_VARC0_SURF(KLUOUT,HNAME,HVAR,       &
                                       HVALUE1,HVALUE2,HVALUE3, &
                                       HVALUE4,HVALUE5,HVALUE6, &
                                       HVALUE7,HVALUE8,HVALUE9, &
                                       HVALUE10  )  
!     #########################################################
!
!!****  *TEST_NAM_VARC0* - routine to test the value of a character var.
!!
!!    PURPOSE
!!    -------
!
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_READ
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
!!      V. MASSON      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     17/04/98
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
INTEGER,          INTENT(IN)           ::KLUOUT   ! output listing logical unit
 CHARACTER(LEN=*) ,INTENT(IN)           ::HNAME    ! name of the variable to test
 CHARACTER(LEN=*) ,INTENT(IN)           ::HVAR     ! variable to test

 CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE1  ! first possible value
 CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE2  ! second possible value
 CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE3  ! third possible value
 CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE4  ! fourth possible value
 CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE5  ! fiveth possible value
 CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE6  ! sixth possible value
 CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE7  ! seventh possible value
 CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE8  ! eightth possible value
 CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE9  ! nineth possible value
 CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE10 ! 10th possible value
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*      0.2   Declarations of local variables
!
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARC0_SURF',0,ZHOOK_HANDLE)
IF ( PRESENT (HVALUE1) ) THEN
  IF ( HVAR==HVALUE1  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARC0_SURF',1,ZHOOK_HANDLE)
  IF ( HVAR==HVALUE1 ) RETURN
END IF
!
IF ( PRESENT (HVALUE2) ) THEN
  IF ( HVAR==HVALUE2  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARC0_SURF',1,ZHOOK_HANDLE)
  IF ( HVAR==HVALUE2 ) RETURN
END IF
!
IF ( PRESENT (HVALUE3) ) THEN
  IF ( HVAR==HVALUE3  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARC0_SURF',1,ZHOOK_HANDLE)
  IF ( HVAR==HVALUE3 ) RETURN
END IF
!
IF ( PRESENT (HVALUE4) ) THEN
  IF ( HVAR==HVALUE4  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARC0_SURF',1,ZHOOK_HANDLE)
  IF ( HVAR==HVALUE4 ) RETURN
END IF
!
IF ( PRESENT (HVALUE5) ) THEN
  IF ( HVAR==HVALUE5  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARC0_SURF',1,ZHOOK_HANDLE)
  IF ( HVAR==HVALUE5 ) RETURN
END IF
!
IF ( PRESENT (HVALUE6) ) THEN
  IF ( HVAR==HVALUE6  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARC0_SURF',1,ZHOOK_HANDLE)
  IF ( HVAR==HVALUE6 ) RETURN
END IF
!
IF ( PRESENT (HVALUE7) ) THEN
  IF ( HVAR==HVALUE7  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARC0_SURF',1,ZHOOK_HANDLE)
  IF ( HVAR==HVALUE7 ) RETURN
END IF
!
IF ( PRESENT (HVALUE8) ) THEN
  IF ( HVAR==HVALUE8  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARC0_SURF',1,ZHOOK_HANDLE)
  IF ( HVAR==HVALUE8 ) RETURN
END IF
!
IF ( PRESENT (HVALUE9) ) THEN
  IF ( HVAR==HVALUE9  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARC0_SURF',1,ZHOOK_HANDLE)
  IF ( HVAR==HVALUE9 ) RETURN
END IF
!
IF ( PRESENT (HVALUE10) ) THEN
  IF ( HVAR==HVALUE10  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARC0_SURF',1,ZHOOK_HANDLE)
  IF ( HVAR==HVALUE10 ) RETURN
END IF
!
!-------------------------------------------------------------------------------
!
WRITE (KLUOUT,*) ' '
WRITE (KLUOUT,*) 'FATAL ERROR:'
WRITE (KLUOUT,*) '-----------'
WRITE (KLUOUT,*) ' '
WRITE (KLUOUT,*) 'Value "',HVAR,'" is not allowed for variable ',HNAME
WRITE (KLUOUT,*) ' '
WRITE (KLUOUT,*) 'Possible values are:'
IF ( PRESENT (HVALUE1) ) WRITE (KLUOUT,*) '"',HVALUE1,'"'
IF ( PRESENT (HVALUE2) ) WRITE (KLUOUT,*) '"',HVALUE2,'"'
IF ( PRESENT (HVALUE3) ) WRITE (KLUOUT,*) '"',HVALUE3,'"'
IF ( PRESENT (HVALUE4) ) WRITE (KLUOUT,*) '"',HVALUE4,'"'
IF ( PRESENT (HVALUE5) ) WRITE (KLUOUT,*) '"',HVALUE5,'"'
IF ( PRESENT (HVALUE6) ) WRITE (KLUOUT,*) '"',HVALUE6,'"'
IF ( PRESENT (HVALUE7) ) WRITE (KLUOUT,*) '"',HVALUE7,'"'
IF ( PRESENT (HVALUE8) ) WRITE (KLUOUT,*) '"',HVALUE8,'"'
IF ( PRESENT (HVALUE9) ) WRITE (KLUOUT,*) '"',HVALUE9,'"'
IF ( PRESENT (HVALUE10)) WRITE (KLUOUT,*) '"',HVALUE10,'"'
!
 CALL ABOR1_SFX('TEST_NAM_VAR_SURF: (1) CHARACTER VALUE NOT ALLOWED')
IF (LHOOK) CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARC0_SURF',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE TEST_NAM_VARC0_SURF
!     #########################################################
      SUBROUTINE TEST_NAM_VARL0_SURF(KLUOUT,HNAME,OVAR,OVALUE)
!     #########################################################
!
!!****  *TEST_NAM_VARL0* - routine to test the value of a logical
!!
!!    PURPOSE
!!    -------
!
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_READ
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
!!      P. Le Moigne      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     04/07
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
INTEGER,          INTENT(IN)           ::KLUOUT   ! output listing logical unit
 CHARACTER(LEN=*) ,INTENT(IN)           ::HNAME    ! name of the variable to test
LOGICAL          ,INTENT(IN)           ::OVAR     ! variable to test

LOGICAL          ,INTENT(IN), OPTIONAL ::OVALUE   ! possible value
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*      0.2   Declarations of local variables
!
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARL0_SURF',0,ZHOOK_HANDLE)
IF ( PRESENT (OVALUE) ) THEN
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARL0_SURF',1,ZHOOK_HANDLE)
  IF ( ( OVAR .AND. OVALUE ) .OR. ( .NOT.OVAR .AND. .NOT.OVALUE)  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARL0_SURF',1,ZHOOK_HANDLE)
  IF ( ( OVAR .AND. OVALUE ) .OR. ( .NOT.OVAR .AND. .NOT.OVALUE) ) RETURN
END IF
!
!-------------------------------------------------------------------------------
!
WRITE (KLUOUT,*) ' '
WRITE (KLUOUT,*) 'FATAL ERROR:'
WRITE (KLUOUT,*) '-----------'
WRITE (KLUOUT,*) ' '
WRITE (KLUOUT,*) 'Value "',OVAR,'" is not allowed for variable ',HNAME
WRITE (KLUOUT,*) ' '
WRITE (KLUOUT,*) 'Possible values are:'
IF ( PRESENT (OVALUE) ) WRITE (KLUOUT,*) '"',OVALUE,'"'
!
 CALL ABOR1_SFX('TEST_NAM_VAR_SURF: (2) LOGICAL VALUE NOT ALLOWED')
IF (LHOOK) CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARL0_SURF',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE TEST_NAM_VARL0_SURF
!     #########################################################
      SUBROUTINE TEST_NAM_VARN0_SURF(KLUOUT,HNAME,KVAR,       &
                                       KVALUE1,KVALUE2,KVALUE3, &
                                       KVALUE4,KVALUE5,KVALUE6, &
                                       KVALUE7,KVALUE8,KVALUE9, &
                                       KVALUE10  )  
!     #########################################################
!
!!****  *TEST_NAM_VARN0* - routine to test the value of an integer var.
!!
!!    PURPOSE
!!    -------
!
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_READ
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
!!      P. Le Moigne      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original  04/2007
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
INTEGER,          INTENT(IN)           ::KLUOUT   ! output listing logical unit
 CHARACTER(LEN=*) ,INTENT(IN)           ::HNAME    ! name of the variable to test
INTEGER          ,INTENT(IN)           ::KVAR     ! variable to test

INTEGER          ,INTENT(IN), OPTIONAL ::KVALUE1  ! first possible value
INTEGER          ,INTENT(IN), OPTIONAL ::KVALUE2  ! second possible value
INTEGER          ,INTENT(IN), OPTIONAL ::KVALUE3  ! third possible value
INTEGER          ,INTENT(IN), OPTIONAL ::KVALUE4  ! fourth possible value
INTEGER          ,INTENT(IN), OPTIONAL ::KVALUE5  ! fiveth possible value
INTEGER          ,INTENT(IN), OPTIONAL ::KVALUE6  ! sixth possible value
INTEGER          ,INTENT(IN), OPTIONAL ::KVALUE7  ! seventh possible value
INTEGER          ,INTENT(IN), OPTIONAL ::KVALUE8  ! eightth possible value
INTEGER          ,INTENT(IN), OPTIONAL ::KVALUE9  ! nineth possible value
INTEGER          ,INTENT(IN), OPTIONAL ::KVALUE10 ! 10th possible value
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*      0.2   Declarations of local variables
!
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',0,ZHOOK_HANDLE)
IF ( PRESENT (KVALUE1) ) THEN
  IF ( KVAR==KVALUE1  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
  IF ( KVAR==KVALUE1 ) RETURN
END IF
!
IF ( PRESENT (KVALUE2) ) THEN
  IF ( KVAR==KVALUE2  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
  IF ( KVAR==KVALUE2 ) RETURN
END IF
!
IF ( PRESENT (KVALUE3) ) THEN
  IF ( KVAR==KVALUE3  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
  IF ( KVAR==KVALUE3 ) RETURN
END IF
!
IF ( PRESENT (KVALUE4) ) THEN
  IF ( KVAR==KVALUE4  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
  IF ( KVAR==KVALUE4 ) RETURN
END IF
!
IF ( PRESENT (KVALUE5) ) THEN
  IF ( KVAR==KVALUE5  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
  IF ( KVAR==KVALUE5 ) RETURN
END IF
!
IF ( PRESENT (KVALUE6) ) THEN
  IF ( KVAR==KVALUE6  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
  IF ( KVAR==KVALUE6 ) RETURN
END IF
!
IF ( PRESENT (KVALUE7) ) THEN
  IF ( KVAR==KVALUE7  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
  IF ( KVAR==KVALUE7 ) RETURN
END IF
!
IF ( PRESENT (KVALUE8) ) THEN
  IF ( KVAR==KVALUE8  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
  IF ( KVAR==KVALUE8 ) RETURN
END IF
!
IF ( PRESENT (KVALUE9) ) THEN
  IF ( KVAR==KVALUE9  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
  IF ( KVAR==KVALUE9 ) RETURN
END IF
!
IF ( PRESENT (KVALUE10) ) THEN
  IF ( KVAR==KVALUE10  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
  IF ( KVAR==KVALUE10 ) RETURN
END IF
!
!-------------------------------------------------------------------------------
!
WRITE (KLUOUT,*) ' '
WRITE (KLUOUT,*) 'FATAL ERROR:'
WRITE (KLUOUT,*) '-----------'
WRITE (KLUOUT,*) ' '
WRITE (KLUOUT,*) 'Value "',KVAR,'" is not allowed for variable ',HNAME
WRITE (KLUOUT,*) ' '
WRITE (KLUOUT,*) 'Possible values are:'
IF ( PRESENT (KVALUE1) ) WRITE (KLUOUT,*) '"',KVALUE1,'"'
IF ( PRESENT (KVALUE2) ) WRITE (KLUOUT,*) '"',KVALUE2,'"'
IF ( PRESENT (KVALUE3) ) WRITE (KLUOUT,*) '"',KVALUE3,'"'
IF ( PRESENT (KVALUE4) ) WRITE (KLUOUT,*) '"',KVALUE4,'"'
IF ( PRESENT (KVALUE5) ) WRITE (KLUOUT,*) '"',KVALUE5,'"'
IF ( PRESENT (KVALUE6) ) WRITE (KLUOUT,*) '"',KVALUE6,'"'
IF ( PRESENT (KVALUE7) ) WRITE (KLUOUT,*) '"',KVALUE7,'"'
IF ( PRESENT (KVALUE8) ) WRITE (KLUOUT,*) '"',KVALUE8,'"'
IF ( PRESENT (KVALUE9) ) WRITE (KLUOUT,*) '"',KVALUE9,'"'
IF ( PRESENT (KVALUE10)) WRITE (KLUOUT,*) '"',KVALUE10,'"'
!
 CALL ABOR1_SFX('TEST_NAM_VAR_SURF: (3) INTEGER VALUE NOT ALLOWED')
IF (LHOOK) CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE TEST_NAM_VARN0_SURF

!     #########################################################
      SUBROUTINE TEST_NAM_VARX0_SURF(KLUOUT,HNAME,PVAR,       &
                                     PVALUE1,PVALUE2,PVALUE3, &
                                     PVALUE4,PVALUE5,PVALUE6, &
                                     PVALUE7,PVALUE8,PVALUE9, & 
                                     PVALUE10  )  
!     #########################################################
!
!!****  *TEST_NAM_VARN0* - routine to test the value of an integer var.
!!
!!    PURPOSE
!!    -------
!
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_READ
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
!!      B. Decharme      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original  01/2014
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
INTEGER,          INTENT(IN)           ::KLUOUT   ! output listing logical unit
 CHARACTER(LEN=*) ,INTENT(IN)           ::HNAME    ! name of the variable to test
REAL             ,INTENT(IN)           ::PVAR     ! variable to test

REAL             ,INTENT(IN), OPTIONAL ::PVALUE1  ! first possible value
REAL             ,INTENT(IN), OPTIONAL ::PVALUE2  ! second possible value
REAL             ,INTENT(IN), OPTIONAL ::PVALUE3  ! third possible value
REAL             ,INTENT(IN), OPTIONAL ::PVALUE4  ! fourth possible value
REAL             ,INTENT(IN), OPTIONAL ::PVALUE5  ! fiveth possible value
REAL             ,INTENT(IN), OPTIONAL ::PVALUE6  ! sixth possible value
REAL             ,INTENT(IN), OPTIONAL ::PVALUE7  ! seventh possible value
REAL             ,INTENT(IN), OPTIONAL ::PVALUE8  ! eightth possible value
REAL             ,INTENT(IN), OPTIONAL ::PVALUE9  ! nineth possible value
REAL             ,INTENT(IN), OPTIONAL ::PVALUE10 ! 10th possible value
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*      0.2   Declarations of local variables
!
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',0,ZHOOK_HANDLE)
IF ( PRESENT (PVALUE1) ) THEN
  IF ( PVAR==PVALUE1  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
  IF ( PVAR==PVALUE1 ) RETURN
END IF
!
IF ( PRESENT (PVALUE2) ) THEN
  IF ( PVAR==PVALUE2  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
  IF ( PVAR==PVALUE2 ) RETURN
END IF
!
IF ( PRESENT (PVALUE3) ) THEN
  IF ( PVAR==PVALUE3  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
  IF ( PVAR==PVALUE3 ) RETURN
END IF
!
IF ( PRESENT (PVALUE4) ) THEN
  IF ( PVAR==PVALUE4  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
  IF ( PVAR==PVALUE4 ) RETURN
END IF
!
IF ( PRESENT (PVALUE5) ) THEN
  IF ( PVAR==PVALUE5  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
  IF ( PVAR==PVALUE5 ) RETURN
END IF
!
IF ( PRESENT (PVALUE6) ) THEN
  IF ( PVAR==PVALUE6  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
  IF ( PVAR==PVALUE6 ) RETURN
END IF
!
IF ( PRESENT (PVALUE7) ) THEN
  IF ( PVAR==PVALUE7  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
  IF ( PVAR==PVALUE7 ) RETURN
END IF
!
IF ( PRESENT (PVALUE8) ) THEN
  IF ( PVAR==PVALUE8  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
  IF ( PVAR==PVALUE8 ) RETURN
END IF
!
IF ( PRESENT (PVALUE9) ) THEN
  IF ( PVAR==PVALUE9  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
  IF ( PVAR==PVALUE9 ) RETURN
END IF
!
IF ( PRESENT (PVALUE10) ) THEN
  IF ( PVAR==PVALUE10  .AND. LHOOK) &
        CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
  IF ( PVAR==PVALUE10 ) RETURN
END IF
!
!-------------------------------------------------------------------------------
!
WRITE (KLUOUT,*) ' '
WRITE (KLUOUT,*) 'FATAL ERROR:'
WRITE (KLUOUT,*) '-----------'
WRITE (KLUOUT,*) ' '
WRITE (KLUOUT,*) 'Value "',PVAR,'" is not allowed for variable ',HNAME
WRITE (KLUOUT,*) ' '
WRITE (KLUOUT,*) 'Possible values are:'
IF ( PRESENT (PVALUE1) ) WRITE (KLUOUT,*) '"',PVALUE1,'"'
IF ( PRESENT (PVALUE2) ) WRITE (KLUOUT,*) '"',PVALUE2,'"'
IF ( PRESENT (PVALUE3) ) WRITE (KLUOUT,*) '"',PVALUE3,'"'
IF ( PRESENT (PVALUE4) ) WRITE (KLUOUT,*) '"',PVALUE4,'"'
IF ( PRESENT (PVALUE5) ) WRITE (KLUOUT,*) '"',PVALUE5,'"'
IF ( PRESENT (PVALUE6) ) WRITE (KLUOUT,*) '"',PVALUE6,'"'
IF ( PRESENT (PVALUE7) ) WRITE (KLUOUT,*) '"',PVALUE7,'"'
IF ( PRESENT (PVALUE8) ) WRITE (KLUOUT,*) '"',PVALUE8,'"'
IF ( PRESENT (PVALUE9) ) WRITE (KLUOUT,*) '"',PVALUE9,'"'
IF ( PRESENT (PVALUE10)) WRITE (KLUOUT,*) '"',PVALUE10,'"'
!
 CALL ABOR1_SFX('TEST_NAM_VAR_SURF: (4) REAL VALUE NOT ALLOWED')
IF (LHOOK) CALL DR_HOOK('MODI_TEST_NAM_VAR_SURF:TEST_NAM_VARN0_SURF',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE TEST_NAM_VARX0_SURF
