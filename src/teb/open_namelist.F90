!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE OPEN_NAMELIST(HPROGRAM,KLUNAM,HFILE)
!     #######################################################
!
!!****  *OPEN_NAMELIST* - routine to open a namelist file
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
#ifdef SFX_OL
USE MODI_OPEN_NAMELIST_OL
#endif
#ifdef SFX_ASC
USE MODI_OPEN_NAMELIST_ASC
#endif
#ifdef SFX_FA
USE MODI_OPEN_NAMELIST_FA
#endif
#ifdef SFX_LFI
USE MODI_OPEN_NAMELIST_LFI
#endif
#ifdef SFX_NC
USE MODI_OPEN_NAMELIST_NC
#endif
#ifdef SFX_MNH
USE MODI_MNHOPEN_NAMELIST
#endif
#ifdef SFX_ARO
USE MODI_AROOPEN_NAMELIST
#endif
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! main program
INTEGER,           INTENT(OUT) :: KLUNAM   ! logical unit of namelist
 CHARACTER(LEN=28), INTENT(IN), OPTIONAL :: HFILE ! ASCII file to open

!
!*       0.2   Declarations of local variables
!              -------------------------------
!
 CHARACTER(LEN=28) :: YFILE
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('OPEN_NAMELIST',0,ZHOOK_HANDLE)
IF (PRESENT(HFILE)) THEN
  YFILE = HFILE
ELSE
  YFILE = '                            '
END IF

IF (HPROGRAM=='MESONH') THEN
#ifdef SFX_MNH
  CALL MNHOPEN_NAMELIST(HPROGRAM,KLUNAM,YFILE)
#endif
ELSE IF (HPROGRAM=='OFFLIN') THEN
#ifdef SFX_OL
  CALL OPEN_NAMELIST_OL(HPROGRAM,KLUNAM,YFILE)
#endif
ELSE IF (HPROGRAM=='ASCII ') THEN
#ifdef SFX_ASC
  CALL OPEN_NAMELIST_ASC(HPROGRAM,KLUNAM,YFILE)
#endif
ELSE IF (HPROGRAM=='AROME ') THEN
#ifdef SFX_ARO
  CALL AROOPEN_NAMELIST(HPROGRAM,KLUNAM,YFILE)
#endif
ELSE IF (HPROGRAM=='FA    ') THEN
#ifdef SFX_FA
  CALL OPEN_NAMELIST_FA(HPROGRAM,KLUNAM,YFILE)
#endif
ELSE IF (HPROGRAM=='LFI   ') THEN
#ifdef SFX_LFI
  CALL OPEN_NAMELIST_LFI(HPROGRAM,KLUNAM,YFILE)
#endif
ELSE IF (HPROGRAM=='NC    ') THEN
#ifdef SFX_NC
  CALL OPEN_NAMELIST_NC(HPROGRAM,KLUNAM,YFILE)
#endif
END IF
!
IF (LHOOK) CALL DR_HOOK('OPEN_NAMELIST',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE OPEN_NAMELIST
