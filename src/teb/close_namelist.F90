!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE CLOSE_NAMELIST(HPROGRAM,KLUNAM)
!     #######################################################
!
!!****  *CLOSE_NAMELIST* - generic routine to close a namelist file
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
!!      Original    01/2004 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
#ifdef SFX_ASC
USE MODI_CLOSE_NAMELIST_ASC
#endif
#ifdef SFX_FA
USE MODI_CLOSE_NAMELIST_FA
#endif
#ifdef SFX_LFI
USE MODI_CLOSE_NAMELIST_LFI
#endif
#ifdef SFX_OL
USE MODI_CLOSE_NAMELIST_OL
#endif
#ifdef SFX_NC
USE MODI_CLOSE_NAMELIST_NC
#endif
#ifdef SFX_MNH
USE MODI_MNHCLOSE_NAMELIST
#endif
#ifdef SFX_ARO
USE MODI_AROCLOSE_NAMELIST
#endif
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! main program
INTEGER,           INTENT(IN)  :: KLUNAM   ! logical unit of namelist
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('CLOSE_NAMELIST',0,ZHOOK_HANDLE)
IF (HPROGRAM=='MESONH') THEN
#ifdef SFX_MNH
  CALL MNHCLOSE_NAMELIST(HPROGRAM,KLUNAM)
#endif
ELSE IF (HPROGRAM=='OFFLIN') THEN
#ifdef SFX_OL
  CALL CLOSE_NAMELIST_OL(HPROGRAM,KLUNAM)
#endif
ELSE IF (HPROGRAM=='ASCII ') THEN
#ifdef SFX_ASC
  CALL CLOSE_NAMELIST_ASC(HPROGRAM,KLUNAM)
#endif
ELSE IF (HPROGRAM=='AROME ') THEN
#ifdef SFX_ARO
  CALL AROCLOSE_NAMELIST(HPROGRAM,KLUNAM)
#endif
ELSE IF (HPROGRAM=='FA    ') THEN
#ifdef SFX_FA
  CALL CLOSE_NAMELIST_FA(HPROGRAM,KLUNAM)
#endif
ELSE IF (HPROGRAM=='LFI   ') THEN
#ifdef SFX_LFI
  CALL CLOSE_NAMELIST_LFI(HPROGRAM,KLUNAM)
#endif
ELSE IF (HPROGRAM=='NC    ') THEN
#ifdef SFX_NC
  CALL CLOSE_NAMELIST_NC(HPROGRAM,KLUNAM)
#endif
END IF
IF (LHOOK) CALL DR_HOOK('CLOSE_NAMELIST',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE CLOSE_NAMELIST
