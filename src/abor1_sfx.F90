!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL licence version 2.1.
! See the following links for details:
! https://cecill.info/licences/Licence_CeCILL_V2.1-en.txt
! https://cecill.info/licences/Licence_CeCILL_V2.1-fr.txt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     #############################################################
      SUBROUTINE ABOR1_SFX(YTEXT)
!     #############################################################
!
!!****  *ABOR1_SFX* - abor1 subroutine
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
!!	P. Le Moigne   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/2008 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODI_CLOSE_FILE
USE MODD_SURF_CONF, ONLY : CPROGNAME
!      
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=*),  INTENT(IN)  :: YTEXT
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
CHARACTER(LEN=6)  :: YPROGRAM      
INTEGER           :: ILUOUT = 25        ! logical unit of output file      
!
!-------------------------------------------------------------------------------
!
!* get output listing file logical unit
!
YPROGRAM = CPROGNAME
!      
!
IF (YPROGRAM=='ASCII ' .OR. YPROGRAM=='TEXTE ' .OR. YPROGRAM=='BINARY') THEN
   WRITE(*,*)YTEXT
   WRITE(*,*)'---------------------------------------------------------------------------'
   WRITE(*,*) 'MORE DETAILS ABOUT THE CRASH IN THE OUTPUT LISTING: SEE THE FILE NAMED    '
   WRITE(*,*) 'LISTING_[NAME OF THE RUNNING .EXE: PGD, PREP, OFFLINE].txt              '
   WRITE(*,*)'---------------------------------------------------------------------------'
ENDIF
!
WRITE(ILUOUT,*) '---------------------------------------------------------------------------'
WRITE(ILUOUT,*) '---------------------------------------------------------------------------'
WRITE(ILUOUT,*) '--------------------   FATAL ERROR in SURFEX  -----------------------------'
WRITE(ILUOUT,*) '---------------------------------------------------------------------------'
WRITE(ILUOUT,*) '---------------------------------------------------------------------------'
WRITE(ILUOUT,*) '-                                                                         -'
WRITE(ILUOUT,*)YTEXT
WRITE(ILUOUT,*) '-                                                                         -'
WRITE(ILUOUT,*) '---------------------------------------------------------------------------'
WRITE(ILUOUT,*) '---------------------------------------------------------------------------'
CALL CLOSE_FILE(YPROGRAM,ILUOUT)
!
STOP
!
END SUBROUTINE ABOR1_SFX
