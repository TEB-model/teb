!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL licence version 2.1.
! See the following links for details:
! https://cecill.info/licences/Licence_CeCILL_V2.1-en.txt
! https://cecill.info/licences/Licence_CeCILL_V2.1-fr.txt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     #############################
      MODULE MODI_CLOSE_FILE_ASC
!     #############################
INTERFACE
      SUBROUTINE CLOSE_FILE_ASC(HPROGRAM,KUNIT )
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! main program
INTEGER,           INTENT(IN)  :: KUNIT    ! logical unit of file
!
END SUBROUTINE CLOSE_FILE_ASC
!
END INTERFACE
END MODULE MODI_CLOSE_FILE_ASC
!
!     #######################################################
      SUBROUTINE CLOSE_FILE_ASC(HPROGRAM,KUNIT )
!     #######################################################
!
!!****  *CLOSE_FILE_ASC* - closes files read by surface in MESOHN
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
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003 
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
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! main program
INTEGER,           INTENT(IN)  :: KUNIT    ! logical unit of file
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
!* closes the file
!  -------------------
!
CLOSE(KUNIT)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE CLOSE_FILE_ASC
