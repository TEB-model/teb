!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
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
