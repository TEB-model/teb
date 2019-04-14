!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     #################
      MODULE MODD_TYPE_DATE_SURF
!     #################
!
!!****  *MODD_TYPE_DATE* - declaration of temporal types
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to define
!      the time types. 
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE 
!!
!!    REFERENCE
!!    --------- 
!!      Book2 of documentation of Meso-NH (module MODD_TYPE_DATE)
!!       
!!    AUTHOR
!!    ------
!!	P. Jabouille   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    11/08/97                      
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!
TYPE DATE
INTEGER :: YEAR = 0
INTEGER :: MONTH = 0
INTEGER :: DAY = 0
END TYPE DATE
!
TYPE DATE_TIME
TYPE (DATE) :: TDATE
REAL :: TIME = 0.
END TYPE DATE_TIME 
!
END MODULE MODD_TYPE_DATE_SURF
