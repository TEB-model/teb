!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL licence version 2.1.
! See the following links for details:
! https://cecill.info/licences/Licence_CeCILL_V2.1-en.txt
! https://cecill.info/licences/Licence_CeCILL_V2.1-fr.txt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     #####################
      MODULE MODD_SURF_CONF
!     #####################
!
!!****  *MODD_SURF_CONF - surfex configuration
!!
!!    PURPOSE
!!    -------
!     Declaration of program name
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	P. Le Moigne *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       06/2008
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!
!-----------------------------------------------------------------------------------------------------
 CHARACTER(LEN=6) :: CPROGNAME
 CHARACTER(LEN=7) :: CSOFTWARE="       " ! software used: 'PGD    ','PREP  ','OFFLINE','      '
!
END MODULE MODD_SURF_CONF
