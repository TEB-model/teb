!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!####################
MODULE MODD_SURF_PAR
!####################
!
!!****  *MODD_SURF_PAR - declaration of surface parameters
!!
!!    PURPOSE
!!    -------
!     Declaration of surface parameters
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
!!	V. Masson *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       02/2004
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!
!-----------------------------------------------------------------------------------------------------
INTEGER :: NVERSION  ! surface version
INTEGER :: NBUGFIX   ! bugfix number of this version
!
REAL,    PARAMETER :: XUNDEF = 1.E+20! undefined value
INTEGER, PARAMETER :: NUNDEF = 1E+9  ! undefined value
!-----------------------------------------------------------------------------------------------------
!
END MODULE MODD_SURF_PAR
