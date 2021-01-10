!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ####################
      MODULE MODD_DUMMY_SURF_FIELDS_n
!     ####################
!
!!****  *MODD_DUMMY_SURF_FIELDS* - declaration of dummy physiographic data arrays
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     dummy physiographic data arrays.
!
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    03/2004                      
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE DUMMY_SURF_FIELDS_t
!
  INTEGER                                          :: NDUMMY_NBR
!                          ! number of dummy pgd fields chosen by user
  CHARACTER(LEN=3) , DIMENSION(:), POINTER         :: CDUMMY_AREA 
!                          ! areas where dummy pgd fields are defined
!                          ! 'ALL' : everywhere
!                          ! 'SEA' : where sea exists
!                          ! 'LAN' : where land exists
!                          ! 'WAT' : where inland water exists
!                          ! 'NAT' : where natural or agricultural areas exist
!                          ! 'TWN' : where town areas exist
!                          ! 'STR' : where streets are present
!                          ! 'BLD' : where buildings are present
!                          !
  CHARACTER(LEN=20), DIMENSION(:), POINTER         :: CDUMMY_NAME 
!                          ! name of the dummy pgd fields (for information)
  REAL,              DIMENSION(:,:), POINTER       :: XDUMMY_FIELDS
!                          ! dummy pgd fields themselves
!
!-------------------------------------------------------------------------------
!
END TYPE DUMMY_SURF_FIELDS_t

CONTAINS
!
!
SUBROUTINE DUMMY_SURF_FIELDS_INIT(YDUMMY_SURF_FIELDS)
TYPE(DUMMY_SURF_FIELDS_t), INTENT(INOUT) :: YDUMMY_SURF_FIELDS
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DUMMY_SURF_FIELDS_N:DUMMY_SURF_FIELDS_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YDUMMY_SURF_FIELDS%CDUMMY_NAME)
  NULLIFY(YDUMMY_SURF_FIELDS%CDUMMY_AREA)
  NULLIFY(YDUMMY_SURF_FIELDS%XDUMMY_FIELDS)
YDUMMY_SURF_FIELDS%NDUMMY_NBR=0
IF (LHOOK) CALL DR_HOOK("MODD_DUMMY_SURF_FIELDS_N:DUMMY_SURF_FIELDS_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DUMMY_SURF_FIELDS_INIT


END MODULE MODD_DUMMY_SURF_FIELDS_n
