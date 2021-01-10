!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ##################
      MODULE MODD_SURF_ATM_GRID_n
!     ##################
!
!!****  *MODD_SURF_ATM_GRID - declaration of SURF_ATM grid
!!
!!    PURPOSE
!!    -------
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
!!      V. Masson  *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_SFX_GRID_n, ONLY : GRID_t, GRID_INIT
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE SURF_ATM_GRID_t
!-------------------------------------------------------------------------------
!
! Grid definition
!
TYPE(GRID_t) :: G
!                                              ! (depends on value of CGRID)
  REAL, POINTER,     DIMENSION(:) :: XGRID_FULL_PAR   ! lits of parameters used to define the grid
!                                                     ! (depends on value of CGRID)
  INTEGER                         :: NGRID_FULL_PAR   ! size of XGRID_FULL_PAR
!
  INTEGER, POINTER, DIMENSION(:,:) :: NNEAR
!-------------------------------------------------------------------------------
!
! General surface parameters:
!
  REAL, POINTER, DIMENSION(:) :: XJPDIR      ! heading of J direction (deg from N clockwise)
!-------------------------------------------------------------------------------
!
END TYPE SURF_ATM_GRID_t
!
 CONTAINS
!
SUBROUTINE SURF_ATM_GRID_INIT(YSURF_ATM_GRID)
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: YSURF_ATM_GRID
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_SURF_ATM_GRID_N:SURF_ATM_GRID_INIT",0,ZHOOK_HANDLE)
  CALL GRID_INIT(YSURF_ATM_GRID%G)
  NULLIFY(YSURF_ATM_GRID%NNEAR)
  NULLIFY(YSURF_ATM_GRID%XGRID_FULL_PAR)
  NULLIFY(YSURF_ATM_GRID%XJPDIR) 
YSURF_ATM_GRID%NGRID_FULL_PAR=0
IF (LHOOK) CALL DR_HOOK("MODD_SURF_ATM_GRID_N:SURF_ATM_GRID_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE SURF_ATM_GRID_INIT
!
END MODULE MODD_SURF_ATM_GRID_n
