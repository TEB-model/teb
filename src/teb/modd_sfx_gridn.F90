!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ##################
      MODULE MODD_SFX_GRID_n
!     ##################
!
!!****  *MODD_ISBA - declaration of grid for ISBA scheme
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE GRID_t
!-------------------------------------------------------------------------------
!
! Grid definition
!
  INTEGER                         :: NDIM        ! number of points
  CHARACTER(LEN=10)               :: CGRID       ! grid type
!                                              ! "NONE        " : no grid computations
!                                              ! "CONF PROJ   " : conformal projection
!                                              ! "SURF ATM    " : nature points of surf. atm. grid
!
  INTEGER                         :: NGRID_PAR   ! size of XGRID_PAR
  REAL, POINTER,     DIMENSION(:) :: XGRID_PAR   ! lits of parameters used to define the grid
!                                              ! (depends on value of CGRID)
!
!-------------------------------------------------------------------------------
!
! General surface parameters:
!
  REAL, POINTER, DIMENSION(:) :: XLAT        ! latitude (degrees +North)               (-)
  REAL, POINTER, DIMENSION(:) :: XLON        ! longitude (degrees +East)               (-)
  REAL, POINTER, DIMENSION(:) :: XMESH_SIZE  ! mesh size                               (m2)
!-------------------------------------------------------------------------------
!
END TYPE GRID_t
!
TYPE GRID_NP_t
!
TYPE(GRID_t), DIMENSION(:), POINTER :: AL=>NULL()
!
END TYPE GRID_NP_t
!
CONTAINS
!
SUBROUTINE GRID_INIT(G)
TYPE(GRID_t), INTENT(INOUT) :: G
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_SFX_GRID_N:IGRID_INIT",0,ZHOOK_HANDLE)
  NULLIFY(G%XGRID_PAR)
  NULLIFY(G%XLAT)
  NULLIFY(G%XLON)
  NULLIFY(G%XMESH_SIZE)
G%NDIM=0
G%CGRID=' '
IF (LHOOK) CALL DR_HOOK("MODD_SFX_GRID_N:GRID_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE GRID_INIT
!
SUBROUTINE GRID_NP_INIT(GP,KPATCH)
TYPE(GRID_NP_t), INTENT(INOUT) :: GP
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_SFX_GRID_N:GRID_NP_INIT",0,ZHOOK_HANDLE)

IF (ASSOCIATED(GP%AL)) THEN
  DO JP=1,KPATCH
    CALL GRID_INIT(GP%AL(JP))
  ENDDO  
  DEALLOCATE(GP%AL)
ELSE
  ALLOCATE(GP%AL(KPATCH))
  DO JP=1,KPATCH
    CALL GRID_INIT(GP%AL(JP))
  ENDDO
ENDIF  
!
IF (LHOOK) CALL DR_HOOK("MODD_SFX_GRID_N:GRID_NP_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE GRID_NP_INIT
!
END MODULE MODD_SFX_GRID_n
