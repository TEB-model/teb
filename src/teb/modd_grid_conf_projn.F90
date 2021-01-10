!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ################
      MODULE MODD_GRID_CONF_PROJ_n
!     ################
!
!!****  *MODD_GRID_CONF_PROJ_n - declaration of Arome gris characteristics
!!
!!    PURPOSE
!!    -------
!     Used if CINGRID_TYPE = 'CONF PROJ '
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
!!      V. Masson    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       20/09/02
!
!*       0.   DECLARATIONS
!             ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
REAL, DIMENSION(:), ALLOCATABLE    :: XX  ! X coordinate (meters)
REAL, DIMENSION(:), ALLOCATABLE    :: XY  ! Y coordinate (meters)
!
REAL, DIMENSION(:,:), ALLOCATABLE    :: XCX
REAL, DIMENSION(:,:), ALLOCATABLE    :: XCY
INTEGER, DIMENSION(:,:), ALLOCATABLE :: NCIJ
!
TYPE GRID_CONF_PROJ_t
!
INTEGER :: NX  ! number of points in X direction
INTEGER :: NY  ! number of points in Y direction
!
REAL    :: XLAT0  ! reference latitude
REAL    :: XLON0  ! reference longitude
REAL    :: XLATORI! origin latitude
REAL    :: XLONORI! origin longitude
REAL    :: XRPK   ! projection parameter for the conformal projection
REAL    :: XBETA  ! rotation   parameter for the conformal projection
REAL    :: XLATC ! centre latitude
REAL    :: XLONC ! centre longitude
!
END TYPE GRID_CONF_PROJ_t
!
CONTAINS
!
SUBROUTINE GRID_CONF_PROJ_INIT(YGRIDCONF)
TYPE(GRID_CONF_PROJ_t), INTENT(INOUT) :: YGRIDCONF
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_GRID_CONF_PROJ_n:GRID_CONF_PROJ_INIT",0,ZHOOK_HANDLE)
YGRIDCONF%NX = 0
YGRIDCONF%NY = 0
YGRIDCONF%XLAT0 = 0.
YGRIDCONF%XLON0 = 0.
YGRIDCONF%XLATORI = 0.
YGRIDCONF%XLONORI = 0.
YGRIDCONF%XRPK = 0.
YGRIDCONF%XBETA = 0.
YGRIDCONF%XLATC = 0.
YGRIDCONF%XLONC = 0.
IF (LHOOK) CALL DR_HOOK("MODD_GRID_CONF_PROJ_n:GRID_CONF_PROJ_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE GRID_CONF_PROJ_INIT
!
!
END MODULE MODD_GRID_CONF_PROJ_n
