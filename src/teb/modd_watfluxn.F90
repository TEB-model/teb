!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ####################
      MODULE MODD_WATFLUX_n
!     ####################
!
!!****  *MODD_WATFLUX_n - declaration of surface parameters for an inland water surface
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
!!      V. Masson  *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_DATE_SURF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE WATFLUX_t
!
! General surface: 
!
  REAL, POINTER, DIMENSION(:)   :: XZS     ! orography                     (m)
  REAL, POINTER, DIMENSION(:,:) :: XCOVER  ! fraction of each ecosystem    (-)
  LOGICAL, POINTER, DIMENSION(:):: LCOVER  ! GCOVER(i)=T --> ith cover field is not 0.
  LOGICAL                       :: LSBL    ! T: SBL scheme within the Surface Boundary Layer
!                                          ! F: no atmospheric layers below forcing level
  CHARACTER(LEN=4)              :: CWAT_ALB ! type of albedo
!
  LOGICAL                       :: LINTERPOL_TS ! Interpotalation of monthly TS
  CHARACTER(LEN=6)              :: CINTERPOL_TS ! Interpotalation of monthly TS
!
! Inland water:
!
  REAL, POINTER, DIMENSION(:) :: XTS       ! water surface temperature               (K)
  REAL, POINTER, DIMENSION(:) :: XTICE     ! water ice temperature
  REAL, POINTER, DIMENSION(:) :: XZ0       ! water surface roughness length          (-)
  REAL, POINTER, DIMENSION(:) :: XEMIS     ! water surface emissivity                (-)
  REAL, POINTER, DIMENSION(:) :: XDIR_ALB  ! water surface direct albedo             (-)
  REAL, POINTER, DIMENSION(:) :: XSCA_ALB  ! water surface diffuse albedo            (-)
  REAL, POINTER, DIMENSION(:) :: XICE_ALB  ! water ice albedo (for ESM coupling)     (-)
!
  REAL, POINTER, DIMENSION(:,:) :: XTS_MTH   ! Monthly water surface temperature               (K)
!
  REAL, POINTER, DIMENSION(:) :: XCPL_WATER_WIND ! 10m wind speed for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_WATER_FWSU ! zonal wind stress for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_WATER_FWSV ! meridian wind stress for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_WATER_SNET ! Solar net heat flux
  REAL, POINTER, DIMENSION(:) :: XCPL_WATER_HEAT ! Non solar net heat flux
  REAL, POINTER, DIMENSION(:) :: XCPL_WATER_EVAP ! Evaporation for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_WATER_RAIN ! Rainfall for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_WATER_SNOW ! Snowfall for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_WATER_FWSM ! wind stress module for ESM coupling
!
  REAL, POINTER, DIMENSION(:) :: XCPL_WATERICE_SNET ! solar net heat flux
  REAL, POINTER, DIMENSION(:) :: XCPL_WATERICE_HEAT ! Non solar net heat flux
  REAL, POINTER, DIMENSION(:) :: XCPL_WATERICE_EVAP ! Sublimation for ESM coupling
!
! Date:
!
  TYPE (DATE_TIME)                  :: TTIME         ! current date and time
  TYPE (DATE_TIME)                  :: TZTIME  
!
! Time-step:
!
  REAL                              :: XTSTEP        ! time step
!
  REAL                              :: XOUT_TSTEP    ! output writing time step
!
!
END TYPE WATFLUX_t



CONTAINS

!




SUBROUTINE WATFLUX_INIT(YWATFLUX)
TYPE(WATFLUX_t), INTENT(INOUT) :: YWATFLUX
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_WATFLUX_N:WATFLUX_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YWATFLUX%XZS)
  NULLIFY(YWATFLUX%XCOVER)
  NULLIFY(YWATFLUX%LCOVER)
  NULLIFY(YWATFLUX%XTS)
  NULLIFY(YWATFLUX%XTICE)
  NULLIFY(YWATFLUX%XZ0)
  NULLIFY(YWATFLUX%XEMIS)
  NULLIFY(YWATFLUX%XDIR_ALB)
  NULLIFY(YWATFLUX%XSCA_ALB)
  NULLIFY(YWATFLUX%XICE_ALB)
  NULLIFY(YWATFLUX%XTS_MTH)
  NULLIFY(YWATFLUX%XCPL_WATER_WIND)
  NULLIFY(YWATFLUX%XCPL_WATER_FWSU)
  NULLIFY(YWATFLUX%XCPL_WATER_FWSV)
  NULLIFY(YWATFLUX%XCPL_WATER_SNET)
  NULLIFY(YWATFLUX%XCPL_WATER_HEAT)
  NULLIFY(YWATFLUX%XCPL_WATER_EVAP)
  NULLIFY(YWATFLUX%XCPL_WATER_RAIN)
  NULLIFY(YWATFLUX%XCPL_WATER_SNOW)
  NULLIFY(YWATFLUX%XCPL_WATER_FWSM)
  NULLIFY(YWATFLUX%XCPL_WATERICE_SNET)
  NULLIFY(YWATFLUX%XCPL_WATERICE_HEAT)
  NULLIFY(YWATFLUX%XCPL_WATERICE_EVAP)
YWATFLUX%LSBL=.FALSE.
YWATFLUX%CWAT_ALB=' '
YWATFLUX%LINTERPOL_TS=.FALSE.
YWATFLUX%CINTERPOL_TS=' '
YWATFLUX%XTSTEP=0.
YWATFLUX%XOUT_TSTEP=0.
IF (LHOOK) CALL DR_HOOK("MODD_WATFLUX_N:WATFLUX_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE WATFLUX_INIT


END MODULE MODD_WATFLUX_n
