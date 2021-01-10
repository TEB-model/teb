!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ####################
MODULE MODD_FLAKE_n
!     ####################
!
!!****  *MODD_FLAKE_n - declaration of surface parameters for the FLake model 
!!                      for inland water surfaces
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
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
TYPE FLAKE_t 
!
!-------------------------------------------------------------------------------------
! General surface: 
!-------------------------------------------------------------------------------------
!
  REAL, POINTER, DIMENSION(:) :: XZS       ! orography                     (m)
  REAL, POINTER, DIMENSION(:) :: XZ0       ! roughness length              (m)
  REAL, POINTER, DIMENSION(:) :: XUSTAR    ! air friction velocity         (m/s)
  REAL, POINTER, DIMENSION(:) :: XEMIS     ! water surface emissivity (NOT USED BY FLAKE)
!
  REAL, POINTER, DIMENSION(:,:) :: XCOVER  ! fraction of each ecosystem    (-)
!                                          ! F: no atmospheric layers below forcing level  
!
  LOGICAL, POINTER, DIMENSION(:) :: LCOVER ! GCOVER(i)=T --> ith cover field is not 0.
  LOGICAL                        :: LSBL   ! T: SBL scheme within the Surface Boundary Layer
!
!-------------------------------------------------------------------------------------
! Date and time:
!-------------------------------------------------------------------------------------
!
  TYPE (DATE_TIME)                  :: TTIME         ! current date and time
!
  REAL                              :: XTSTEP        ! time step
!
  REAL                              :: XOUT_TSTEP    ! output writing time step
!
!-------------------------------------------------------------------------------------
! FLake switches
!-------------------------------------------------------------------------------------
!
  LOGICAL            :: LSEDIMENTS  ! flag to use or not the bottom sediments
  LOGICAL            :: LSKINTEMP   ! flag to use or not the skin temperature computation
  CHARACTER(LEN=3)   :: CSNOW_FLK   ! FLake snow scheme
  CHARACTER(LEN=5)   :: CFLK_FLUX   ! Type of flux computation
  CHARACTER(LEN=4)   :: CFLK_ALB    ! Type of albedo
!
!-------------------------------------------------------------------------------------
! FLake parameters and variables
!-------------------------------------------------------------------------------------
!
  REAL, POINTER, DIMENSION(:) :: XWATER_DEPTH  ! Lake depth (m)
  REAL, POINTER, DIMENSION(:) :: XWATER_FETCH  ! Lake fetch (m)
  REAL, POINTER, DIMENSION(:) :: XT_BS         ! Temperature at the outer edge of the thermally 
                                               !       active layer of the bottom sediments [K]
  REAL, POINTER, DIMENSION(:) :: XDEPTH_BS     ! Depth of the thermally active layer of the
                                               !       bottom sediments [m]
  REAL, POINTER, DIMENSION(:) :: XCORIO        ! The Coriolis parameter [s^{-1}]
  REAL, POINTER, DIMENSION(:) :: XDIR_ALB      ! Water surface direct albedo
  REAL, POINTER, DIMENSION(:) :: XSCA_ALB      ! Water surface diffuse albedo
  REAL, POINTER, DIMENSION(:) :: XICE_ALB      ! Ice surface albedo (for ESM coupling)
  REAL, POINTER, DIMENSION(:) :: XSNOW_ALB     ! Snow surface albedo
  REAL, POINTER, DIMENSION(:) :: XEXTCOEF_WATER ! Extinction coefficient for the water [m^{-1}]
  REAL, POINTER, DIMENSION(:) :: XEXTCOEF_ICE   ! Extinction coefficient for the ice [m^{-1}]
  REAL, POINTER, DIMENSION(:) :: XEXTCOEF_SNOW  ! Extinction coefficient for the snow [m^{-1}] 
  REAL, POINTER, DIMENSION(:) :: XT_SNOW       ! Temperature at the air-snow interface [K]    
  REAL, POINTER, DIMENSION(:) :: XT_ICE        ! Temperature at the snow-ice or air-ice 
                                               !        interface [K]
  REAL, POINTER, DIMENSION(:) :: XT_MNW        ! Mean temperature of the water column [K]
  REAL, POINTER, DIMENSION(:) :: XT_WML        ! Mixed-layer temperature [K]
  REAL, POINTER, DIMENSION(:) :: XT_BOT        ! Temperature at the water-bottom sediment 
                                               !        interface [K]
  REAL, POINTER, DIMENSION(:) :: XT_B1         ! Temperature at the bottom of the upper 
                                               !        layer of the sediments [K]
  REAL, POINTER, DIMENSION(:) :: XCT           ! Shape factor (thermocline)
  REAL, POINTER, DIMENSION(:) :: XH_SNOW       ! Snow thickness [m]
  REAL, POINTER, DIMENSION(:) :: XH_ICE        ! Ice thickness [m]
  REAL, POINTER, DIMENSION(:) :: XH_ML         ! Thickness of the mixed-layer [m]
  REAL, POINTER, DIMENSION(:) :: XH_B1         ! Thickness of the upper layer of bottom sediments [m]                                    
!
  REAL, POINTER, DIMENSION(:) :: XTS  ! surface temperature  (K)
                                      ! (water or ice or snow)
!
!-------------------------------------------------------------------------------------
! Coupling field for Earth system model
!-------------------------------------------------------------------------------------
!
  REAL, POINTER, DIMENSION(:) :: XCPL_FLAKE_EVAP ! Evaporation for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_FLAKE_RAIN ! Rainfall for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_FLAKE_SNOW ! Snowfall for ESM coupling
!
END TYPE FLAKE_t
!
!-------------------------------------------------------------------------------------
!


CONTAINS
!
!-------------------------------------------------------------------------------------
!
!
!


!
!
!-------------------------------------------------------------------------------------
!
SUBROUTINE FLAKE_INIT(YFLAKE)
TYPE(FLAKE_t), INTENT(INOUT) :: YFLAKE
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_FLAKE_N:FLAKE_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YFLAKE%XZS)
  NULLIFY(YFLAKE%XZ0)
  NULLIFY(YFLAKE%XUSTAR)
  NULLIFY(YFLAKE%XCOVER)
  NULLIFY(YFLAKE%LCOVER)
  NULLIFY(YFLAKE%XEMIS)
  NULLIFY(YFLAKE%XWATER_DEPTH)
  NULLIFY(YFLAKE%XWATER_FETCH)
  NULLIFY(YFLAKE%XT_BS)
  NULLIFY(YFLAKE%XDEPTH_BS)
  NULLIFY(YFLAKE%XCORIO)
  NULLIFY(YFLAKE%XDIR_ALB)
  NULLIFY(YFLAKE%XSCA_ALB)
  NULLIFY(YFLAKE%XICE_ALB)
  NULLIFY(YFLAKE%XSNOW_ALB)
  NULLIFY(YFLAKE%XEXTCOEF_WATER)
  NULLIFY(YFLAKE%XEXTCOEF_ICE)
  NULLIFY(YFLAKE%XEXTCOEF_SNOW)
  NULLIFY(YFLAKE%XT_SNOW)
  NULLIFY(YFLAKE%XT_ICE)
  NULLIFY(YFLAKE%XT_MNW)
  NULLIFY(YFLAKE%XT_WML)
  NULLIFY(YFLAKE%XT_BOT)
  NULLIFY(YFLAKE%XT_B1)
  NULLIFY(YFLAKE%XCT)
  NULLIFY(YFLAKE%XH_SNOW)
  NULLIFY(YFLAKE%XH_ICE)
  NULLIFY(YFLAKE%XH_ML)
  NULLIFY(YFLAKE%XH_B1)
  NULLIFY(YFLAKE%XTS)
  NULLIFY(YFLAKE%XCPL_FLAKE_EVAP)
  NULLIFY(YFLAKE%XCPL_FLAKE_RAIN)
  NULLIFY(YFLAKE%XCPL_FLAKE_SNOW)
YFLAKE%LSBL=.FALSE.
YFLAKE%XTSTEP=0.
YFLAKE%XOUT_TSTEP=0.
YFLAKE%LSEDIMENTS=.FALSE.
YFLAKE%LSKINTEMP=.FALSE.
YFLAKE%CSNOW_FLK='   '
YFLAKE%CFLK_ALB='    '
YFLAKE%CFLK_FLUX='     '
IF (LHOOK) CALL DR_HOOK("MODD_FLAKE_N:FLAKE_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE FLAKE_INIT
!
!-------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------
!
END MODULE MODD_FLAKE_n
