!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ################
      MODULE MODD_BLD_DESCRIPTION_n
!     ################
!
!!****  *MODD_BLD_DESCRIPTION_n - declaration of surface parameters for typical
!                               buildings
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
!!      G. Pigeon   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       08/2011
!!       V. Masson     08/2013 adds solar panels
!!       V. Masson     10/2013 adds residential fraction
!!----------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE BLD_DESC_t
!
! Number of layers
!
  INTEGER                       :: NDESC_BLD          ! number of types of buildings
  INTEGER                       :: NDESC_AGE          ! number of building's construction dates ranges
  INTEGER                       :: NDESC_CODE         ! number of codes for buildings (merges type & age)
  INTEGER                       :: NDESC_USE          ! number of types of building's uses
  INTEGER                       :: NDESC_ROOF_LAYER   ! number of layers in roofs
  INTEGER                       :: NDESC_ROAD_LAYER   ! number of layers in roads
  INTEGER                       :: NDESC_WALL_LAYER   ! number of layers in walls
  INTEGER                       :: NDESC_FLOOR_LAYER  ! number of layers in floor
  INTEGER, POINTER, DIMENSION(:):: NDESC_BLD_LIST     ! list of the types for buildings
  INTEGER, POINTER, DIMENSION(:):: NDESC_AGE_DATE     ! list of the contruction dates for buildings
  INTEGER, POINTER, DIMENSION(:):: NDESC_CODE_LIST    ! list of the codes for buildings
  INTEGER, POINTER, DIMENSION(:):: NDESC_AGE_LIST     ! list of the contruction dates' codes
  INTEGER, POINTER, DIMENSION(:):: NDESC_USE_LIST     ! list of the codes for building's uses
  !
  REAL, POINTER, DIMENSION(:)   :: XDESC_ALB_ROOF     ! Roof albedo
  REAL, POINTER, DIMENSION(:)   :: XDESC_ALB_ROAD     ! Road albedo
  REAL, POINTER, DIMENSION(:)   :: XDESC_ALB_WALL     ! Wall albedo
  REAL, POINTER, DIMENSION(:)   :: XDESC_EMIS_ROOF    ! Roof emissivity
  REAL, POINTER, DIMENSION(:)   :: XDESC_EMIS_ROAD    ! Road emissivity
  REAL, POINTER, DIMENSION(:)   :: XDESC_EMIS_WALL    ! Wall emissivity
  REAL, POINTER, DIMENSION(:)   :: XDESC_TCOOL_TARGET ! cooling setpoint of indoor air
  REAL, POINTER, DIMENSION(:)   :: XDESC_THEAT_TARGET ! heating setpoint of indoor air  
  REAL, POINTER, DIMENSION(:)   :: XDESC_F_WASTE_CAN  ! fraction of waste heat into the canyon
  REAL, POINTER, DIMENSION(:)   :: XDESC_EFF_HEAT     ! efficiency of the heating system
  REAL, POINTER, DIMENSION(:,:) :: XDESC_HC_FLOOR     ! heat capacity of floor layers [J m-3 K-1]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_TC_FLOOR     ! thermal conductivity of floor layers [W m-1 K-1]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_D_FLOOR      ! thickness of floor layers [m]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_HC_ROOF      ! heat capacity of roof layers [J m-3 K-1]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_TC_ROOF      ! thermal conductivity of roof layers [W m-1 K-1]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_D_ROOF       ! thickness of roof layers [m]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_HC_ROAD      ! heat capacity of road layers [J m-3 K-1]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_TC_ROAD      ! thermal conductivity of road layers [W m-1 K-1]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_D_ROAD       ! thickness of road layers [m]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_HC_WALL      ! heat capacity of wall layers [J m-3 K-1]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_TC_WALL      ! thermal conductivity of wall layers [W m-1 K-1]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_D_WALL       ! thickness of wall layers [m]
  REAL, POINTER, DIMENSION(:)   :: XDESC_QIN          ! internal heat gains [W m-2(floor)]
  REAL, POINTER, DIMENSION(:)   :: XDESC_QIN_FRAD     ! radiant fraction of internal heat gains
  REAL, POINTER, DIMENSION(:)   :: XDESC_SHGC         ! solar transmitance of windows
  REAL, POINTER, DIMENSION(:)   :: XDESC_U_WIN        ! glazing thermal resistance [K m W-2]
  REAL, POINTER, DIMENSION(:)   :: XDESC_GR           ! glazing ratio
  REAL, POINTER, DIMENSION(:)   :: XDESC_SHGC_SH      ! solar transmitance of windows + shading
  REAL, POINTER, DIMENSION(:)   :: XDESC_FLOOR_HEIGHT ! building floor height [m]
  REAL, POINTER, DIMENSION(:)   :: XDESC_INF          ! infiltration/ventilation flow rate [AC/H]
!
  REAL, POINTER, DIMENSION(:)   :: XDESC_F_WATER_COND ! fraction of evaporation for condensers
  REAL, POINTER, DIMENSION(:)   :: XDESC_SHADE        ! Flag to activate shading devices 0->inactivated , 1->activated
  REAL, POINTER, DIMENSION(:)   :: XDESC_NATVENT      ! Flag to describe bld surventilation solution 0-> NONE ; 1 -> MANU ; 2-> AUTO
  REAL, POINTER, DIMENSION(:)   :: XDESC_QIN_FLAT     ! Latent franction of internal heat gains
  REAL, POINTER, DIMENSION(:)   :: XDESC_HR_TARGET    ! Relative humidity setpoint  
  REAL, POINTER, DIMENSION(:)   :: XDESC_V_VENT       ! Ventilation flow rate [AC/H]
  REAL, POINTER, DIMENSION(:)   :: XDESC_COP_RAT      ! Rated COP of the cooling system
  REAL, POINTER, DIMENSION(:)   :: XDESC_GREENROOF    ! Greenroof fraction
  REAL, POINTER, DIMENSION(:)   :: XDESC_EMIS_PANEL   ! Emissivity of Solar panels
  REAL, POINTER, DIMENSION(:)   :: XDESC_ALB_PANEL    ! Albedo     of Solar panels
  REAL, POINTER, DIMENSION(:)   :: XDESC_EFF_PANEL    ! Efficiency of Solar panels
  REAL, POINTER, DIMENSION(:)   :: XDESC_FRAC_PANEL   ! Fraction   of Solar panels on roofs
  REAL, POINTER, DIMENSION(:)   :: XDESC_RESIDENTIAL  ! Fraction of residential use
!
END TYPE BLD_DESC_t



CONTAINS

!
!


!
!
SUBROUTINE BLD_DESC_INIT(YBLD_DESC)
TYPE(BLD_DESC_t), INTENT(INOUT) :: YBLD_DESC
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_BLD_DESCRIPTION_n:BLD_DESC_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YBLD_DESC%NDESC_BLD_LIST)
  NULLIFY(YBLD_DESC%NDESC_CODE_LIST)
  NULLIFY(YBLD_DESC%NDESC_AGE_LIST)
  NULLIFY(YBLD_DESC%NDESC_AGE_DATE)
  NULLIFY(YBLD_DESC%NDESC_USE_LIST)
  NULLIFY(YBLD_DESC%XDESC_ALB_ROOF)
  NULLIFY(YBLD_DESC%XDESC_ALB_ROAD)
  NULLIFY(YBLD_DESC%XDESC_ALB_WALL)
  NULLIFY(YBLD_DESC%XDESC_EMIS_ROOF)
  NULLIFY(YBLD_DESC%XDESC_EMIS_ROAD)
  NULLIFY(YBLD_DESC%XDESC_EMIS_WALL)
  NULLIFY(YBLD_DESC%XDESC_TCOOL_TARGET)
  NULLIFY(YBLD_DESC%XDESC_THEAT_TARGET)  
  NULLIFY(YBLD_DESC%XDESC_F_WASTE_CAN)  
  NULLIFY(YBLD_DESC%XDESC_EFF_HEAT)
  NULLIFY(YBLD_DESC%XDESC_HC_FLOOR)
  NULLIFY(YBLD_DESC%XDESC_TC_FLOOR)
  NULLIFY(YBLD_DESC%XDESC_D_FLOOR)
  NULLIFY(YBLD_DESC%XDESC_HC_ROOF)
  NULLIFY(YBLD_DESC%XDESC_TC_ROOF)
  NULLIFY(YBLD_DESC%XDESC_D_ROOF)
  NULLIFY(YBLD_DESC%XDESC_HC_ROAD)
  NULLIFY(YBLD_DESC%XDESC_TC_ROAD)
  NULLIFY(YBLD_DESC%XDESC_D_ROAD)
  NULLIFY(YBLD_DESC%XDESC_HC_WALL)
  NULLIFY(YBLD_DESC%XDESC_TC_WALL)
  NULLIFY(YBLD_DESC%XDESC_D_WALL)
  NULLIFY(YBLD_DESC%XDESC_QIN)
  NULLIFY(YBLD_DESC%XDESC_QIN_FRAD)
  NULLIFY(YBLD_DESC%XDESC_SHGC) 
  NULLIFY(YBLD_DESC%XDESC_U_WIN)
  NULLIFY(YBLD_DESC%XDESC_GR)
  NULLIFY(YBLD_DESC%XDESC_FLOOR_HEIGHT)
  NULLIFY(YBLD_DESC%XDESC_INF)
  NULLIFY(YBLD_DESC%XDESC_F_WATER_COND)
  NULLIFY(YBLD_DESC%XDESC_QIN_FLAT)
  NULLIFY(YBLD_DESC%XDESC_HR_TARGET)  
  NULLIFY(YBLD_DESC%XDESC_V_VENT)
  NULLIFY(YBLD_DESC%XDESC_COP_RAT)
  NULLIFY(YBLD_DESC%XDESC_GREENROOF)
  NULLIFY(YBLD_DESC%XDESC_SHADE)
  NULLIFY(YBLD_DESC%XDESC_NATVENT)
  NULLIFY(YBLD_DESC%XDESC_EMIS_PANEL)
  NULLIFY(YBLD_DESC%XDESC_ALB_PANEL)
  NULLIFY(YBLD_DESC%XDESC_EFF_PANEL)
  NULLIFY(YBLD_DESC%XDESC_FRAC_PANEL)
  NULLIFY(YBLD_DESC%XDESC_RESIDENTIAL)
YBLD_DESC%NDESC_BLD=0
YBLD_DESC%NDESC_AGE=0
YBLD_DESC%NDESC_CODE=0
YBLD_DESC%NDESC_USE=0
YBLD_DESC%NDESC_ROOF_LAYER=0
YBLD_DESC%NDESC_ROAD_LAYER=0
YBLD_DESC%NDESC_WALL_LAYER=0
YBLD_DESC%NDESC_FLOOR_LAYER=0
IF (LHOOK) CALL DR_HOOK("MODD_BLD_DESCRIPTION_n:BLD_DESC_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE BLD_DESC_INIT
!
!
END MODULE MODD_BLD_DESCRIPTION_n
