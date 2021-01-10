!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ################
      MODULE MODD_DATA_TEB_n
!     ################
!
!!****  *MODD_DATA_TEB_n - declaration of surface parameters for urban surface
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004
!!      Modified       08/2012 G. Pigeon ROUGH_WALL, ROUGH_ROOF
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE DATA_TEB_t
!
  LOGICAL :: LDATA_BLDTYPE 
  LOGICAL :: LDATA_BLD_AGE
  LOGICAL :: LDATA_USETYPE 
  LOGICAL :: LDATA_GARDEN
  LOGICAL :: LDATA_GREENROOF
  LOGICAL :: LDATA_ROAD_DIR
  LOGICAL :: LDATA_BLD
  LOGICAL :: LDATA_BLD_HEIGHT
  LOGICAL :: LDATA_WALL_O_HOR
  LOGICAL :: LDATA_Z0_TOWN
  LOGICAL :: LDATA_ALB_ROOF
  LOGICAL :: LDATA_EMIS_ROOF
  LOGICAL :: LDATA_HC_ROOF
  LOGICAL :: LDATA_TC_ROOF
  LOGICAL :: LDATA_D_ROOF
  LOGICAL :: LDATA_ALB_ROAD
  LOGICAL :: LDATA_EMIS_ROAD
  LOGICAL :: LDATA_HC_ROAD
  LOGICAL :: LDATA_TC_ROAD
  LOGICAL :: LDATA_D_ROAD
  LOGICAL :: LDATA_ALB_WALL
  LOGICAL :: LDATA_EMIS_WALL
  LOGICAL :: LDATA_HC_WALL
  LOGICAL :: LDATA_TC_WALL
  LOGICAL :: LDATA_D_WALL
  LOGICAL :: LDATA_H_TRAFFIC
  LOGICAL :: LDATA_LE_TRAFFIC
  LOGICAL :: LDATA_LE_INDUSTRY
  LOGICAL :: LDATA_H_INDUSTRY
  LOGICAL :: LDATA_ROUGH_ROOF
  LOGICAL :: LDATA_ROUGH_WALL
  LOGICAL :: LDATA_RESIDENTIAL
  LOGICAL :: LDATA_EMIS_PANEL
  LOGICAL :: LDATA_ALB_PANEL
  LOGICAL :: LDATA_EFF_PANEL
  LOGICAL :: LDATA_FRAC_PANEL
!
! Number of layers in the specification of thermal characteristics
!
  INTEGER                    :: NPAR_ROOF_LAYER   ! number of layers in roofs
  INTEGER                    :: NPAR_ROAD_LAYER   ! number of layers in roads
  INTEGER                    :: NPAR_WALL_LAYER   ! number of layers in walls
!
!
! Geometric Parameters:
!
  INTEGER, POINTER, DIMENSION(:):: NPAR_BLDTYPE      ! type of buidlings
  INTEGER, POINTER, DIMENSION(:):: NPAR_BLD_AGE      ! date of construction of buildings
  INTEGER, POINTER, DIMENSION(:):: NPAR_BLDCODE      ! code for buildings (type+age)
  INTEGER, POINTER, DIMENSION(:):: NPAR_USETYPE      ! type of use in the buildings
  REAL, POINTER, DIMENSION(:)   :: XPAR_GARDEN       ! fraction of veg in the streets   (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_GREENROOF    ! fraction of greenroofs on roofs  (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_ROAD_DIR     ! road direction (deg from North, clockwise)
  REAL, POINTER, DIMENSION(:)   :: XPAR_BLD          ! fraction of buildings            (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_BLD_HEIGHT   ! buildings height 'h'             (m)
  REAL, POINTER, DIMENSION(:)   :: XPAR_WALL_O_HOR   ! wall surf. / hor. surf.          (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_Z0_TOWN      ! roughness length for momentum    (m)
!
! Roof parameters
!
  REAL, POINTER, DIMENSION(:)   :: XPAR_ALB_ROOF     ! roof albedo                      (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_EMIS_ROOF    ! roof emissivity                  (-)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_HC_ROOF      ! roof layers heat capacity        (J/K/m3)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_TC_ROOF      ! roof layers thermal conductivity (W/K/m)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_D_ROOF       ! depth of roof layers             (m)
  REAL, POINTER, DIMENSION(:)   :: XPAR_ROUGH_ROOF   ! outside roof roughness coef 
!
!
! Road parameters
!
  REAL, POINTER, DIMENSION(:)   :: XPAR_ALB_ROAD     ! road albedo                      (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_EMIS_ROAD    ! road emissivity                  (-)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_HC_ROAD      ! road layers heat capacity        (J/K/m3)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_TC_ROAD      ! road layers thermal conductivity (W/K/m)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_D_ROAD       ! depth of road layers             (m)
!
! Wall parameters
!
  REAL, POINTER, DIMENSION(:)   :: XPAR_ALB_WALL     ! wall albedo                      (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_EMIS_WALL    ! wall emissivity                  (-)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_HC_WALL      ! wall layers heat capacity        (J/K/m3)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_TC_WALL      ! wall layers thermal conductivity (W/K/m)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_D_WALL       ! depth of wall layers             (m)
  REAL, POINTER, DIMENSION(:)   :: XPAR_ROUGH_WALL   ! outside wall roughness coef
!
! anthropogenic fluxes
!
  REAL, POINTER, DIMENSION(:)   :: XPAR_H_TRAFFIC    ! anthropogenic sensible
!                                                  ! heat fluxes due to traffic       (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XPAR_LE_TRAFFIC   ! anthropogenic latent
!                                                  ! heat fluxes due to traffic       (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XPAR_H_INDUSTRY   ! anthropogenic sensible                   
!                                                  ! heat fluxes due to factories     (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XPAR_LE_INDUSTRY  ! anthropogenic latent
!                                                  ! heat fluxes due to factories     (W/m2)
! solar panels
!
  REAL, POINTER, DIMENSION(:)   :: XPAR_EMIS_PANEL   ! emissivity of solar panel      (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_ALB_PANEL    ! albedo     of solar panel      (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_EFF_PANEL    ! efficiency of solar panel      (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_FRAC_PANEL   ! fraction   of solar panel      (-)
!
! fraction of residential use
!
  REAL, POINTER, DIMENSION(:)   :: XPAR_RESIDENTIAL  ! fraction of residential use    (-)
!
END TYPE DATA_TEB_t



CONTAINS

!




SUBROUTINE DATA_TEB_INIT(YDATA_TEB)
TYPE(DATA_TEB_t), INTENT(INOUT) :: YDATA_TEB
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DATA_TEB_N:DATA_TEB_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YDATA_TEB%NPAR_BLDTYPE)
  NULLIFY(YDATA_TEB%NPAR_BLD_AGE)
  NULLIFY(YDATA_TEB%NPAR_BLDCODE)
  NULLIFY(YDATA_TEB%NPAR_USETYPE)
  NULLIFY(YDATA_TEB%XPAR_GARDEN)
  NULLIFY(YDATA_TEB%XPAR_GREENROOF)
  NULLIFY(YDATA_TEB%XPAR_ROAD_DIR)
  NULLIFY(YDATA_TEB%XPAR_BLD)
  NULLIFY(YDATA_TEB%XPAR_BLD_HEIGHT)
  NULLIFY(YDATA_TEB%XPAR_WALL_O_HOR)
  NULLIFY(YDATA_TEB%XPAR_Z0_TOWN)
  NULLIFY(YDATA_TEB%XPAR_ALB_ROOF)
  NULLIFY(YDATA_TEB%XPAR_EMIS_ROOF)
  NULLIFY(YDATA_TEB%XPAR_HC_ROOF)
  NULLIFY(YDATA_TEB%XPAR_TC_ROOF)
  NULLIFY(YDATA_TEB%XPAR_D_ROOF)
  NULLIFY(YDATA_TEB%XPAR_ALB_ROAD)
  NULLIFY(YDATA_TEB%XPAR_EMIS_ROAD)
  NULLIFY(YDATA_TEB%XPAR_HC_ROAD)
  NULLIFY(YDATA_TEB%XPAR_TC_ROAD)
  NULLIFY(YDATA_TEB%XPAR_D_ROAD)
  NULLIFY(YDATA_TEB%XPAR_ALB_WALL)
  NULLIFY(YDATA_TEB%XPAR_EMIS_WALL)
  NULLIFY(YDATA_TEB%XPAR_HC_WALL)
  NULLIFY(YDATA_TEB%XPAR_TC_WALL)
  NULLIFY(YDATA_TEB%XPAR_D_WALL)
  NULLIFY(YDATA_TEB%XPAR_H_TRAFFIC)
  NULLIFY(YDATA_TEB%XPAR_LE_TRAFFIC)
  NULLIFY(YDATA_TEB%XPAR_H_INDUSTRY)
  NULLIFY(YDATA_TEB%XPAR_LE_INDUSTRY)
  NULLIFY(YDATA_TEB%XPAR_ROUGH_ROOF)
  NULLIFY(YDATA_TEB%XPAR_ROUGH_WALL)
  NULLIFY(YDATA_TEB%XPAR_RESIDENTIAL)
  NULLIFY(YDATA_TEB%XPAR_EMIS_PANEL)
  NULLIFY(YDATA_TEB%XPAR_ALB_PANEL)
  NULLIFY(YDATA_TEB%XPAR_EFF_PANEL)
  NULLIFY(YDATA_TEB%XPAR_EMIS_PANEL)
YDATA_TEB%LDATA_BLDTYPE=.FALSE.
YDATA_TEB%LDATA_BLD_AGE=.FALSE.
YDATA_TEB%LDATA_USETYPE=.FALSE.
YDATA_TEB%LDATA_GARDEN=.FALSE.
YDATA_TEB%LDATA_GREENROOF=.FALSE.
YDATA_TEB%LDATA_ROAD_DIR=.FALSE.
YDATA_TEB%LDATA_BLD=.FALSE.
YDATA_TEB%LDATA_BLD_HEIGHT=.FALSE.
YDATA_TEB%LDATA_WALL_O_HOR=.FALSE.
YDATA_TEB%LDATA_Z0_TOWN=.FALSE.
YDATA_TEB%LDATA_ALB_ROOF=.FALSE.
YDATA_TEB%LDATA_EMIS_ROOF=.FALSE.
YDATA_TEB%LDATA_HC_ROOF=.FALSE.
YDATA_TEB%LDATA_TC_ROOF=.FALSE.
YDATA_TEB%LDATA_D_ROOF=.FALSE.
YDATA_TEB%LDATA_ALB_ROAD=.FALSE.
YDATA_TEB%LDATA_EMIS_ROAD=.FALSE.
YDATA_TEB%LDATA_HC_ROAD=.FALSE.
YDATA_TEB%LDATA_TC_ROAD=.FALSE.
YDATA_TEB%LDATA_D_ROAD=.FALSE.
YDATA_TEB%LDATA_ALB_WALL=.FALSE.
YDATA_TEB%LDATA_EMIS_WALL=.FALSE.
YDATA_TEB%LDATA_HC_WALL=.FALSE.
YDATA_TEB%LDATA_TC_WALL=.FALSE.
YDATA_TEB%LDATA_D_WALL=.FALSE.
YDATA_TEB%LDATA_H_TRAFFIC=.FALSE.
YDATA_TEB%LDATA_LE_TRAFFIC=.FALSE.
YDATA_TEB%LDATA_H_INDUSTRY=.FALSE.
YDATA_TEB%LDATA_LE_INDUSTRY=.FALSE.
YDATA_TEB%LDATA_ROUGH_ROOF=.FALSE.
YDATA_TEB%LDATA_ROUGH_WALL=.FALSE.
YDATA_TEB%LDATA_RESIDENTIAL=.FALSE.
YDATA_TEB%LDATA_EMIS_PANEL=.FALSE.
YDATA_TEB%LDATA_ALB_PANEL=.FALSE.
YDATA_TEB%LDATA_EFF_PANEL=.FALSE.
YDATA_TEB%LDATA_FRAC_PANEL=.FALSE.
YDATA_TEB%NPAR_ROOF_LAYER=0
YDATA_TEB%NPAR_ROAD_LAYER=0
YDATA_TEB%NPAR_WALL_LAYER=0
IF (LHOOK) CALL DR_HOOK("MODD_DATA_TEB_N:DATA_TEB_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DATA_TEB_INIT


END MODULE MODD_DATA_TEB_n
