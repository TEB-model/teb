!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ##################
      MODULE MODD_DATA_ISBA_n
!     ##################
!
!!****  *MODD_DATA_ISBA - declaration of DATA surface parameters for ISBA scheme
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
!!      Original       05/2005
!!      P Samuelsson   02/2012  MEB
!!      B. Decharme    01/16 : Bug when vegetation veg, z0 and emis are imposed whith interactive vegetation
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE_ECOSG, NVEGTYPE_OLD
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE DATA_ISBA_t
!-------------------------------------------------------------------------------
!
! Mask and number of grid elements containing patches/tiles:
!
  REAL, POINTER, DIMENSION(:,:)  :: XPAR_VEGTYPE   ! fraction of each vegetation type for
!                                                  ! each grid mesh                          (-)
!
!-------------------------------------------------------------------------------
!
  INTEGER                       :: NVEGTYPE
!
  INTEGER                       :: NTIME               ! number of time data
!                                                      ! for VEG, LAI, EMIS, Z0
  LOGICAL                        :: LDATA_MIXPAR
!
  LOGICAL                        :: LDATA_VEGTYPE
!
  LOGICAL, POINTER, DIMENSION(:)   :: LDATA_LAI
  LOGICAL, POINTER, DIMENSION(:)   :: LDATA_VEG
  LOGICAL, POINTER, DIMENSION(:)   :: LDATA_Z0
  LOGICAL, POINTER, DIMENSION(:)   :: LDATA_EMIS
  LOGICAL, POINTER, DIMENSION(:)   :: LDATA_ALBNIR_VEG
  LOGICAL, POINTER, DIMENSION(:)   :: LDATA_ALBVIS_VEG
  LOGICAL, POINTER, DIMENSION(:)   :: LDATA_ALBUV_VEG 
  LOGICAL, POINTER, DIMENSION(:)   :: LDATA_IRRIG  
  LOGICAL, POINTER, DIMENSION(:)   :: LDATA_WATSUP  
  LOGICAL, POINTER, DIMENSION(:)   :: LDATA_ALBNIR_SOIL
  LOGICAL, POINTER, DIMENSION(:)   :: LDATA_ALBVIS_SOIL
  LOGICAL, POINTER, DIMENSION(:)   :: LDATA_ALBUV_SOIL
! - For multi-energy balance (MEB)
!
  LOGICAL, POINTER, DIMENSION(:)   :: LDATA_GNDLITTER
  LOGICAL, POINTER, DIMENSION(:)   :: LDATA_Z0LITTER
  LOGICAL, POINTER, DIMENSION(:)   :: LDATA_H_VEG  
  !
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_DG
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_ROOTFRAC
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_H_TREE
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_DICE  
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_GROUND_DEPTH
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_ROOT_DEPTH
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_ROOT_EXTINCTION
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_ROOT_LIN
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_RSMIN
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_GAMMA
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_WRMAX_CF
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_CV
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_Z0_O_Z0H
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_RGL
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_BSLAI
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_LAIMIN
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_SEFOLD
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_GMES
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_RE25
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_STRESS
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_F2I
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_GC
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_DMAX
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_CE_NITRO
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_CF_NITRO
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_CNA_NITRO
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_SEED_M
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_SEED_D
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_REAP_M
  LOGICAL, DIMENSION(MAX(NVEGTYPE_OLD,NVEGTYPE_ECOSG))   :: LDATA_REAP_D
!
  LOGICAL :: LDATA_CONDSAT
  LOGICAL :: LDATA_MPOTSAT
  LOGICAL :: LDATA_BCOEF
  LOGICAL :: LDATA_WWILT
  LOGICAL :: LDATA_WFC
  LOGICAL :: LDATA_WSAT
!
  LOGICAL   :: LIMP_VEG
  LOGICAL   :: LIMP_Z0
  LOGICAL   :: LIMP_EMIS
!  
! Input Parameters, per patch:
!
! - vegetation + bare soil:
!
  REAL, POINTER, DIMENSION(:,:) :: XPAR_Z0_O_Z0H         ! ratio of surface roughness lengths
!                                                      ! (momentum to heat)                      (-)
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_EMIS             ! surface emissivity                      (-)
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_Z0               ! surface roughness length                (m)
!
! - vegetation:
!
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_ALBNIR_VEG       ! vegetation near-infra-red albedo        (-)
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_ALBVIS_VEG       ! vegetation visible albedo               (-)
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_ALBUV_VEG        ! vegetation UV albedo                    (-)
!
! - vegetation: default option (Jarvis) and general parameters:
!
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_VEG            ! vegetation cover fraction               (-)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_WRMAX_CF         ! coefficient for maximum water 
!                                                      ! interception 
!                                                      ! storage capacity on the vegetation      (-)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_RSMIN            ! minimum stomatal resistance             (s/m)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_GAMMA            ! coefficient for the calculation
!                                                      ! of the surface stomatal
!                                                      ! resistance
  REAL, POINTER, DIMENSION(:,:) :: XPAR_CV               ! vegetation thermal inertia coefficient  (K m2/J)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_RGL              ! maximum solar radiation
!                                                      ! usable in photosynthesis                (W/m2)
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_ROOTFRAC       ! root fraction profile ('DIF' option)
!
  REAL, POINTER, DIMENSION(:,:) :: XPAR_ROOT_DEPTH       ! root depth ('DIF' option)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_ROOT_EXTINCTION  ! root extinction parameter ('DIF' option)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_ROOT_LIN         ! root linear parameter ('DIF' option)
!
! - For multi-energy balance (MEB)
!
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_GNDLITTER      ! ground litter fraction
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_Z0LITTER       ! ground litter roughness length
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_H_VEG          ! height of canopy vegetation
!  
!-------------------------------------------------------------------------------
!
! - vegetation: Ags parameters ('AGS', 'LAI', 'AST', 'LST', 'NIT', 'NCB' options)
!
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_BSLAI        ! ratio d(biomass)/d(lai)                 (kg/m2)
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_LAIMIN       ! minimum LAI (Leaf Area Index)           (m2/m2)
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_SEFOLD       ! e-folding time for senescence           (s)
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_H_TREE       ! height of trees                         (m)
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_GMES         ! mesophyll conductance                   (m s-1)
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_RE25         ! Ecosystem respiration parameter         (kg m2 s-1)
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags Stress parameters ('AST', 'LST', 'NIT', 'NCB' options)
!
  LOGICAL, POINTER, DIMENSION(:,:) :: LPAR_STRESS       ! vegetation response type to water
!                                                     ! stress (true:defensive false:offensive) (-)
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_F2I          ! critical normilized soil water 
!                                                     ! content for stress parameterisation
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_GC           ! cuticular conductance                   (m s-1)
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_DMAX         ! maximum air saturation deficit
!                                                     ! tolerate by vegetation                  (kg/kg)
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags Nitrogen-model parameters ('NIT', 'NCB' option)
!
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_CE_NITRO       ! leaf aera ratio sensitivity to 
!                                                       ! nitrogen concentration                (m2/kg)
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_CF_NITRO       ! lethal minimum value of leaf area
!                                                       ! ratio                                 (m2/kg)
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_CNA_NITRO      ! nitrogen concentration of active 
!                                                       ! biomass                               (kg/kg)
!
!-------------------------------------------------------------------------------
!
! - soil: primary parameters
!
  REAL, POINTER, DIMENSION(:,:,:)  :: XPAR_DG          ! soil layer depth                        (m)
!                                                      ! NOTE: in Force-Restore mode, the 
!                                                      ! uppermost layer thickness is superficial
!                                                      ! and is only explicitly used for soil 
!                                                      ! water phase changes                     (m)
!
  REAL, POINTER,DIMENSION(:,:)     :: XPAR_GROUND_DEPTH ! ground depth (DIF option)
!
  REAL, POINTER,DIMENSION(:,:)     :: XPAR_DICE        ! depth of the soil column for the calculation
!                                                        of the frozen soil fraction (m) (Force restore)
!
! - bare soil albedo
!
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_ALBNIR_SOIL      ! soil near-infra-red albedo        (-)
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_ALBVIS_SOIL      ! soil visible albedo               (-)
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_ALBUV_SOIL       ! soil UV albedo                    (-)
!
  REAL, POINTER, DIMENSION(:,:) :: XPAR_CONDSAT
  REAL, POINTER, DIMENSION(:,:) :: XPAR_MPOTSAT
  REAL, POINTER, DIMENSION(:,:) :: XPAR_BCOEF
  REAL, POINTER, DIMENSION(:,:) :: XPAR_WWILT
  REAL, POINTER, DIMENSION(:,:) :: XPAR_WFC
  REAL, POINTER, DIMENSION(:,:) :: XPAR_WSAT
!
!-------------------------------------------------------------------------------
!
! - Vegetation: Ags Prognostic (YPHOTO = ('NIT' or 'NCB') or prescribed (YPHOTO='NON' or 'AST')
!
  REAL, POINTER, DIMENSION(:,:,:)     :: XPAR_LAI          ! Leaf Area Index                         (m2/m2)
!
  REAL, POINTER, DIMENSION(:,:,:)     :: XPAR_IRRIG
  REAL, POINTER, DIMENSION(:,:,:)     :: XPAR_WATSUP
  REAL, POINTER, DIMENSION(:,:)     :: XPAR_SEED_M
  REAL, POINTER, DIMENSION(:,:)     :: XPAR_SEED_D
  REAL, POINTER, DIMENSION(:,:)     :: XPAR_REAP_M
  REAL, POINTER, DIMENSION(:,:)     :: XPAR_REAP_D  
!
!-------------------------------------------------------------------------------
!
  REAL, POINTER, DIMENSION(:)   :: XPAR_FRAC_HVEG        ! fraction of high vegetation
  REAL, POINTER, DIMENSION(:)   :: XPAR_FRAC_LVEG        ! fraction of low  vegetation
  REAL, POINTER, DIMENSION(:)   :: XPAR_FRAC_NVEG        ! fraction of bare soil
  REAL, POINTER, DIMENSION(:,:) :: XPAR_LAI_HVEG         ! LAI      of high vegetation
  REAL, POINTER, DIMENSION(:,:) :: XPAR_LAI_LVEG         ! LAI      of low  vegetation
  REAL, POINTER, DIMENSION(:)   :: XPAR_H_HVEG           ! height of trees
!
END TYPE DATA_ISBA_t



CONTAINS

!




SUBROUTINE DATA_ISBA_INIT(YDATA_ISBA)
TYPE(DATA_ISBA_t), INTENT(INOUT) :: YDATA_ISBA
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DATA_ISBA_N:DATA_ISBA_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YDATA_ISBA%XPAR_VEGTYPE)
  NULLIFY(YDATA_ISBA%XPAR_Z0_O_Z0H)
  NULLIFY(YDATA_ISBA%XPAR_EMIS)
  NULLIFY(YDATA_ISBA%XPAR_Z0)
  NULLIFY(YDATA_ISBA%XPAR_ALBNIR_VEG)
  NULLIFY(YDATA_ISBA%XPAR_ALBVIS_VEG)
  NULLIFY(YDATA_ISBA%XPAR_ALBUV_VEG)
  NULLIFY(YDATA_ISBA%XPAR_VEG)
  NULLIFY(YDATA_ISBA%XPAR_WRMAX_CF)
  NULLIFY(YDATA_ISBA%XPAR_RSMIN)
  NULLIFY(YDATA_ISBA%XPAR_GAMMA)
  NULLIFY(YDATA_ISBA%XPAR_CV)
  NULLIFY(YDATA_ISBA%XPAR_RGL)
  NULLIFY(YDATA_ISBA%XPAR_ROOTFRAC)
  NULLIFY(YDATA_ISBA%XPAR_BSLAI)
  NULLIFY(YDATA_ISBA%XPAR_LAIMIN)
  NULLIFY(YDATA_ISBA%XPAR_SEFOLD)
  NULLIFY(YDATA_ISBA%XPAR_H_TREE)
  NULLIFY(YDATA_ISBA%XPAR_GMES)
  NULLIFY(YDATA_ISBA%XPAR_RE25)
  NULLIFY(YDATA_ISBA%LPAR_STRESS)
  NULLIFY(YDATA_ISBA%XPAR_F2I)
  NULLIFY(YDATA_ISBA%XPAR_GC)
  NULLIFY(YDATA_ISBA%XPAR_DMAX)
  NULLIFY(YDATA_ISBA%XPAR_CE_NITRO)
  NULLIFY(YDATA_ISBA%XPAR_CF_NITRO)
  NULLIFY(YDATA_ISBA%XPAR_CNA_NITRO)
  NULLIFY(YDATA_ISBA%XPAR_DG)
  NULLIFY(YDATA_ISBA%XPAR_DICE)
  NULLIFY(YDATA_ISBA%XPAR_GROUND_DEPTH)  
  NULLIFY(YDATA_ISBA%XPAR_ROOT_DEPTH)  
  NULLIFY(YDATA_ISBA%XPAR_ROOT_EXTINCTION)  
  NULLIFY(YDATA_ISBA%XPAR_ROOT_LIN)
  NULLIFY(YDATA_ISBA%XPAR_ALBNIR_SOIL)
  NULLIFY(YDATA_ISBA%XPAR_ALBVIS_SOIL)
  NULLIFY(YDATA_ISBA%XPAR_ALBUV_SOIL)
  NULLIFY(YDATA_ISBA%XPAR_LAI)
  NULLIFY(YDATA_ISBA%XPAR_IRRIG)
  NULLIFY(YDATA_ISBA%XPAR_WATSUP)
  NULLIFY(YDATA_ISBA%XPAR_SEED_M)
  NULLIFY(YDATA_ISBA%XPAR_SEED_D)
  NULLIFY(YDATA_ISBA%XPAR_REAP_M)
  NULLIFY(YDATA_ISBA%XPAR_REAP_D)  
  NULLIFY(YDATA_ISBA%XPAR_GNDLITTER)
  NULLIFY(YDATA_ISBA%XPAR_Z0LITTER)
  NULLIFY(YDATA_ISBA%XPAR_H_VEG)
  NULLIFY(YDATA_ISBA%XPAR_FRAC_HVEG)
  NULLIFY(YDATA_ISBA%XPAR_FRAC_LVEG)
  NULLIFY(YDATA_ISBA%XPAR_FRAC_NVEG)
  NULLIFY(YDATA_ISBA%XPAR_LAI_HVEG)
  NULLIFY(YDATA_ISBA%XPAR_LAI_LVEG)
  NULLIFY(YDATA_ISBA%XPAR_H_HVEG)
  NULLIFY(YDATA_ISBA%XPAR_CONDSAT)
  NULLIFY(YDATA_ISBA%XPAR_MPOTSAT)
  NULLIFY(YDATA_ISBA%XPAR_BCOEF)
  NULLIFY(YDATA_ISBA%XPAR_WWILT)
  NULLIFY(YDATA_ISBA%XPAR_WFC)
  NULLIFY(YDATA_ISBA%XPAR_WSAT)
!
YDATA_ISBA%NVEGTYPE=0
YDATA_ISBA%NTIME=0
YDATA_ISBA%LDATA_MIXPAR=.FALSE.
YDATA_ISBA%LDATA_VEGTYPE=.FALSE. 
!
NULLIFY(YDATA_ISBA%LDATA_LAI)
NULLIFY(YDATA_ISBA%LDATA_VEG)
NULLIFY(YDATA_ISBA%LDATA_Z0)
NULLIFY(YDATA_ISBA%LDATA_EMIS)
NULLIFY(YDATA_ISBA%LDATA_ALBNIR_VEG)
NULLIFY(YDATA_ISBA%LDATA_ALBVIS_VEG)
NULLIFY(YDATA_ISBA%LDATA_ALBUV_VEG)
NULLIFY(YDATA_ISBA%LDATA_ALBNIR_SOIL)
NULLIFY(YDATA_ISBA%LDATA_ALBVIS_SOIL)
NULLIFY(YDATA_ISBA%LDATA_ALBUV_SOIL)
NULLIFY(YDATA_ISBA%LDATA_GNDLITTER)
NULLIFY(YDATA_ISBA%LDATA_Z0LITTER)
NULLIFY(YDATA_ISBA%LDATA_H_VEG)
NULLIFY(YDATA_ISBA%LDATA_IRRIG)
NULLIFY(YDATA_ISBA%LDATA_WATSUP)
!
YDATA_ISBA%LDATA_DG(:) = .FALSE.
YDATA_ISBA%LDATA_ROOTFRAC(:) = .FALSE.
YDATA_ISBA%LDATA_H_TREE(:) = .FALSE.
YDATA_ISBA%LDATA_DICE(:) = .FALSE.
YDATA_ISBA%LDATA_GROUND_DEPTH(:) = .FALSE.
YDATA_ISBA%LDATA_ROOT_DEPTH(:) = .FALSE.
YDATA_ISBA%LDATA_ROOT_EXTINCTION(:) = .FALSE.
YDATA_ISBA%LDATA_ROOT_LIN(:) = .FALSE.
YDATA_ISBA%LDATA_RSMIN(:) = .FALSE.
YDATA_ISBA%LDATA_GAMMA(:) = .FALSE.
YDATA_ISBA%LDATA_WRMAX_CF(:) = .FALSE.
YDATA_ISBA%LDATA_CV(:) = .FALSE.
YDATA_ISBA%LDATA_RGL(:) = .FALSE.
YDATA_ISBA%LDATA_Z0_O_Z0H(:) = .FALSE.
YDATA_ISBA%LDATA_BSLAI(:) = .FALSE.
YDATA_ISBA%LDATA_LAIMIN(:) = .FALSE.
YDATA_ISBA%LDATA_SEFOLD(:) = .FALSE.
YDATA_ISBA%LDATA_GMES(:) = .FALSE.
YDATA_ISBA%LDATA_RE25(:) = .FALSE.
YDATA_ISBA%LDATA_STRESS(:) = .FALSE.
YDATA_ISBA%LDATA_F2I(:) = .FALSE.
YDATA_ISBA%LDATA_GC(:) = .FALSE.
YDATA_ISBA%LDATA_DMAX(:) = .FALSE.
YDATA_ISBA%LDATA_CE_NITRO(:) = .FALSE.
YDATA_ISBA%LDATA_CF_NITRO(:) = .FALSE.
YDATA_ISBA%LDATA_CNA_NITRO(:) = .FALSE.
YDATA_ISBA%LDATA_SEED_M(:) = .FALSE.
YDATA_ISBA%LDATA_SEED_D(:) = .FALSE.
YDATA_ISBA%LDATA_REAP_M(:) = .FALSE.
YDATA_ISBA%LDATA_REAP_D(:) = .FALSE.
!
YDATA_ISBA%LDATA_CONDSAT = .FALSE.
YDATA_ISBA%LDATA_MPOTSAT = .FALSE.
YDATA_ISBA%LDATA_BCOEF = .FALSE.
YDATA_ISBA%LDATA_WWILT = .FALSE.
YDATA_ISBA%LDATA_WFC = .FALSE.
YDATA_ISBA%LDATA_WSAT = .FALSE.
!
YDATA_ISBA%LIMP_VEG=.FALSE.
YDATA_ISBA%LIMP_Z0=.FALSE.
YDATA_ISBA%LIMP_EMIS=.FALSE.
!
IF (LHOOK) CALL DR_HOOK("MODD_DATA_ISBA_N:DATA_ISBA_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DATA_ISBA_INIT


END MODULE MODD_DATA_ISBA_n
