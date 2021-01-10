!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!##################
MODULE MODD_ISBA_n
!##################
!
!!****  *MODD_ISBA - declaration of packed surface parameters for ISBA scheme
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
!!      A. Boone   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       20/09/02
!!      A.L. Gibelin    04/2009 : BIOMASS and RESP_BIOMASS arrays 
!!      A.L. Gibelin    04/2009 : TAU_WOOD for NCB option 
!!      A.L. Gibelin    05/2009 : Add carbon spinup
!!      A.L. Gibelin    06/2009 : Soil carbon variables for CNT option
!!      A.L. Gibelin    07/2009 : Suppress RDK and transform GPP as a diagnostic
!!      A.L. Gibelin    07/2009 : Suppress PPST and PPSTF as outputs
!!      P. Samuelsson   02/2012 : MEB
!!      B. Decharme    10/2016  bug surface/groundwater coupling 
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_SNOW
USE MODD_TYPE_DATE_SURF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE ISBA_S_t
!
! General surface parameters:
!
REAL, POINTER, DIMENSION(:)   :: XZS               ! relief                                  (m)
REAL, POINTER, DIMENSION(:,:) :: XCOVER            ! fraction of each ecosystem              (-)
LOGICAL, POINTER, DIMENSION(:):: LCOVER            ! GCOVER(i)=T --> ith cover field is not 0.
!
! Topmodel statistics
!
REAL, POINTER, DIMENSION(:)      :: XTI_MIN,XTI_MAX,XTI_MEAN,XTI_STD,XTI_SKEW
!
REAL, POINTER, DIMENSION(:,:)    :: XSOC           ! soil organic carbon content             (kg/m2)
REAL, POINTER, DIMENSION(:)      :: XPH            ! soil pH
REAL, POINTER, DIMENSION(:)      :: XFERT          ! soil fertilisation rate (kgN/ha/h)
!
!
REAL, POINTER, DIMENSION(:)      :: XABC           ! abscissa needed for integration
!                                                  ! of net assimilation and stomatal
!                                                  ! conductance over canopy depth           (-)
REAL, POINTER, DIMENSION(:)      :: XPOI           ! Gaussian weights for integration
!                                                  ! of net assimilation and stomatal
!                                                  ! conductance over canopy depth           (-)
!
TYPE (DATE_TIME) :: TTIME
!
REAL, POINTER, DIMENSION(:,:) :: XTAB_FSAT !Satured fraction array
REAL, POINTER, DIMENSION(:,:) :: XTAB_WTOP !Active TOPMODEL-layer array
REAL, POINTER, DIMENSION(:,:) :: XTAB_QTOP !Subsurface flow TOPMODEL array
!
REAL, POINTER, DIMENSION(:) :: XF_PARAM
REAL, POINTER, DIMENSION(:) :: XC_DEPTH_RATIO
!
! - Coupling with river routing model
!  
REAL, POINTER, DIMENSION(:)  :: XCPL_DRAIN   ! Surface runoff
REAL, POINTER, DIMENSION(:)  :: XCPL_RUNOFF  ! Deep drainage or gourdwater recharge
REAL, POINTER, DIMENSION(:)  :: XCPL_ICEFLUX ! Calving flux
REAL, POINTER, DIMENSION(:)  :: XCPL_EFLOOD  ! floodplains evaporation
REAL, POINTER, DIMENSION(:)  :: XCPL_PFLOOD  ! floodplains precipitation interception
REAL, POINTER, DIMENSION(:)  :: XCPL_IFLOOD  ! floodplains infiltration
!
!  - Random perturbations
!
REAL, POINTER, DIMENSION(:)     :: XPERTVEG
REAL, POINTER, DIMENSION(:)     :: XPERTLAI
REAL, POINTER, DIMENSION(:)     :: XPERTCV
REAL, POINTER, DIMENSION(:)     :: XPERTALB
REAL, POINTER, DIMENSION(:)     :: XPERTZ0
!
REAL, POINTER, DIMENSION(:)   :: XTSRAD_NAT        ! patch averaged radiative temperature    (K)
!
REAL, POINTER, DIMENSION(:)   :: XEMIS_NAT         ! patch averaged emissivity               (-)
!
!  - Assimilation: ENKF
!
REAL, POINTER, DIMENSION(:,:) :: XFRACSOC ! Fraction of organic carbon in each soil layer
!
REAL, POINTER, DIMENSION(:,:) :: XVEGTYPE
!
REAL, POINTER, DIMENSION(:,:)    :: XPATCH         ! fraction of each tile/patch   (-)
!
! Mask and number of grid elements containing patches/tiles:
!
REAL, POINTER, DIMENSION(:,:,:)  :: XVEGTYPE_PATCH ! fraction of each vegetation type for
!
REAL, POINTER, DIMENSION(:,:)     :: XINNOV
REAL, POINTER, DIMENSION(:,:)     :: XRESID
!
REAL, POINTER, DIMENSION(:,:) :: XWORK_WR
!
REAL, POINTER, DIMENSION(:,:,:) :: XWSN_WR
REAL, POINTER, DIMENSION(:,:,:) :: XRHO_WR
REAL, POINTER, DIMENSION(:,:,:) :: XHEA_WR
REAL, POINTER, DIMENSION(:,:,:) :: XAGE_WR
REAL, POINTER, DIMENSION(:,:,:) :: XSG1_WR
REAL, POINTER, DIMENSION(:,:,:) :: XSG2_WR
REAL, POINTER, DIMENSION(:,:,:) :: XHIS_WR
REAL, POINTER, DIMENSION(:,:,:) :: XT_WR
REAL, POINTER, DIMENSION(:,:) :: XALB_WR
!
TYPE(DATE_TIME), POINTER, DIMENSION(:,:) :: TDATE_WR
!
END TYPE ISBA_S_t
!
!
TYPE ISBA_K_t
!
REAL, POINTER, DIMENSION(:,:)    :: XSAND          ! sand fraction                           (-)
REAL, POINTER, DIMENSION(:,:)    :: XCLAY          ! clay fraction                           (-)
!
REAL, POINTER, DIMENSION(:)      :: XPERM          ! permafrost distribution                 (-)
!
REAL, POINTER, DIMENSION(:)      :: XRUNOFFB       ! sub-grid dt92 surface runoff slope parameter (-)  
REAL, POINTER, DIMENSION(:)      :: XWDRAIN        ! continuous drainage parameter           (-)
!
!
REAL, POINTER, DIMENSION(:)      :: XTDEEP         ! prescribed deep soil temperature 
!                                                  ! (optional)
REAL, POINTER, DIMENSION(:)      :: XGAMMAT        ! 'Force-Restore' timescale when using a
!                                                  ! prescribed lower boundary temperature   (1/days)
!
REAL, POINTER, DIMENSION(:,:)    :: XMPOTSAT       ! matric potential at saturation          (m)
REAL, POINTER, DIMENSION(:,:)    :: XBCOEF         ! soil water CH78 b-parameter             (-)
REAL, POINTER, DIMENSION(:,:)    :: XWWILT         ! wilting point volumetric water content 
!                                                  ! profile                                 (m3/m3)
REAL, POINTER, DIMENSION(:,:)    :: XWFC           ! field capacity volumetric water content
!                                                  ! profile                                 (m3/m3)
REAL, POINTER, DIMENSION(:,:)    :: XWSAT          ! porosity profile                        (m3/m3) 
!
REAL, POINTER, DIMENSION(:)      :: XCGSAT         ! soil thermal inertia coefficient at 
!                                                  ! saturation
REAL, POINTER, DIMENSION(:)      :: XC4B           ! 'Force-Restore' sub-surface vertical 
!                                                  ! diffusion coefficient (slope parameter) (-)
REAL, POINTER, DIMENSION(:)      :: XACOEF         ! 'Force-Restore' surface vertical 
!                                                  ! diffusion coefficient                   (-)
REAL, POINTER, DIMENSION(:)      :: XPCOEF         ! 'Force-Restore' surface vertical 
!                                                  ! diffusion coefficient                   (-)
REAL, POINTER, DIMENSION(:,:)    :: XHCAPSOIL      ! soil heat capacity                      (J/K/m3)
REAL, POINTER, DIMENSION(:,:)    :: XCONDDRY       ! soil dry thermal conductivity           (W/m/K)
REAL, POINTER, DIMENSION(:,:)    :: XCONDSLD       ! soil solids thermal conductivity        (W/m/K)
!
! - Water table depth coupling
!  
REAL, POINTER, DIMENSION(:)  :: XFWTD         ! grid-cell fraction of water table rise
REAL, POINTER, DIMENSION(:)  :: XWTD          ! water table depth (negative below soil surface) (m)
!
! * Physiographic radiative fields
!
REAL, POINTER, DIMENSION(:)   :: XALBNIR_DRY       ! dry soil near-infra-red albedo          (-)
REAL, POINTER, DIMENSION(:)   :: XALBVIS_DRY       ! dry soil visible albedo                 (-)
REAL, POINTER, DIMENSION(:)   :: XALBUV_DRY        ! dry soil UV albedo                      (-)
REAL, POINTER, DIMENSION(:)   :: XALBNIR_WET       ! wet soil near-infra-red albedo          (-)
REAL, POINTER, DIMENSION(:)   :: XALBVIS_WET       ! wet soil visible albedo                 (-)
REAL, POINTER, DIMENSION(:)   :: XALBUV_WET        ! wet soil UV albedo                      (-)
!
! * SGH initializations
!
REAL, POINTER, DIMENSION(:,:)  :: XWD0     ! water content equivalent to TOPMODEL maximum deficit
REAL, POINTER, DIMENSION(:,:)  :: XKANISO  ! Anisotropy coeficient for hydraulic conductivity
!
REAL, POINTER, DIMENSION(:)    :: XMUF     ! fraction of the grid cell reached by the rainfall
REAL, POINTER, DIMENSION(:)    :: XFSAT    ! Topmodel or dt92 saturated fracti
!
REAL, POINTER, DIMENSION(:)  :: XFFLOOD      ! Grid-cell flood fraction
REAL, POINTER, DIMENSION(:)  :: XPIFLOOD     ! flood potential infiltration (kg/m2/s)
!
! - Flood scheme
!
REAL, POINTER, DIMENSION(:) :: XFF         ! Total Flood fraction  
REAL, POINTER, DIMENSION(:) :: XFFG        ! Flood fraction over ground
REAL, POINTER, DIMENSION(:) :: XFFV        ! Flood fraction over vegetation
REAL, POINTER, DIMENSION(:) :: XFFROZEN    ! Fraction of frozen floodplains
REAL, POINTER, DIMENSION(:) :: XALBF       ! Flood albedo
REAL, POINTER, DIMENSION(:) :: XEMISF      ! Flood emissivity
!
! - Snow and flood fractions and total albedo at time t:                             (-)
!
REAL, POINTER, DIMENSION(:,:) :: XDIR_ALB_WITH_SNOW ! total direct albedo by bands
REAL, POINTER, DIMENSION(:,:) :: XSCA_ALB_WITH_SNOW ! total diffuse albedo by bands
!
REAL, POINTER, DIMENSION(:,:) :: XVEGTYPE
!
END TYPE ISBA_K_t
!
!
TYPE ISBA_P_t
!
INTEGER  :: NSIZE_P ! number of sub-patchs/tiles              (-)
!
REAL, POINTER, DIMENSION(:)    :: XPATCH         ! fraction of each tile/patch   (-)
!
! Mask and number of grid elements containing patches/tiles:
!
REAL, POINTER, DIMENSION(:,:)  :: XVEGTYPE_PATCH ! fraction of each vegetation type for
!
INTEGER, POINTER, DIMENSION(:) :: NR_P    ! patch/tile mask  
!
REAL, POINTER, DIMENSION(:)    :: XPATCH_OLD     ! fraction of each tile/patchfor land use (-)
!
!
REAL, POINTER, DIMENSION(:)    :: XANMAX         ! maximum photosynthesis rate             (
REAL, POINTER, DIMENSION(:)    :: XFZERO         ! ideal value of F, no photo- 
!                                                ! respiration or saturation deficit  
REAL, POINTER, DIMENSION(:)    :: XEPSO          ! maximum initial quantum use             
!                                                ! efficiency                              (mg J-1 PAR)
REAL, POINTER, DIMENSION(:)    :: XGAMM          ! CO2 conpensation concentration          (ppm)
REAL, POINTER, DIMENSION(:)    :: XQDGAMM        ! Log of Q10 function for CO2 conpensation 
!                                               ! concentration                           (-)
REAL, POINTER, DIMENSION(:)    :: XQDGMES        ! Log of Q10 function for mesophyll conductance  (-)
REAL, POINTER, DIMENSION(:)    :: XT1GMES        ! reference temperature for computing 
!                                                ! compensation concentration function for 
!                                                ! mesophyll conductance: minimum
!                                                ! temperature                             (K)
REAL, POINTER, DIMENSION(:)    :: XT2GMES        ! reference temperature for computing 
!                                                ! compensation concentration function for 
!                                                ! mesophyll conductance: maximum
!                                                ! temperature                             (K)
REAL, POINTER, DIMENSION(:)    :: XAMAX          ! leaf photosynthetic capacity            (mg m-2 s-1)
REAL, POINTER, DIMENSION(:)    :: XQDAMAX        ! Log of Q10 function for leaf photosynthetic 
!                                                ! capacity                                (-)
REAL, POINTER, DIMENSION(:)    :: XT1AMAX        ! reference temperature for computing 
!                                                ! compensation concentration function for 
!                                                ! leaf photosynthetic capacity: minimum
!                                                ! temperature                             (K)
REAL, POINTER, DIMENSION(:)    :: XT2AMAX        ! reference temperature for computing 
!                                                ! compensation concentration function for 
!                                                ! leaf photosynthetic capacity: maximum
!                                                ! temperature                             (K)
REAL, POINTER, DIMENSION(:)    :: XAH            ! coefficients for herbaceous water stress 
!                                                ! response (offensive or defensive)       (log(mm/s))
REAL, POINTER, DIMENSION(:)    :: XBH            ! coefficients for herbaceous water stress 
!                                                ! response (offensive or defensive)
REAL, POINTER, DIMENSION(:)    :: XTAU_WOOD      ! residence time in woody biomass         (s)
REAL, POINTER, DIMENSION(:,:)   :: XINCREASE     ! biomass increase                     (kg/m2/day)
REAL, POINTER, DIMENSION(:,:)   :: XTURNOVER     ! turnover rates from biomass to litter (gC/m2/s)
!
! *Soil hydraulic characteristics
!
REAL, POINTER, DIMENSION(:,:)  :: XCONDSAT       ! hydraulic conductivity at saturation    (m/s)
!
REAL, POINTER, DIMENSION(:)      :: XTAUICE        ! soil freezing characteristic timescale  (s)
!
REAL, POINTER, DIMENSION(:)    :: XC1SAT         ! 'Force-Restore' C1 coefficient at 
!                                                ! saturation                              (-)
REAL, POINTER, DIMENSION(:)    :: XC2REF         ! 'Force-Restore' reference value of C2   (-)
REAL, POINTER, DIMENSION(:,:)  :: XC3            ! 'Force-Restore' C3 drainage coefficient (m)
REAL, POINTER, DIMENSION(:)    :: XC4REF         ! 'Force-Restore' sub-surface vertical 
!                                                ! for lateral drainage ('DIF' option)
!
REAL, POINTER, DIMENSION(:)    :: XBSLAI_NITRO   ! biomass/LAI ratio from nitrogen 
!                                                  ! decline theory                        (kg/m2)
! * Soil thermal characteristics
!
REAL, POINTER, DIMENSION(:)    :: XCPS
REAL, POINTER, DIMENSION(:)    :: XLVTT
REAL, POINTER, DIMENSION(:)    :: XLSTT
!
! * Initialize hydrology
!
  REAL, POINTER, DIMENSION(:)    :: XRUNOFFD     ! depth over which sub-grid runoff is
!                                                ! computed: in Force-Restore this is the
!                                                ! total soil column ('2-L'), or root zone
!                                                ! ('3-L'). For the 'DIF' option, it can
!                                                ! be any depth within soil column         (m)
!
REAL, POINTER, DIMENSION(:,:)  :: XDZG           ! soil layers thicknesses (DIF option)
REAL, POINTER, DIMENSION(:,:)  :: XDZDIF         ! distance between consecuative layer mid-points (DIF option)
REAL, POINTER, DIMENSION(:,:)  :: XSOILWGHT      ! VEG-DIF: weights for vertical
!                                                ! integration of soil water and properties
!
REAL, POINTER, DIMENSION(:)  :: XKSAT_ICE        ! hydraulic conductivity at saturation
!                                                    over frozen area (m s-1)
!
REAL, POINTER, DIMENSION(:,:) :: XTOPQS  ! Topmodel subsurface flow by layer (m/s)
!
REAL, POINTER, DIMENSION(:,:) :: XDG           ! soil layer depth                  (m)
!                                                ! NOTE: in Force-Restore mode, the 
!                                                ! uppermost layer depth is superficial
!                                                ! and is only explicitly used for soil 
!                                                ! water phase changes                     (m)
!
REAL, POINTER, DIMENSION(:,:)  :: XDG_OLD      ! For land use
REAL, POINTER, DIMENSION(:)    :: XDG2
INTEGER, POINTER, DIMENSION(:) :: NWG_LAYER    ! Number of soil moisture layers for DIF
REAL, POINTER, DIMENSION(:)    :: XDROOT       ! effective root depth for DIF (m)
REAL, POINTER, DIMENSION(:,:)  :: XROOTFRAC    ! root fraction profile ('DIF' option)
!
REAL, POINTER, DIMENSION(:) :: XD_ICE          ! depth of the soil column for the calculation
!                                                 of the frozen soil fraction (m)
!
REAL, POINTER, DIMENSION(:) :: XH_TREE         ! height of trees                         (m)
!
REAL, POINTER, DIMENSION(:) :: XZ0_O_Z0H       ! ratio of surface roughness lengths
!                                                ! (momentum to heat)                      (-)
!
REAL, POINTER, DIMENSION(:) :: XRE25           ! Ecosystem respiration parameter         (kg/kg.m.s-1)
!
REAL, POINTER, DIMENSION(:) :: XDMAX           ! maximum air saturation deficit
!                                                ! tolerate by vegetation
!                                                (kg/kg)
!
REAL, POINTER, DIMENSION(:,:)     :: XRED_NOISE
REAL, POINTER, DIMENSION(:,:)     :: XINCR
REAL, POINTER, DIMENSION(:,:,:)   :: XHO
!
END TYPE ISBA_P_t
!
TYPE ISBA_PE_t
!
! Prognostic variables:
!
! - Soil and vegetation heat and water:
!
REAL, POINTER, DIMENSION(:,:) :: XWG           ! soil volumetric water content profile   (m3/m3)
REAL, POINTER, DIMENSION(:,:) :: XWGI          ! soil liquid water equivalent volumetric 
!                                                ! ice content profile                     (m3/m3)
REAL, POINTER, DIMENSION(:)   :: XWR           ! liquid water retained on the
!                                                ! foliage of the vegetation
!                                                ! canopy                                  (kg/m2)
REAL, POINTER, DIMENSION(:,:) :: XTG           ! surface and sub-surface soil 
!                                                ! temperature profile                     (K)
!
! - Snow Cover:
!
TYPE(SURF_SNOW) :: TSNOW                         ! snow state: 
!                                                ! scheme type/option                      (-)
!                                                ! number of layers                        (-)
!                                                ! snow (& liq. water) content             (kg/m2)
!                                                ! heat content                            (J/m2)
!                                                ! temperature                             (K)
!                                                ! density                                 (kg m-3)
!
REAL, POINTER, DIMENSION(:) :: XICE_STO        ! Glacier ice storage reservoir
!
! - For multi-energy balance:
!
REAL, POINTER, DIMENSION(:) :: XWRL            ! liquid water retained on litter          (kg/m2)
REAL, POINTER, DIMENSION(:) :: XWRLI           ! ice retained on litter          (kg/m2)
REAL, POINTER, DIMENSION(:) :: XWRVN           ! snow retained on the foliage
!                                                ! of the canopy vegetation                  (kg/m2)
REAL, POINTER, DIMENSION(:) :: XTV             ! canopy vegetation temperature             (K)
REAL, POINTER, DIMENSION(:) :: XTL             ! litter temperature             (K)
REAL, POINTER, DIMENSION(:) :: XTC             ! canopy air temperature                    (K)
REAL, POINTER, DIMENSION(:) :: XQC             ! canopy air specific humidity              (kg/kg)
!
! * Half prognostic fields
!
REAL, POINTER, DIMENSION(:)     :: XRESA         ! aerodynamic resistance                  (s/m)
!
! - Vegetation: Ags Prognostic (YPHOTO = 'AST', 'NIT', 'NCB')
!
REAL, POINTER, DIMENSION(:) :: XAN           ! net CO2 assimilation                    (mg/m2/s)
REAL, POINTER, DIMENSION(:) :: XANDAY        ! daily net CO2 assimilation              (mg/m2)
REAL, POINTER, DIMENSION(:) :: XANFM         ! maximum leaf assimilation               (mg/m2/s)
REAL, POINTER, DIMENSION(:) :: XLE           ! evapotranspiration                      (W/m2)
!
REAL, POINTER, DIMENSION(:) :: XFAPARC       ! Fapar of vegetation (cumul)
REAL, POINTER, DIMENSION(:) :: XFAPIRC       ! Fapir of vegetation (cumul)
REAL, POINTER, DIMENSION(:) :: XLAI_EFFC     ! Effective LAI (cumul)
REAL, POINTER, DIMENSION(:) :: XMUS          ! cos zenithal angle (cumul)
!
REAL, POINTER, DIMENSION(:,:) :: XRESP_BIOMASS    ! daily cumulated respiration of 
!                                                   ! biomass                              (kg/m2/s)
REAL, POINTER, DIMENSION(:,:) :: XBIOMASS         ! biomass of previous day              (kg/m2) 
!
! - Soil carbon (ISBA-CC, YRESPSL = 'CNT')
!
REAL, POINTER, DIMENSION(:,:,:) :: XLITTER          ! litter pools                         (gC/m2)
REAL, POINTER, DIMENSION(:,:)   :: XSOILCARB        ! soil carbon pools                    (gC/m2) 
REAL, POINTER, DIMENSION(:,:)   :: XLIGNIN_STRUC    ! ratio Lignin/Carbon in structural
!                                                       litter                               (gC/m2)
!
REAL, POINTER, DIMENSION(:) :: XPSNG         ! Snow fraction over ground
REAL, POINTER, DIMENSION(:) :: XPSNV         ! Snow fraction over vegetation
REAL, POINTER, DIMENSION(:) :: XPSNV_A       ! Snow fraction over vegetation
REAL, POINTER, DIMENSION(:) :: XPSN
!
REAL, POINTER, DIMENSION(:)   :: XSNOWFREE_ALB     ! snow free albedo                        (-)
REAL, POINTER, DIMENSION(:)   :: XSNOWFREE_ALB_VEG ! snow free albedo for vegetation         (-)
REAL, POINTER, DIMENSION(:)   :: XSNOWFREE_ALB_SOIL! snow free albedo for soil
!
REAL, POINTER, DIMENSION(:) :: XVEG            ! vegetation cover fraction               (-)
!
REAL, POINTER, DIMENSION(:) :: XLAI          ! Leaf Area Index                         (m2/m2)
!
REAL, POINTER, DIMENSION(:) :: XEMIS         ! surface emissivity                      (-)
REAL, POINTER, DIMENSION(:) :: XZ0           ! surface roughness length                (m)
!
REAL, POINTER, DIMENSION(:) :: XRSMIN        ! minimum stomatal resistance             (s/m)
REAL, POINTER, DIMENSION(:) :: XGAMMA        ! coefficient for the calculation
!                                              ! of the surface stomatal
!                                              ! resistance
REAL, POINTER, DIMENSION(:) :: XWRMAX_CF     ! coefficient for maximum water 
!                                              ! interception 
!                                              ! storage capacity on the vegetation      (-)
REAL, POINTER, DIMENSION(:) :: XRGL          ! maximum solar radiation
!                                              ! usable in photosynthesis      
REAL, POINTER, DIMENSION(:) :: XCV           ! vegetation thermal inertia coefficient  (K m2/J)
REAL, POINTER, DIMENSION(:)    :: XLAIMIN    ! minimum LAI (Leaf Area Index)           (m2/m2)
REAL, POINTER, DIMENSION(:)    :: XSEFOLD    ! e-folding time for senescence           (s)
REAL, POINTER, DIMENSION(:)    :: XGMES      ! mesophyll conductance                   (m s-1)
REAL, POINTER, DIMENSION(:)    :: XGC        ! cuticular conductance                   (m s-1)
REAL, POINTER, DIMENSION(:)    :: XF2I       ! critical normilized soil water 
!                                              ! content for stress parameterisation
REAL, POINTER, DIMENSION(:)    :: XBSLAI     ! ratio d(biomass)/d(lai)                 (kg/m2)
!
REAL, POINTER, DIMENSION(:)    :: XCE_NITRO  ! leaf aera ratio sensitivity to 
!                                            ! nitrogen concentration                (m2/kg)
REAL, POINTER, DIMENSION(:)    :: XCF_NITRO  ! lethal minimum value of leaf area
!                                              ! ratio                                 (m2/kg)
REAL, POINTER, DIMENSION(:)    :: XCNA_NITRO ! nitrogen concentration of active 
!                                              ! biomass       
LOGICAL, POINTER, DIMENSION(:) :: LSTRESS    ! vegetation response type to water
!                                              ! stress (true:defensive false:offensive) (-)
!
REAL, POINTER, DIMENSION(:) :: XALBNIR_VEG       ! vegetation near-infra-red albedo        (-)
REAL, POINTER, DIMENSION(:) :: XALBVIS_VEG       ! vegetation visible albedo               (-)
REAL, POINTER, DIMENSION(:) :: XALBUV_VEG        ! vegetation UV albedo                    (-)
!
REAL, POINTER, DIMENSION(:) :: XALBNIR       ! near-infra-red albedo                   (-)
REAL, POINTER, DIMENSION(:) :: XALBVIS       ! visible albedo                          (-)
REAL, POINTER, DIMENSION(:) :: XALBUV        ! UV albedo
!
REAL, POINTER, DIMENSION(:) :: XGNDLITTER        ! ground litter fraction                  (-)
REAL, POINTER, DIMENSION(:) :: XH_VEG            ! height of vegetation                           (m)
REAL, POINTER, DIMENSION(:) :: XZ0LITTER         ! ground litter roughness length                 (m)
!
REAL, POINTER, DIMENSION(:) :: XALBNIR_SOIL      ! soil near-infra-red albedo              (-)
REAL, POINTER, DIMENSION(:) :: XALBVIS_SOIL      ! soil visible albedo                     (-)
REAL, POINTER, DIMENSION(:) :: XALBUV_SOIL       ! soil UV albedo
!
TYPE (DATE_TIME), POINTER, DIMENSION(:)  :: TSEED          ! date of seeding
TYPE (DATE_TIME), POINTER, DIMENSION(:)  :: TREAP          ! date of reaping
REAL, POINTER, DIMENSION(:)         :: XWATSUP        ! water supply during irrigation process (mm)
REAL, POINTER, DIMENSION(:)         :: XIRRIG         ! flag for irrigation (irrigation if >0.)
!
!
END TYPE ISBA_PE_t
!
TYPE ISBA_NK_t
!
TYPE(ISBA_K_t), DIMENSION(:), POINTER :: AL=>NULL()
!
END TYPE ISBA_NK_t
!
TYPE ISBA_NP_t
!
TYPE(ISBA_P_t), DIMENSION(:), POINTER :: AL=>NULL()
!
END TYPE ISBA_NP_t
!
TYPE ISBA_NPE_t
!
TYPE(ISBA_PE_t), DIMENSION(:), POINTER :: AL=>NULL()
!
END TYPE ISBA_NPE_t
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
CONTAINS
!
SUBROUTINE ISBA_S_INIT(YISBA_S)
TYPE(ISBA_S_t), INTENT(INOUT) :: YISBA_S
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_N:ISBA_S_INIT",0,ZHOOK_HANDLE)
!
NULLIFY(YISBA_S%XZS)
NULLIFY(YISBA_S%XCOVER)
NULLIFY(YISBA_S%LCOVER)
!
NULLIFY(YISBA_S%XTI_MIN)
NULLIFY(YISBA_S%XTI_MAX)
NULLIFY(YISBA_S%XTI_MEAN)
NULLIFY(YISBA_S%XTI_STD)
NULLIFY(YISBA_S%XTI_SKEW)
!
NULLIFY(YISBA_S%XSOC)
NULLIFY(YISBA_S%XPH)
NULLIFY(YISBA_S%XFERT)
!
NULLIFY(YISBA_S%XABC)
NULLIFY(YISBA_S%XPOI) 
!
NULLIFY(YISBA_S%XFRACSOC)
NULLIFY(YISBA_S%XTAB_FSAT)
NULLIFY(YISBA_S%XTAB_WTOP)
NULLIFY(YISBA_S%XTAB_QTOP)
NULLIFY(YISBA_S%XF_PARAM)
NULLIFY(YISBA_S%XC_DEPTH_RATIO)
NULLIFY(YISBA_S%XCPL_DRAIN)
NULLIFY(YISBA_S%XCPL_RUNOFF)
NULLIFY(YISBA_S%XCPL_ICEFLUX)
NULLIFY(YISBA_S%XCPL_EFLOOD)
NULLIFY(YISBA_S%XCPL_PFLOOD)
NULLIFY(YISBA_S%XCPL_IFLOOD)
NULLIFY(YISBA_S%XPERTVEG)
NULLIFY(YISBA_S%XPERTLAI)
NULLIFY(YISBA_S%XPERTCV)
NULLIFY(YISBA_S%XPERTALB)
NULLIFY(YISBA_S%XPERTZ0)
NULLIFY(YISBA_S%XEMIS_NAT)
!
NULLIFY(YISBA_S%XTSRAD_NAT)
!
NULLIFY(YISBA_S%XINNOV)
NULLIFY(YISBA_S%XRESID)
!
NULLIFY(YISBA_S%XWORK_WR)
!
NULLIFY(YISBA_S%XWSN_WR)
NULLIFY(YISBA_S%XRHO_WR)
NULLIFY(YISBA_S%XALB_WR)
NULLIFY(YISBA_S%XHEA_WR)
NULLIFY(YISBA_S%XAGE_WR)
NULLIFY(YISBA_S%XSG1_WR)
NULLIFY(YISBA_S%XSG2_WR)
NULLIFY(YISBA_S%XHIS_WR)
!
NULLIFY(YISBA_S%TDATE_WR)
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_N:ISBA_S_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_S_INIT
!
SUBROUTINE ISBA_K_INIT(YISBA_K)
TYPE(ISBA_K_t), INTENT(INOUT) :: YISBA_K
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_N:ISBA_K_INIT",0,ZHOOK_HANDLE)
!
NULLIFY(YISBA_K%XSAND)
NULLIFY(YISBA_K%XCLAY)
NULLIFY(YISBA_K%XPERM)
NULLIFY(YISBA_K%XRUNOFFB)
NULLIFY(YISBA_K%XWDRAIN)
!
NULLIFY(YISBA_K%XTDEEP)  
NULLIFY(YISBA_K%XGAMMAT)
NULLIFY(YISBA_K%XMPOTSAT)
NULLIFY(YISBA_K%XBCOEF)
NULLIFY(YISBA_K%XWWILT)
NULLIFY(YISBA_K%XWFC)  
NULLIFY(YISBA_K%XWSAT)
NULLIFY(YISBA_K%XCGSAT)
NULLIFY(YISBA_K%XC4B)
NULLIFY(YISBA_K%XACOEF)
NULLIFY(YISBA_K%XPCOEF)  
NULLIFY(YISBA_K%XHCAPSOIL)
NULLIFY(YISBA_K%XCONDDRY)
NULLIFY(YISBA_K%XCONDSLD)
NULLIFY(YISBA_K%XFWTD)
NULLIFY(YISBA_K%XWTD)   
NULLIFY(YISBA_K%XALBNIR_DRY)
NULLIFY(YISBA_K%XALBVIS_DRY)
NULLIFY(YISBA_K%XALBUV_DRY)
NULLIFY(YISBA_K%XALBNIR_WET)
NULLIFY(YISBA_K%XALBVIS_WET)
NULLIFY(YISBA_K%XALBUV_WET)
NULLIFY(YISBA_K%XWD0)
NULLIFY(YISBA_K%XKANISO)
NULLIFY(YISBA_K%XMUF)
NULLIFY(YISBA_K%XFSAT)
NULLIFY(YISBA_K%XFFLOOD)
NULLIFY(YISBA_K%XPIFLOOD)  
NULLIFY(YISBA_K%XFF)
NULLIFY(YISBA_K%XFFG)
NULLIFY(YISBA_K%XFFV)
NULLIFY(YISBA_K%XFFROZEN)
NULLIFY(YISBA_K%XALBF)
NULLIFY(YISBA_K%XEMISF)
NULLIFY(YISBA_K%XDIR_ALB_WITH_SNOW)
NULLIFY(YISBA_K%XSCA_ALB_WITH_SNOW)
!
NULLIFY(YISBA_K%XVEGTYPE)
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_N:ISBA_K_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_K_INIT
!
SUBROUTINE ISBA_P_INIT(YISBA_P)
TYPE(ISBA_P_t), INTENT(INOUT) :: YISBA_P
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_N:ISBA_P_INIT",0,ZHOOK_HANDLE)
!
YISBA_P%NSIZE_P = 0
NULLIFY(YISBA_P%XPATCH)
NULLIFY(YISBA_P%XVEGTYPE_PATCH)
NULLIFY(YISBA_P%NR_P)
NULLIFY(YISBA_P%XPATCH_OLD)  
NULLIFY(YISBA_P%XANMAX)
NULLIFY(YISBA_P%XFZERO)
NULLIFY(YISBA_P%XEPSO)
NULLIFY(YISBA_P%XGAMM)
NULLIFY(YISBA_P%XQDGAMM)
NULLIFY(YISBA_P%XQDGMES)
NULLIFY(YISBA_P%XT1GMES)
NULLIFY(YISBA_P%XT2GMES)
NULLIFY(YISBA_P%XAMAX)
NULLIFY(YISBA_P%XQDAMAX)
NULLIFY(YISBA_P%XT1AMAX)
NULLIFY(YISBA_P%XT2AMAX)
NULLIFY(YISBA_P%XAH)
NULLIFY(YISBA_P%XBH)
NULLIFY(YISBA_P%XTAU_WOOD)
NULLIFY(YISBA_P%XINCREASE)
NULLIFY(YISBA_P%XTURNOVER) 
NULLIFY(YISBA_P%XCONDSAT)
NULLIFY(YISBA_P%XTAUICE)
NULLIFY(YISBA_P%XC1SAT)
NULLIFY(YISBA_P%XC2REF)
NULLIFY(YISBA_P%XC3)
NULLIFY(YISBA_P%XC4REF)
NULLIFY(YISBA_P%XCPS)
NULLIFY(YISBA_P%XLVTT)
NULLIFY(YISBA_P%XLSTT)
NULLIFY(YISBA_P%XRUNOFFD)
NULLIFY(YISBA_P%XDZG)
NULLIFY(YISBA_P%XDZDIF)
NULLIFY(YISBA_P%XSOILWGHT)
NULLIFY(YISBA_P%XKSAT_ICE)
NULLIFY(YISBA_P%XBSLAI_NITRO)
NULLIFY(YISBA_P%XTOPQS)
!
NULLIFY(YISBA_P%XDG)
NULLIFY(YISBA_P%XDG_OLD)
NULLIFY(YISBA_P%NWG_LAYER)
NULLIFY(YISBA_P%XDROOT)
NULLIFY(YISBA_P%XDG2)
NULLIFY(YISBA_P%XROOTFRAC)
NULLIFY(YISBA_P%XD_ICE)
NULLIFY(YISBA_P%XH_TREE)
NULLIFY(YISBA_P%XZ0_O_Z0H)
NULLIFY(YISBA_P%XRE25)  
NULLIFY(YISBA_P%XDMAX)
!
NULLIFY(YISBA_P%XRED_NOISE)
NULLIFY(YISBA_P%XINCR)
NULLIFY(YISBA_P%XHO)
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_N:ISBA_P_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_P_INIT
!
SUBROUTINE ISBA_PE_INIT(YISBA_PE)
TYPE(ISBA_PE_t), INTENT(INOUT) :: YISBA_PE
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_N:ISBA_PE_INIT",0,ZHOOK_HANDLE)
!
NULLIFY(YISBA_PE%XLAI)  
NULLIFY(YISBA_PE%XVEG)
NULLIFY(YISBA_PE%XEMIS)
NULLIFY(YISBA_PE%XZ0)
NULLIFY(YISBA_PE%XRSMIN)
NULLIFY(YISBA_PE%XGAMMA)
NULLIFY(YISBA_PE%XWRMAX_CF)
NULLIFY(YISBA_PE%XRGL)  
NULLIFY(YISBA_PE%XCV)
NULLIFY(YISBA_PE%XLAIMIN)
NULLIFY(YISBA_PE%XSEFOLD)
NULLIFY(YISBA_PE%XGMES)
NULLIFY(YISBA_PE%XGC)
NULLIFY(YISBA_PE%XF2I)
NULLIFY(YISBA_PE%XBSLAI)
NULLIFY(YISBA_PE%XCE_NITRO)
NULLIFY(YISBA_PE%XCF_NITRO)
NULLIFY(YISBA_PE%XCNA_NITRO)
NULLIFY(YISBA_PE%LSTRESS)
NULLIFY(YISBA_PE%XALBNIR_VEG)
NULLIFY(YISBA_PE%XALBVIS_VEG)
NULLIFY(YISBA_PE%XALBUV_VEG)
NULLIFY(YISBA_PE%XALBNIR)
NULLIFY(YISBA_PE%XALBVIS)
NULLIFY(YISBA_PE%XALBUV)
!
NULLIFY(YISBA_PE%XGNDLITTER)
NULLIFY(YISBA_PE%XH_VEG)
NULLIFY(YISBA_PE%XZ0LITTER)
!
NULLIFY(YISBA_PE%XALBNIR_SOIL)
NULLIFY(YISBA_PE%XALBVIS_SOIL)
NULLIFY(YISBA_PE%XALBUV_SOIL)
!
NULLIFY(YISBA_PE%XWATSUP)
NULLIFY(YISBA_PE%XIRRIG)
!
NULLIFY(YISBA_PE%XWG)
NULLIFY(YISBA_PE%XWGI)
NULLIFY(YISBA_PE%XWR)
NULLIFY(YISBA_PE%XTG)
NULLIFY(YISBA_PE%XICE_STO)
NULLIFY(YISBA_PE%XWRL)
NULLIFY(YISBA_PE%XWRLI)
NULLIFY(YISBA_PE%XWRVN)
NULLIFY(YISBA_PE%XTV)
NULLIFY(YISBA_PE%XTL)
NULLIFY(YISBA_PE%XTC)
NULLIFY(YISBA_PE%XQC)
NULLIFY(YISBA_PE%XRESA)
NULLIFY(YISBA_PE%XAN)
NULLIFY(YISBA_PE%XANDAY)
NULLIFY(YISBA_PE%XANFM)
NULLIFY(YISBA_PE%XLE)
NULLIFY(YISBA_PE%XFAPARC)
NULLIFY(YISBA_PE%XFAPIRC)
NULLIFY(YISBA_PE%XLAI_EFFC)  
NULLIFY(YISBA_PE%XMUS)   
NULLIFY(YISBA_PE%XRESP_BIOMASS)
NULLIFY(YISBA_PE%XBIOMASS)
NULLIFY(YISBA_PE%XLITTER)
NULLIFY(YISBA_PE%XSOILCARB)
NULLIFY(YISBA_PE%XLIGNIN_STRUC)
NULLIFY(YISBA_PE%XPSNG)
NULLIFY(YISBA_PE%XPSNV)
NULLIFY(YISBA_PE%XPSNV_A)
NULLIFY(YISBA_PE%XSNOWFREE_ALB)
NULLIFY(YISBA_PE%XSNOWFREE_ALB_VEG)
NULLIFY(YISBA_PE%XSNOWFREE_ALB_SOIL)
NULLIFY(YISBA_PE%XPSN)
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_N:ISBA_PE_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_PE_INIT
!
SUBROUTINE ISBA_NK_INIT(YISBA_NK,KPATCH)
TYPE(ISBA_NK_t), INTENT(INOUT) :: YISBA_NK
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_N:ISBA_NK_INIT",0,ZHOOK_HANDLE)
!
IF (ASSOCIATED(YISBA_NK%AL)) THEN
  DO JP = 1,KPATCH
    CALL ISBA_K_INIT(YISBA_NK%AL(JP))
  ENDDO
  DEALLOCATE(YISBA_NK%AL)
ELSE
  ALLOCATE(YISBA_NK%AL(KPATCH))
  DO JP = 1,KPATCH
    CALL ISBA_K_INIT(YISBA_NK%AL(JP))
  ENDDO
ENDIF
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_N:ISBA_NK_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_NK_INIT
!
SUBROUTINE ISBA_NP_INIT(YISBA_NP,KPATCH)
TYPE(ISBA_NP_t), INTENT(INOUT) :: YISBA_NP
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_N:ISBA_NP_INIT",0,ZHOOK_HANDLE)
!
IF (ASSOCIATED(YISBA_NP%AL)) THEN
  DO JP = 1,KPATCH
    CALL ISBA_P_INIT(YISBA_NP%AL(JP))
  ENDDO
  DEALLOCATE(YISBA_NP%AL)
ELSE
  ALLOCATE(YISBA_NP%AL(KPATCH))
  DO JP = 1,KPATCH
    CALL ISBA_P_INIT(YISBA_NP%AL(JP))
  ENDDO
ENDIF
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_N:ISBA_NP_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_NP_INIT
!
SUBROUTINE ISBA_NPE_INIT(YISBA_NPE,KPATCH)
TYPE(ISBA_NPE_t), INTENT(INOUT) :: YISBA_NPE
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_N:ISBA_NPE_INIT",0,ZHOOK_HANDLE)
!
IF (ASSOCIATED(YISBA_NPE%AL)) THEN
  DO JP = 1,KPATCH
    CALL ISBA_PE_INIT(YISBA_NPE%AL(JP))
  ENDDO
  DEALLOCATE(YISBA_NPE%AL)
ELSE
  ALLOCATE(YISBA_NPE%AL(KPATCH))
  DO JP = 1,KPATCH
    CALL ISBA_PE_INIT(YISBA_NPE%AL(JP))
  ENDDO
ENDIF
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_N:ISBA_NPE_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_NPE_INIT

END MODULE MODD_ISBA_n
