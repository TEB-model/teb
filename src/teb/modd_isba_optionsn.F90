!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!##################
MODULE MODD_ISBA_OPTIONS_n
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
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE ISBA_OPTIONS_t
!
! *  General PGD options
!
LOGICAL :: LECOCLIMAP         ! T: parameters computed from ecoclimap
!                             ! F: they are read in the file
!
LOGICAL :: LPAR               ! T: parameters computed from ecoclimap
!                             ! F: they are read in the file
!
INTEGER :: NPATCH             ! maximum number of sub-tiles (patches)
!                             ! used at any grid point within a 
!                             ! natural surface fraction
INTEGER :: NGROUND_LAYER      ! number of ground layers
!
 CHARACTER(LEN=3) :: CISBA    ! type of ISBA version:
!                             ! '2-L' (default)
!                             ! '3-L'
!                             ! 'DIF'
 CHARACTER(LEN=4) :: CPEDOTF  ! NOTE: Only used when HISBA = DIF
!                             ! 'CH78' = Clapp and Hornberger 1978 for BC (Default)
!                             ! 'CO84' = Cosby et al. 1988 for BC
 CHARACTER(LEN=3) :: CPHOTO   ! type of photosynthesis
!                             ! 'NON'
!                             ! 'AST'
!                             ! 'NIT'
!                             ! 'NCB'
!
REAL, POINTER, DIMENSION(:)  :: XSOILGRID   ! Soil layer grid as reference for DIF
!
LOGICAL :: LTR_ML             ! new radiative transfert
!
REAL :: XRM_PATCH             ! threshold to remove little fractions of patches 
!
LOGICAL :: LSOCP              ! Soil organic carbon profile data
LOGICAL :: LCTI               ! Topographic index data
LOGICAL :: LPERM              ! Permafrost distribution data
LOGICAL :: LNOF  
!
! Type of vegetation (simplification of vegetation charaterization)
!
CHARACTER(LEN=4)             :: CTYPE_HVEG   ! type of high vegetation
CHARACTER(LEN=4)             :: CTYPE_LVEG   ! type of low vegetation
CHARACTER(LEN=4)             :: CTYPE_NVEG   ! type of bare soil (no vegetation)
CHARACTER(LEN=5)             :: CTYP_COV     ! type of green roof
!
! * AGS and Carbon PGD options
!
INTEGER :: NNBIOMASS    ! number of biomass pools
INTEGER :: NNLITTER     ! number of litter pools
INTEGER :: NNLITTLEVS   ! number of litter levels
INTEGER :: NNSOILCARB   ! number of soil carbon pools  
!
! * PGD MEB OPTIONS
!
LOGICAL, POINTER, DIMENSION(:) :: LMEB_PATCH  ! Vector with T/F values
                                              ! True = treat patch with multi-energy balance 
                                              ! False = treat patch with classical ISBA 
LOGICAL                        :: LFORC_MEASURE ! True = Forcing data from observations
!                                               ! False = Forcing data from atmospheric model (default)
LOGICAL                        :: LMEB_LITTER ! Activate Litter
LOGICAL                        :: LMEB_GNDRES ! Activate Ground Resistance
!
! * General PREP options
!
LOGICAL :: LCANOPY ! T: SBL scheme within the canopy
!                  ! F: no atmospheric layers below forcing level
!
! * Carbon PREP options
 CHARACTER(LEN=3) :: CRESPSL  ! Soil respiration
!                             ! 'DEF' = Default: Norman (1992)
!                             ! 'PRM' = New Parameterization
!                             ! 'CNT' = CENTURY model (Gibelin 2007)
!
!
! * General MODEL options
!
 CHARACTER(LEN=4) :: CC1DRY   ! C1 formulation for dry soils
!                             ! 'DEF ' = DEFault: Giard-Bazile formulation
!                             ! 'GB93' = Giordani 1993, Braud 1993 
!                             ! discontinuous at WILT
 CHARACTER(LEN=4) :: CSCOND   ! Thermal conductivity
!                             ! 'DEF ' = DEFault: NP89 implicit method
!                             ! 'PL98' = Peters-Lidard et al. 1998 used
!                             ! for explicit computation of CG
 CHARACTER(LEN=3) :: CSOILFRZ ! soil freezing-physics option
!                             ! 'DEF' = Default (Boone et al. 2000; 
!                             !        Giard and Bazile 2000)
!                             ! 'LWT' = Phase changes as above,
!                             !         but relation between unfrozen 
!                             !         water and temperature considered
 CHARACTER(LEN=4) :: CDIFSFCOND ! Mulch effects
!                               ! 'MLCH' = include the insulating effect of
!                               ! leaf litter/mulch on the surf. thermal cond.
!                               ! 'DEF ' = no mulch effect
!                               NOTE: Only used when YISBA = DIF
 CHARACTER(LEN=3) :: CSNOWRES   ! Turbulent exchanges over snow
!                               ! 'DEF' = Default: Louis (ISBA)
!                               ! 'RIL' = Maximum Richardson number limit
!                               !         for stable conditions ISBA-SNOW3L
!                               !         turbulent exchange option
 CHARACTER(LEN=4) :: CALBEDO    ! albedo type
!                               ! 'DRY ' 
!                               ! 'EVOL' 
!                               ! 'WET ' 
!                               ! 'USER' 
!                            NOTE that when using the YISBA='DIF' multi-layer soil option,
!                            the 'LWT' method is used. It is only an option
!                            when using the force-restore soil method ('2-L' or '3-L')
 CHARACTER(LEN=3) :: CCPSURF    ! specific heat at surface
!                               ! 'DRY' = default value (dry Cp)
!                               ! 'HUM' = Cp as a fct of specific humidity  
!
REAL :: XOUT_TSTEP              ! ISBA output writing time step
REAL :: XTSTEP                  ! ISBA time step
REAL :: XCGMAX                  ! maximum soil heat capacity
REAL :: XCDRAG                  ! drag coefficient in canopy
!
LOGICAL :: LGLACIER             ! True = Over permanent snow and ice, 
!                                 initialise WGI=WSAT,
!                                 Hsnow>=10m and allow 0.8<SNOALB<0.85
                                ! False = No specific treatment
LOGICAL :: LCANOPY_DRAG         ! T: drag activated in SBL scheme within the canopy
!                               ! F: no drag activated in SBL atmospheric layers
LOGICAL :: LVEGUPD              ! True = update vegetation parameters every decade
LOGICAL :: LPERTSURF            ! True  = apply random perturbations for ensemble prediction
                                ! False = no random perturbation (default)
!
!
! * SGH model options
!
! - Adjustable physical parameters
!
INTEGER :: NLAYER_HORT
INTEGER :: NLAYER_DUN
!
! - Sub-grid hydrology and vertical hydrology
!                                                     
 CHARACTER(LEN=4) :: CRUNOFF ! surface runoff formulation
!                            ! 'WSAT'
!                            ! 'DT92'
!                            ! 'SGH ' Topmodel
!                                                     
 CHARACTER(LEN=3) :: CKSAT   ! ksat
!                            ! 'DEF' = default value 
!                            ! 'SGH' = profil exponentiel
!                                           
LOGICAL :: LSOC              ! soil organic carbon effect
!                            ! False = default value 
!                            ! True  = soil SOC profil
!
 CHARACTER(LEN=3) :: CRAIN   ! Rainfall spatial distribution
                             ! 'DEF' = No rainfall spatial distribution
                             ! 'SGH' = Rainfall exponential spatial distribution
                             ! 
!
 CHARACTER(LEN=3) :: CHORT   ! Horton runoff
                             ! 'DEF' = no Horton runoff
                             ! 'SGH' = Horton runoff
!
! * AGS and carbon options
!
LOGICAL :: LNITRO_DILU          ! nitrogen dilution fct of CO2 (Calvet et al. 2008)
                                ! False = keep vegetation parameters constant in time
!
LOGICAL  :: LSPINUPCARBS  ! T: do the soil carb spinup, F: no
LOGICAL  :: LSPINUPCARBW  ! T: do the wood carb spinup, F: no  
REAL :: XSPINMAXS         ! max number of times CARBON_SOIL subroutine is
                          ! called for each timestep in simulation during
                          ! acceleration procedure number                             
REAL :: XSPINMAXW         ! max number of times the wood is accelerated  
REAL :: XCO2_START        ! Pre-industrial CO2 concentration
REAL :: XCO2_END          ! Begin-transient CO2 concentration
INTEGER :: NNBYEARSPINS   ! nbr years needed to reaches soil equilibrium 
INTEGER :: NNBYEARSPINW   ! nbr years needed to reaches wood equilibrium
INTEGER :: NNBYEARSOLD    ! nbr years executed at curent time step
INTEGER :: NSPINS         ! number of times the soil is accelerated
INTEGER :: NSPINW         ! number of times the wood is accelerated
!
LOGICAL :: LAGRI_TO_GRASS ! During soil carbon spinup with ISBA-CC, 
                          ! grass parameters are attributed to all agricultural PFT
!
! * Snow model options
!                                           
LOGICAL :: LSNOWDRIFT, LSNOWDRIFT_SUBLIM ! Logicals for snowdrift and sublimation

LOGICAL :: LSNOW_ABS_ZENITH  ! if True modify solar absorption as a function of solar zenithal angle
                             ! (physically wrong but better results in polar regions when CSNOWRAD=B92)
! Scheme of snow metamorphism (Crocus)
 CHARACTER(3) :: CSNOWMETAMO ! B92 (historical version, Brun et al 92), C13, T07, F06 (see Carmagnola et al 2014)
!
! radiative transfer scheme in snow (Crocus)
 CHARACTER(3) :: CSNOWRAD    ! B92 (historical version, Brun et al 92), TAR, TA1, TA2 (see Libois et al 2013)
!
! * Other options
!
LOGICAL :: LFLOOD       ! Activation of the flooding scheme
LOGICAL :: LWTD         ! Activation of Water table depth coupling
LOGICAL :: LCPL_RRM     ! Activation of the coupling
!
LOGICAL :: LTEMP_ARP      ! True  = time-varying force-restore soil temperature (as in ARPEGE)
                          ! False = No time-varying force-restore soil temperature (Default
INTEGER :: NTEMPLAYER_ARP ! Number of force-restore soil temperature layer, including Ts (Default = 4)
                          ! Only used if LTEMP_ARP=True
REAL, POINTER, DIMENSION(:) ::  XSODELX       ! Pulsation for each layer (Only used if LTEMP_ARP=True)
!
END TYPE ISBA_OPTIONS_t
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
CONTAINS

SUBROUTINE ISBA_OPTIONS_INIT(IO)
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_OPTIONS_N:ISBA_OPTIONS_INIT",0,ZHOOK_HANDLE)
!
NULLIFY(IO%XSOILGRID)
NULLIFY(IO%LMEB_PATCH)
!  
IO%CTYPE_HVEG='    '
IO%CTYPE_LVEG='    '
IO%CTYPE_NVEG='    '
IO%CTYP_COV='     '
IO%LPAR=.FALSE.
IO%LECOCLIMAP=.FALSE.
IO%NPATCH=0
IO%NGROUND_LAYER=0
IO%CISBA=' '
IO%CPEDOTF=' '
IO%CPHOTO=' '
IO%LTR_ML=.FALSE.
IO%XRM_PATCH=0.0
IO%LSOCP=.FALSE.
IO%LCTI=.FALSE.
IO%LPERM=.FALSE.
IO%LNOF=.FALSE.
IO%NNBIOMASS=0
IO%NNLITTER=0
IO%NNLITTLEVS=0
IO%NNSOILCARB=0
IO%LFORC_MEASURE=.FALSE.
IO%LMEB_LITTER=.FALSE.
IO%LMEB_GNDRES=.FALSE.
!
IO%LCANOPY=.FALSE.
IO%CRESPSL=' '
!
IO%CC1DRY=' '
IO%CSCOND=' '
IO%CSOILFRZ=' '
IO%CDIFSFCOND=' '
IO%CSNOWRES=' '
IO%CALBEDO=' '
IO%CCPSURF=' '
IO%XOUT_TSTEP=0.
IO%XTSTEP=0.
IO%XCGMAX=0.
IO%XCDRAG=0.
IO%LGLACIER=.FALSE.
IO%LCANOPY_DRAG=.FALSE.
IO%LVEGUPD=.FALSE.
IO%LPERTSURF=.FALSE.
IO%NLAYER_HORT=0
IO%NLAYER_DUN=0
IO%CRUNOFF=' '
IO%CKSAT=' '
IO%CRAIN=' '
IO%CHORT=' '
IO%LSOC=.FALSE.
IO%LNITRO_DILU=.FALSE.
IO%LSPINUPCARBS=.FALSE.
IO%LSPINUPCARBW=.FALSE.
IO%XSPINMAXS=0.
IO%XSPINMAXW=0.
IO%XCO2_START=0.
IO%XCO2_END=0.
IO%NNBYEARSPINS=0
IO%NNBYEARSPINW=0
IO%NNBYEARSOLD=0
IO%NSPINS=1
IO%NSPINW=1
IO%LAGRI_TO_GRASS=.FALSE.
IO%LSNOWDRIFT=.TRUE.
IO%LSNOWDRIFT_SUBLIM=.FALSE.
IO%LSNOW_ABS_ZENITH=.FALSE.
IO%CSNOWMETAMO='B92'
IO%CSNOWRAD='B92'
IO%LFLOOD=.FALSE.
IO%LWTD=.FALSE.
IO%LCPL_RRM=.FALSE.
IO%LTEMP_ARP=.FALSE.
IO%NTEMPLAYER_ARP=0
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_OPTIONS_N:ISBA_OPTIONS_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_OPTIONS_INIT

END MODULE MODD_ISBA_OPTIONS_n
