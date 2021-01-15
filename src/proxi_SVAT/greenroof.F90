!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL version 2.1
!SFX_LIC version 1. See LICENSE, Licence_CeCILL_V2.1-en.txt and Licence_CeCILL_V2.1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE GREENROOF (DTCO, G, T, TOP, TIR, DTV, GB, DK, DEK, DMK, GRO, S, K, P, PEK,    &
                          HIMPLICIT_WIND, TPTIME, PTSUN, PPEW_A_COEF, PPEW_B_COEF,  &
                          PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,       &
                          PTSTEP, PZREF, PUREF, PALB_GD, PTA, PQA, PEXNS, PEXNA, PRHOA,      &
                          PCO2, PPS, PRR, PSR, PZENITH, PSW, PLW, PVMOD,            &
                          PALBNIR_TVEG, PALBVIS_TVEG, PALBNIR_TSOIL, PALBVIS_TSOIL, &                
                          PRN, PH, PLE, PGFLUX, PSFCO2, PEVAP, PUW, PRUNOFF, PDRAIN,&
                          PAC, PQSAT, PTSRAD, PAC_AGG, PHU_AGG, PDEEP_FLUX, PIRRIG )  
!   ##################################################################################
!
!!****  *GREENROOF*  
!!
!!    PURPOSE
!!    -------
!!
!!    call the vegetation scheme (ISBA) inside TEB for greenroofs
!!     
!!**  METHOD
!!     ------
!!    based on subroutine "garden" 
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      
!!    REFERENCE
!!    ---------
!!    Based on subroutine "garden"
!!      
!!    AUTHOR
!!    ------
!!
!!      C. de Munck & A. Lemonsu          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!     Original    09/2011 
!     C. de Munck   02/2013  irrigation (drip irrigation)
!     B. decharme 04/2013 : Variables required in TEB to allow coupling with AROME/ALADIN/ARPEGE
!                           phasing call isba
!                           calculation of vegetation CO2 flux
!                           dummy for water table / surface coupling
!!    P. Samuelsson  10/2014  Introduced dummy variables in call to ISBA for MEB
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_SFX_GRID_n, ONLY : GRID_t
USE MODD_SSO_n, ONLY : SSO_t, SSO_INIT
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_IRRIG_n, ONLY : TEB_IRRIG_t
!
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_t
!
USE MODD_DIAG_n, ONLY : DIAG_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_S_t, ISBA_K_t, ISBA_P_t, ISBA_PE_t
!
USE MODD_AGRI_n, ONLY : AGRI_t, AGRI_INIT
!
USE MODD_SURF_PAR,             ONLY: XUNDEF
USE MODD_TYPE_DATE_SURF,       ONLY: DATE_TIME
USE MODD_CSTS,                 ONLY: XCPD, XLVTT, XKARMAN
!
USE MODI_ISBA
USE MODI_VEGETATION_UPDATE
USE MODI_VEGETATION_EVOL
USE MODI_CARBON_EVOL
USE MODE_THERMOS
USE MODI_ROOF_IMPL_COEF
USE MODI_TEB_IRRIG
USE MODI_FLAG_TEB_VEG_n

USE PROXI_SVAT_CONSTS, ONLY: GREENROOF_BR
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    Declarations of arguments
!
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(GRID_t), INTENT(INOUT) :: G
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_IRRIG_t), INTENT(INOUT) :: TIR
!
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTV
TYPE(GR_BIOG_t), INTENT(INOUT) :: GB
!
TYPE(DIAG_t), INTENT(INOUT) :: DK
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DEK
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DMK
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: GRO
TYPE(ISBA_S_t), INTENT(INOUT) :: S
TYPE(ISBA_K_t), INTENT(INOUT) :: K
TYPE(ISBA_P_t), INTENT(INOUT) :: P
TYPE(ISBA_PE_t), INTENT(INOUT) :: PEK
!
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
!                                                     ! 'OLD' = direct
!                                                     ! 'NEW' = Taylor serie, order 1
TYPE(DATE_TIME)     , INTENT(IN)    :: TPTIME             ! current date and time from teb
REAL, DIMENSION(:)  , INTENT(IN)    :: PTSUN              ! solar time      (s from midnight)
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEW_A_COEF        ! implicit coefficients
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEW_B_COEF        ! for wind coupling
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEQ_A_COEF        ! implicit coefficients
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEQ_B_COEF        ! for humidity
REAL, DIMENSION(:)  , INTENT(IN)    :: PPET_A_COEF        ! implicit coefficients
REAL, DIMENSION(:)  , INTENT(IN)    :: PPET_B_COEF        ! for temperature
REAL                , INTENT(IN)    :: PTSTEP             ! time step
REAL, DIMENSION(:)  , INTENT(IN)    :: PZREF              ! height of the first atmospheric level
REAL, DIMENSION(:)  , INTENT(IN)    :: PUREF              ! reference height for the wind
REAL, DIMENSION(:)  , INTENT(IN)    :: PALB_GD            ! green areas albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PTA                ! temperature at first atm. level 
REAL, DIMENSION(:)  , INTENT(IN)    :: PQA                ! specific humidity at first atm. level
REAL, DIMENSION(:)  , INTENT(IN)    :: PPS                ! pressure at the surface
REAL, DIMENSION(:)  , INTENT(IN)    :: PEXNA              ! Exner function at first atm. level
REAL, DIMENSION(:)  , INTENT(IN)    :: PEXNS              ! surface Exner function
REAL, DIMENSION(:)  , INTENT(IN)    :: PRHOA              ! air density at the lowest level
REAL, DIMENSION(:)  , INTENT(IN)    :: PCO2               ! CO2 concentration in the air    (kg/m3)
REAL, DIMENSION(:)  , INTENT(IN)    :: PRR                ! rain rate
REAL, DIMENSION(:)  , INTENT(IN)    :: PSR                ! snow rate
REAL, DIMENSION(:)  , INTENT(IN)    :: PZENITH            ! solar zenithal angle
REAL, DIMENSION(:)  , INTENT(IN)    :: PSW                ! incoming total solar rad on an horizontal surface
REAL, DIMENSION(:)  , INTENT(IN)    :: PLW                ! atmospheric infrared radiation
REAL, DIMENSION(:)  , INTENT(IN)    :: PVMOD              ! module of horizontal wind near first atm. level
REAL, DIMENSION(:)  , INTENT(IN)    :: PALBNIR_TVEG       ! nearIR  veg tot albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PALBVIS_TVEG       ! visible veg tot albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PALBNIR_TSOIL      ! nearIR  soil tot albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PALBVIS_TSOIL      ! visible soil tot albedo
!
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN         ! net radiation over greenroofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH          ! sensible heat flux over greenroofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE         ! latent heat flux over greenroofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX      ! flux through the greenroofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PSFCO2      ! flux of greenroof CO2       (m/s*kg_CO2/kg_air)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PEVAP       ! total evaporation over greenroofs (kg/m2/s)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PUW         ! friction flux (m2/s2)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRUNOFF     ! greenroof surface runoff
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDRAIN      ! greenroof surface drainage
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC         ! greenroof aerodynamical conductance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PQSAT       ! saturation humidity
REAL, DIMENSION(:)  , INTENT(OUT)   :: PTSRAD      ! greenroof radiative surface temp. (snow free)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_AGG     ! aggreg. aeodynamic resistance for greenroofs for latent heat flux
REAL, DIMENSION(:)  , INTENT(OUT)   :: PHU_AGG     ! aggreg. relative humidity for greenroofs for latent heat flux
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDEEP_FLUX  ! Heat Flux at the bottom layer of the greenroof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PIRRIG      ! greenroof summer irrigation rate
!
!
!*      0.2    Declarations of local variables
!
!
!-------------------------------------------------------------------------------
!
!*      1.     Proxi model based on a fixed Bowen ratio
!              ----------------------------------------
!
!
PRN(:) = (1.-PALB_GD) * PSW(:)
!
PH (:) = GREENROOF_BR * PRN(:)
PLE(:) = (1 - GREENROOF_BR) * PRN(:)
!
!* Conduction heat flux is neglected
PGFLUX(:) = 0.
!
!* CO2 flux is neglected
PSFCO2(:) = 0.
!
!* evaporation
PEVAP(:) = PLE(:) / XLVTT
!
!* Friction flux: assumes neutral formulation with roughness length of 0.01m
PUW(:) = - (XKARMAN/LOG(PUREF(:)/0.01))**2 * PVMOD(:)**2
!
!* Aerodynamical conductance: neglected because used further only for
!  implicitation of canyon air temperature when the heat flux depends on the
!  surface temperature
!
PAC(:) = 0.
!
!* surface saturation humidity
PQSAT(:) = QSAT(PTA(:),PPS(:))
!
!* Surface temperature : set equal to air temperature
PTSRAD(:) = PTA(:)
!
!* aerocynamical conductance for latent heat and surface humidity
PAC_AGG(:) = 0.    ! neglected (latent flux does not depend on surface humidity)
PHU_AGG(:) = 0.3   ! surface humidity set to 30%
!
!* Heat Flux at the bottom layer of the greenroof
PDEEP_FLUX(:) = 0.
!
!* greenroof hydrological diagnostics
PRUNOFF(:) = 0.    ! greenroof surface runoff
PDRAIN (:) = 0.    ! greenroof total (vertical) drainage
PIRRIG (:) = 0.    ! greenroof irrigation during time step
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE GREENROOF
