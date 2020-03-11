!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL licence version 2.1.
! See the following links for details:
! https://cecill.info/licences/Licence_CeCILL_V2.1-en.txt
! https://cecill.info/licences/Licence_CeCILL_V2.1-fr.txt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     #########
    SUBROUTINE GREENROOF(HIMPLICIT_WIND, TPTIME, PTSUN, PPEW_A_COEF, PPEW_B_COEF,    &
                PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,                  &
                PTSTEP, PZREF, PUREF,                                                &
                PTA, PQA, PEXNS, PEXNA,PRHOA, PCO2, PPS, PRR, PSR, PZENITH,          &
                PSW,PLW, PVMOD,                                                      &
                PRN_GREENROOF,PH_GREENROOF,PLE_GREENROOF,PGFLUX_GREENROOF,           &
                PSFCO2,PEVAP_GREENROOF, PUW_GREENROOF,                               &
                PAC_GREENROOF,PQSAT_GREENROOF,PTS_GREENROOF,                         &
                PAC_AGG_GREENROOF, PHU_AGG_GREENROOF,PDEEP_FLUX,                     &
                PRUNOFF_GREENROOF, PDRAIN_GREENROOF, PIRRIG_GREENROOF                )  
!   ##################################################################################
!
!!****  *GREENROOF*  
!!
!!    PURPOSE
!!    -------
!
!!call  a proxi of green roof scheme inside TEB
!
!!========================================================================
!!========================================================================
!!========================================================================
!!
!! ==> YOU ARE (MORE THAN) WELCOME TO USE YOUR OWN GREEN ROOF SCHEME HERE
!!
!!========================================================================
!!========================================================================
!!========================================================================
!
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
!!	C. de Munck & A. Lemonsu          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!     Original    09/2011 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS, ONLY : XLVTT , &   ! Latent heat constant for evaporation
                      XKARMAN     ! Von Karman constant
USE MODE_THERMOS                  ! Function to compute humidity at saturation
USE MODD_TYPE_DATE_SURF,    ONLY: DATE_TIME
!
IMPLICIT NONE
!
!*      0.1    Declarations of arguments
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
REAL, DIMENSION(:)  , INTENT(IN)    :: PZREF              ! height of the first atmospheric level                                                !
REAL, DIMENSION(:)  , INTENT(IN)    :: PUREF              ! reference height for the wind
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

REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_GREENROOF         ! net radiation over greenroofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_GREENROOF          ! sensible heat flux over greenroofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_GREENROOF         ! latent heat flux over greenroofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX_GREENROOF      ! flux through the greenroofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PSFCO2                ! flux of greenroof CO2       (kg/m2/s)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PEVAP_GREENROOF       ! total evaporation over greenroofs (kg/m2/s)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PUW_GREENROOF         ! friction flux (m2/s2)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_GREENROOF         ! greenroof aerodynamical conductance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PQSAT_GREENROOF       ! saturation humidity
REAL, DIMENSION(:)  , INTENT(INOUT) :: PTS_GREENROOF         ! greenroof radiative surface temp. (snow free)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_AGG_GREENROOF     ! aggreg. aeodynamic resistance for greenroofs for latent heat flux
REAL, DIMENSION(:)  , INTENT(OUT)   :: PHU_AGG_GREENROOF     ! aggreg. relative humidity for greenroofs for latent heat flux
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDEEP_FLUX            ! Heat Flux at the bottom layer of the greenroof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRUNOFF_GREENROOF     ! greenroof surface runoff
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDRAIN_GREENROOF      ! greenroof total (vertical) drainage
REAL, DIMENSION(:)  , INTENT(OUT)   :: PIRRIG_GREENROOF      ! greenroof irrigation during time step
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
!* albedo fixed to 0.15
PRN_GREENROOF(:) = (1.-0.15) * PSW(:)
!
!* Bowen ratio fixed to 1.
PH_GREENROOF (:) = 0.5 * PRN_GREENROOF(:)
PLE_GREENROOF(:) = 0.5 * PRN_GREENROOF(:)
!
!* Conduction heat flux is neglected
PGFLUX_GREENROOF(:) = 0.
!
!* CO2 flux is neglected
PSFCO2(:) = 0.
!
!* evaporation
PEVAP_GREENROOF(:) = PLE_GREENROOF(:) / XLVTT
!
!* Friction flux: assumes neutral formulation with roughness length of 0.01m
PUW_GREENROOF(:) = - (XKARMAN/LOG(PUREF(:)/0.01))**2 * PVMOD(:)**2
!
!* Aerodynamical conductance: neglected because used further only for
!  implicitation of canyon air temperature when the heat flux depends on the
!  surface temperature
!
PAC_GREENROOF(:) = 0.
!
!* surface saturation humidity
PQSAT_GREENROOF(:) = QSAT(PTA(:),PPS(:))
!
!* Surface temperature : set equal to air temperature
PTS_GREENROOF(:) = PTA(:)
!
!* aerocynamical conductance for latent heat and surface humidity
PAC_AGG_GREENROOF(:) = 0.    ! neglected (latent flux does not depend on surface humidity)
PHU_AGG_GREENROOF(:) = 0.3   ! surface humidity set to 30%
!
!* Heat Flux at the bottom layer of the greenroof
PDEEP_FLUX(:) = 0.
!
!* greenroof hydrological diagnostics
PRUNOFF_GREENROOF(:) = 0.    ! greenroof surface runoff
PDRAIN_GREENROOF (:) = 0.    ! greenroof total (vertical) drainage
PIRRIG_GREENROOF (:) = 0.    ! greenroof irrigation during time step
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE GREENROOF
