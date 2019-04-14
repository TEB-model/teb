!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE SOLAR_PANEL(PTSTEP, PTSUN, PRESIDENTIAL,                  &
                               PEMIT_LW_ROOF, PEMIT_LWDN_PANEL, PLW_RAD, &
                               PABS_SW_PANEL, PTA, PN_FLOOR, PFRAC_PANEL,&
                               PEMIS_PANEL, PALB_PANEL, PEFF_PANEL,      &
                               PABS_LW_PANEL, PH_PANEL, PRN_PANEL,       &
                               PTHER_PRODC_DAY,                          &
                               PTHER_PROD_PANEL, PPHOT_PROD_PANEL,       &
                               PPROD_PANEL,                              &
                               PTHER_PROD_BLD,   PPHOT_PROD_BLD,         &
                               PPROD_BLD                                 )
!   ##########################################################################
!
!!****  *SOLAR_PANEL*  
!!
!!    PURPOSE
!!    -------
!
! Computes the energy budget of the solar panels.
!         
!     
!!**  METHOD
!     ------
!
!
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
!!
!!      
!!    AUTHOR
!!    ------
!!
!!      V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/2013
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_CSTS,     ONLY : XSTEFAN
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
REAL,               INTENT(IN)  :: PTSTEP          ! time step  (s)
REAL, DIMENSION(:), INTENT(IN)  :: PTSUN           ! solar time (s since solar midnight)
REAL, DIMENSION(:), INTENT(IN)  :: PRESIDENTIAL    ! Buildings Residential use fraction        (-)
REAL, DIMENSION(:), INTENT(IN)  :: PEMIT_LW_ROOF   ! Upwards   LW flux from roof               (W/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PEMIT_LWDN_PANEL! Downwards LW flux from panel              (W/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PLW_RAD         ! Incoming Longwave radiation               (W/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PABS_SW_PANEL   ! Absorbed solar energy by the solar panel  (W/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PTA             ! Air temperature                           (K)
REAL, DIMENSION(:), INTENT(IN)  :: PN_FLOOR        ! number of floors                          (-)
REAL, DIMENSION(:), INTENT(IN)  :: PFRAC_PANEL     ! fraction   of solar panel on roofs        (-)
REAL, DIMENSION(:), INTENT(IN)  :: PEMIS_PANEL     ! emissivity of solar panel                 (-)
REAL, DIMENSION(:), INTENT(IN)  :: PALB_PANEL      ! albedo     of solar panel                 (-)
REAL, DIMENSION(:), INTENT(IN)  :: PEFF_PANEL      ! efficiency of solar panel                 (-)
REAL, DIMENSION(:), INTENT(OUT) :: PABS_LW_PANEL   ! Absorbed LW enerby by solar panel         (W/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PH_PANEL        ! Sensible heat released by the solar panel (W/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PRN_PANEL       ! Net radiation     of the solar panel      (W/m2)
REAL, DIMENSION(:), INTENT(INOUT)::PTHER_PRODC_DAY ! Present day integrated thermal production of energy    (J/m2 panel)
REAL, DIMENSION(:), INTENT(OUT) :: PTHER_PROD_PANEL! Thermal      Energy production of the solar panel      (W/m2 panel)
REAL, DIMENSION(:), INTENT(OUT) :: PPHOT_PROD_PANEL! Photovoltaic Energy production of the solar panel      (W/m2 panel)
REAL, DIMENSION(:), INTENT(OUT) :: PPROD_PANEL     ! Averaged     Energy production of the solar panel      (W/m2 panel)
REAL, DIMENSION(:), INTENT(OUT) :: PTHER_PROD_BLD  ! Thermal      Energy production of the solar panel      (W/m2 bld)
REAL, DIMENSION(:), INTENT(OUT) :: PPHOT_PROD_BLD  ! Photovoltaic Energy production of the solar panel      (W/m2 bld)
REAL, DIMENSION(:), INTENT(OUT) :: PPROD_BLD       ! Averaged     Energy production of the solar panel      (W/m2 bld)
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*      0.2    declarations of parameters
!
! coefficient to take into account irradiance to estimate surface temperature of panel
REAL, PARAMETER :: XKT    = 0.05   ! (Km2/W)
!
! coefficient to take into account optimum orientation of the solar panel
REAL, PARAMETER :: XFT    = 1.10   ! (-)
!
! optimum panel temperature for photovoltaic production
REAL, PARAMETER :: XT_OPT = 298.15 ! (K)
!
! coefficient to take into account decrease of production for warm temperatures
REAL, PARAMETER :: XT_LOSS= 0.005  ! (K-1)
!
! difference of temperature between cold and hot water
REAL, PARAMETER :: XWATER_DT = 45. ! (K)
!
! density of m2 of panel per m2 of floor necessary for hot water production
REAL, PARAMETER :: XTHER_FLOOR= 1./30. ! (-)
!
! Annual thermal production for an increase of 1K of the water temperature
REAL, PARAMETER :: XTHER_RATE= 0.5 * 1.16*32. ! (kWh/year/m2 panel)
!
!! Maximum irradiance above which there is no more heating necessary (target temperature reached quickly)
!REAL, PARAMETER :: XTHER_IRR = 500.           ! (W/m2)
! Efficiency for thermal panels
REAL, PARAMETER :: XTHER_EFF = 0.60           ! (-)
!
!
!
!
!*      0.3    declarations of local variables
!
REAL                       :: ZTHER_DAILY_TARGET ! daily target for thermal production         (J/m2)
REAL, DIMENSION(SIZE(PTA)) :: ZIRRADIANCE ! incoming solar radiation normal to the solar panel (W/m2)
REAL, DIMENSION(SIZE(PTA)) :: ZTS_PANEL   ! Surface temperature of the upwards face of the solar panel (K)
REAL, DIMENSION(SIZE(PTA)) :: ZLWU_PANEL  ! Upwards longwave radiation from the solar panel (W/m2)
REAL, DIMENSION(SIZE(PTA)) :: ZTHER_FRAC  ! Fraction of thermal      panels per m2 of roof  (m2 panel/m2 bld)
REAL, DIMENSION(SIZE(PTA)) :: ZPHOT_FRAC  ! Fraction of photovoltaic panels per m2 of roof  (m2 panel/m2 bld)
REAL, DIMENSION(SIZE(PTA)) :: ZTHER_PRODC_DAY ! guess of daily production by thermal panels (J/m2)
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SOLAR_PANEL',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!*      0.4    Default values for output variables
!              -----------------------------------
!
PABS_LW_PANEL   = XUNDEF ! Absorbed LW enerby by solar panel         (W/m2)
PH_PANEL        = XUNDEF ! Sensible heat released by the solar panel (W/m2)
PRN_PANEL       = XUNDEF ! Net radiation     of the solar panel      (W/m2)
PTHER_PROD_PANEL= XUNDEF ! Thermal      Energy production of the solar panel      (W/m2)
PPHOT_PROD_PANEL= XUNDEF ! Photovoltaic Energy production of the solar panel      (W/m2)
PPROD_PANEL     = XUNDEF ! Averaged     Energy production of the solar panel      (W/m2)
PTHER_PROD_BLD  = 0.     ! Thermal      Energy production of the solar panel      (W/m2)
PPHOT_PROD_BLD  = 0.     ! Photovoltaic Energy production of the solar panel      (W/m2)
PPROD_BLD       = 0.     ! Averaged     Energy production of the solar panel      (W/m2)
!
!-------------------------------------------------------------------------------
!
!*      1.     Initializations
!              ---------------
!
!*      1.1    Checks hour for daily thermal production reset
!              ----------------------------------------------
!
!* energy reset between midnight and 1AM.
WHERE (PTSUN(:)>=0. .AND. PTSUN(:)<=3600.)  PTHER_PRODC_DAY(:) = 0.
!
!
!*      1.2    Daily target production for thermal panels
!              ------------------------------------------
!
ZTHER_DAILY_TARGET = (2.*XTHER_RATE) * XWATER_DT * (1000. / 365. * 3600. ) ! (J/m2)
! the factor 2 is to remove the assumption of cloudy days in the annual mean production
!
!-------------------------------------------------------------------------------
!* Note that computations are done only where solar panels are present
WHERE (PFRAC_PANEL(:)>0.)
!-------------------------------------------------------------------------------
!
!*      2.     Irradiance on panel
!              -------------------
!
  ZIRRADIANCE(:) = XFT * PABS_SW_PANEL(:) / (1.-PALB_PANEL(:) ) 
!
!-------------------------------------------------------------------------------
!
!*      3.     Solar panel temperature
!              -----------------------
!
  ZTS_PANEL  (:) = PTA(:) + XKT * ZIRRADIANCE(:)
!
!-------------------------------------------------------------------------------
!
!*      4.     Upwards solar panel LW radiation
!              --------------------------------
!
  ZLWU_PANEL  (:) =       PEMIS_PANEL(:)  * XSTEFAN * ZTS_PANEL(:)**4 &
                    + (1.-PEMIS_PANEL(:)) * PLW_RAD(:)
!
!-------------------------------------------------------------------------------
!
!*      5.     Solar panel LW budget
!              ---------------------
!
  PABS_LW_PANEL(:)= PLW_RAD(:) + PEMIT_LW_ROOF(:) - PEMIT_LWDN_PANEL(:) - ZLWU_PANEL(:)
!
!-------------------------------------------------------------------------------
!
!*      6.     Solar panel Net radiation
!              -------------------------
!
  PRN_PANEL(:)    = PABS_SW_PANEL(:) + PABS_LW_PANEL(:)
!
!-------------------------------------------------------------------------------
!
!*      7.     Fraction of panel surface types
!              -------------------------------
!
!*      7.1    Panel dedicated to thermal production of hot water 
!              --------------------------------------------------
!
  ZTHER_FRAC(:) = MIN( XTHER_FLOOR * PN_FLOOR(:) * PRESIDENTIAL(:), PFRAC_PANEL(:) )   ! (m2 thermal       panel / m2 roof)
!
!*      7.2    Photovoltaic panel
!              ------------------
!
  ZPHOT_FRAC(:) =  PFRAC_PANEL(:) - ZTHER_FRAC(:)                      ! (m2 photovoltaic panel / m2 roof)
!
!-------------------------------------------------------------------------------
!
!*      8.     Thermal Production of hot water (W/m2 thermal panel)
!              -------------------------------
!
!*      8.1    Instantaneous production
!              ------------------------
!
  PTHER_PROD_PANEL(:)= XTHER_EFF * ZIRRADIANCE(:)                       ! (W/m2)
!  PTHER_PROD_PANEL(:) = XTHER_RATE * XWATER_DT * (1000. / 24. / 365.)
!
!*      8.2    Integrated daily production
!              ---------------------------
!
  ZTHER_PRODC_DAY(:) = PTHER_PRODC_DAY(:) + PTSTEP * PTHER_PROD_PANEL(:) ! (J/m2)
!
!*      8.3    Daily production limited by daily target
!              ----------------------------------------
!
  ZTHER_PRODC_DAY(:) = MIN( ZTHER_PRODC_DAY(:) , ZTHER_DAILY_TARGET )
!
!*      8.4     Instantaneous production taking into account target limit if reached
!               --------------------------------------------------------------------
!
  PTHER_PROD_PANEL(:)= ( ZTHER_PRODC_DAY(:) - PTHER_PRODC_DAY(:) ) / PTSTEP
!
!*      8.5    Updates daily production
!              ------------------------
!
  PTHER_PRODC_DAY(:) = ZTHER_PRODC_DAY(:)

!-------------------------------------------------------------------------------
!
!*      9.     Photovoltaic Production  (W/m2 photovoltaic panel)
!              -----------------------
!
  PPHOT_PROD_PANEL(:) = PEFF_PANEL(:) * ZIRRADIANCE(:) * MIN(1.,1.-XT_LOSS*(ZTS_PANEL(:)-XT_OPT))
!
!-------------------------------------------------------------------------------
!
!*     10.     Averaged Production  (W/m2 panel)
!              -------------------
!
  PPROD_PANEL(:) =  PTHER_PROD_PANEL(:) * (ZTHER_FRAC(:) / PFRAC_PANEL(:)) &
                  + PPHOT_PROD_PANEL(:) * (ZPHOT_FRAC(:) / PFRAC_PANEL(:))
!
!-------------------------------------------------------------------------------
!
!*     11.     Sensible heat flux  (W/m2 panel)
!              ------------------
!
  PH_PANEL(:) = PRN_PANEL(:) - PPROD_PANEL(:)
!
!-------------------------------------------------------------------------------
!
!*     12.     Productions per building (W/m2 bld)
!              ------------------------
!
  PTHER_PROD_BLD(:) = PTHER_PROD_PANEL(:) * ZTHER_FRAC(:)
  PPHOT_PROD_BLD(:) = PPHOT_PROD_PANEL(:) * ZPHOT_FRAC(:)
  PPROD_BLD     (:) = PTHER_PROD_BLD  (:) + PPHOT_PROD_PANEL(:)
!
!-------------------------------------------------------------------------------
END WHERE
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SOLAR_PANEL',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE SOLAR_PANEL
