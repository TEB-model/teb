!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ####################
      MODULE MODD_CO2V_PAR
!     ####################
!
!!*****MODD_CO2V_PAR*  
!!
!!    PURPOSE
!!    -------
!!    Supply constants and parameter values as
!!    a function of vegetation class to describe
!!    the vgegtation cover and for the
!!    stomatal conductance model of Jacobs.
!!    (Calvet et al. 1997, Agrig. and For. Met.)
!!     
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!      
!!    REFERENCE
!!    ---------
!!      
!!    AUTHOR
!!    ------
!!      A. Boone           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/10/97 
!!      Modified    12/03/04  by P LeMoigne: XAMMIN in (kg m-2 s-1)
!!      Modified    12/03/04  by P LeMoigne: XFZERO taken from Calvet 98
!!      Modified    18/10/04  by P LeMoigne: XFZERO splitted into 2
!!      P Le Moigne 09/2005 AGS modifs of L. Jarlan
!       S. Lafont    03/2009 : change unit of EPSO GAMM ANMAX 
!!      A.L. Gibelin 04/2009 : correction of some AGS parameters and suppress
!!                             useless parameters
!!      A.L. Gibelin 04/2009 : add parameters for PHOTO='NCB'
!!      A.L. Gibelin 06/2009 : add parameters for RESPSL='CNT'
!!      B. Decharme  05/2012 : Add XCC_NITRO and XBIOMASST_LIM (optimization)
!!      R. Alkama    05/2012 : parameters for 19 vegtype rather than 12
!!Seferian & Delire  06/2015 : updating XAMAX peding on TRY database Kattge et al. GCB 2011 
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!                           *** First, CO2 parameters ***
!                               ---------------------
!
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE_OLD, NVEGTYPE_ECOSG
!
IMPLICIT NONE
!
REAL, PARAMETER                      :: XSPIN_CO2 = 0.1
!                                       fraction of the total spinup period used
!                                       to ramp up CO2 concentration from XCO2_START to XCO2_END
!
REAL, PARAMETER                      :: XMCO2 = 44.0E-3, XMC   = 12.0E-3
!                                       molecular mass of CO2, 
!                                       and C (Carbon), respectively
!                                       (used for conversions in kg) 
!
REAL, PARAMETER                      :: XPARCF = 0.48
!                                       coefficient: PAR fraction of incoming solar radiation
!
REAL, PARAMETER                      :: XRACCF = 1.00
!                                       Factor for aerodynamic resistance for CO2
!
REAL, PARAMETER                      :: XPCCO2 = 0.40
!                                       proportion of Carbon in dry plant biomass (kgC/kgDM)
!
REAL, PARAMETER                      :: XIAOPT = 500.,  XDSPOPT = 0.0
!                                       optimum/initial values for absorbed global 
!                                       radiation, saturation deficit, respectively.

REAL, PARAMETER                      :: XXGT = 0.5
!                                       spherical distribution of leaves
!
REAL, PARAMETER                      :: XDIFRACF = 0.25
!                                       used in computation of fraction of diffusion
!
REAL, PARAMETER                      :: XXBOMEGA = 0.9442719
!                                       Factor Assuming the albedo from simple diffusion 
!                                       of the leaf (PAR) or 'omega'=0.2
!
REAL, PARAMETER                      :: XRDCF = 1./9.
!                                       Dark respiration factor/coefficient
!
REAL, PARAMETER                      :: XCONDCTMIN = 0.0002
!                                       minimum canopy conductance
!                                       or (1/5000)
!
REAL, PARAMETER                      :: XCONDSTMIN = 1.0E-5
!                                       minimum stomatal conductance for CO2
!                                       or (1/100000)
!
REAL, PARAMETER                      :: XANFMINIT = 1.E-8
!                                       initial maximum leaf assimilation
!
REAL, PARAMETER                      :: XAIRTOH2O = 1.606
!                                       ratio XMD (AIR)/XMV (H2O)
!
!
REAL, PARAMETER                      :: XCO2TOH2O = 1.6
!                                       ratio of the binary diffusivities of CO2 and H2O in air
!
!
REAL, PARAMETER                      :: XAW = 4.7, XASW = 2.8, XBW = 7.0
!                                       coefficient of stress universal relationship
!                                       for Woody species 
!
REAL, PARAMETER                      :: XDMAXN = 3.0E-2, XDMAXX = 3.0E-1
!                                       minimum and maximum air deficit stress parameters
!                                       used for Herbaceous offensive and defensive strategy
!
!
!                      Parameter values [C3,C4] at 25 C:
!
REAL, PARAMETER, DIMENSION(2) :: XTOPT   = (/25.0, 35.0/)
!                                       optimum/initial temperatures for evaluating 
!                                       compensation points
!
REAL, PARAMETER, DIMENSION(2) :: XFZERO1  = (/0.85, 0.50/)   ! AGS LAI
REAL, PARAMETER, DIMENSION(2) :: XFZERO2  = (/0.95, 0.60/)   ! AST, LST, NIT, NCB 
!                                       CO2_atm/CO2_int with no photorespiration or saturation deficit, 
!                                       used only for crops & herbaceous plants       
REAL, PARAMETER, DIMENSION(2) :: XFZEROTROP  = (/0.74, 0.74/) ! Tropical forests with and without 
!                                       Carrer et al. radiative transfer                                                                   
!
REAL, PARAMETER, DIMENSION(2) :: XEPSO   = (/0.017E-6,0.014E-6/) 
!                                       maximum initial quantum use efficiency 
!                                       (kgCO2 J-1 PAR )
!
REAL, PARAMETER, DIMENSION(2) :: XGAMM   = (/45.0, 2.8/)
!                                       CO2 conpensation concentration (ppmv)
!
REAL, PARAMETER, DIMENSION(2) :: XQDGAMM = (/1.5, 1.5/)
!                                       Q10 function for CO2 conpensation concentration
!
REAL, PARAMETER, DIMENSION(2) :: XQDGMES = (/2.0, 2.0/)
!                                       Q10 function for mesophyll conductance 
!
REAL, PARAMETER, DIMENSION(2) :: XT1GMES = (/5.0, 13.0/)
!                                       reference temperature for computing compensation
!                                       concentration function for mesophyll conductance:
!                                       minimum temperature 
!
REAL, PARAMETER, DIMENSION(2) :: XT2GMES = (/36.0, 36.0/)
!                                       reference temperature for computing compensation
!                                       concentration function for mesophyll conductance: 
!                                       maximum temperature
!
REAL, PARAMETER, DIMENSION(2) :: XQDAMAX = (/2.0, 2.0/)
!                                       Q10 function for leaf photosynthetic capacity
!
REAL, PARAMETER, DIMENSION(2) :: XT1AMAX = (/8.0, 13.0/)
!                                       reference temperature for computing compensation
!                                       concentration function for leaf photosynthetic 
!                                       capacity: minimum temperature
!
REAL, PARAMETER, DIMENSION(2) :: XT2AMAX = (/38.0, 38.0/)
!                                       reference temperature for computing compensation
!                                       concentration function for leaf photosynthetic
!                                       capacity: maximum temperature
!
REAL, PARAMETER, DIMENSION(2) :: XAH     = (/2.381, 5.323/)
!
!
REAL, PARAMETER, DIMENSION(2) :: XBH     = (/-0.6103, -0.8929/)
!                                       coefficient of stress universal relationship
!                                       for herbaceous plants
!
!
! Nitrogen option (YPHOTO='NIT' OR YPHOTO='NCB') parameters:
!
REAL, PARAMETER               :: XRESPFACTOR_NIT = 1.16E-7
!                                maintenance respiration rate (1% per day)
!                                of structural biomass (Faurie, 1994) [s-1]
!
REAL, PARAMETER               :: XCA_NIT = 0.38 
!                                rate of nitrogen dilution of above-ground biomass at all [CO2] 
!                                (Calvet and Soussana 2001, Gibelin et al. 2006)
!
REAL, PARAMETER               :: XCC_NIT = 0.753846
!                                Proportion of active biomass for 1t ha-1
!                                of total above-ground biomass [-]
!
REAL, PARAMETER               :: XCC_NITRO = 0.31425531725
!                                coef c for biomass in kg/m2
!                                Old : XCC_NIT/10.**XCA_NIT in nitro_decline.F90
!
REAL, PARAMETER               :: XBIOMASST_LIM = 4.7540042445E-2
!                                threshold value for leaf biomass and total 
!                                above ground biomass in nitrogen dilution theory
!                                Old : XCC_NITRO**(1.0/XCA_NIT) in nitro_decline.F90
!
!
! Nitrogen option with allocation of carbon (YPHOTO='NCB') parameters:
!
REAL, PARAMETER               :: XCOEFF_MAINT_RESP_ZERO = 1.19E-4/86400.
!                                value of coefficient for maintenance respiration at 0C (g/g/s)
!
REAL, PARAMETER               :: XSLOPE_MAINT_RESP = 0.16
!                                slope for maintenance respiration for temperature dependance (1/C)
!
REAL, PARAMETER, DIMENSION(NVEGTYPE_OLD) :: XTAU_WOOD = &
        (/ 0., 0., 0., 40.*365.*86400., 50.*365.*86400., 30.*365.*86400., 0., 0., 0., 0., &
           0., 0., 40.*365.*86400., 40.*365.*86400., 50.*365.*86400., 40.*365.*86400.,    &
           50.*365.*86400., 0., 40.*365.*86400. /)
!                                 Residence time in woody pools (s) (YPHOTO='NCB')
!
!
! Soil carbon (YPHOTO='NCB' and YRESPSL='CNT') parameters:
!
REAL, DIMENSION(6)      :: XCN
!                                       Biomass Carbon/Nitrogen ratio
!
REAL, DIMENSION(6)      :: XLC
!                                       Biomass Lignin/Carbon ratio
!
REAL, DIMENSION(6,2)    :: XFRAC_LITTER
!                                       Fraction of biomass pools going into litter pools
!
REAL, DIMENSION(2)      :: XTAU_LITTER
!                                       Residence times in litter pools (s)
!
REAL, DIMENSION(2,3,6)  :: XFRAC_SOILCARB
!                                       Fraction of litter decomposition flux that goes into soil
!
REAL, DIMENSION(3)      :: XTAU_SOILCARB
!                                       Residence times in carbon pools (s)
!
! Radiative transfer parameters
!
! single scattering albedo
REAL, PARAMETER                      :: XSSA_SUP = 0.15 ! single scatering albedo (PAR) for upper layer 
REAL, PARAMETER                      :: XSSA_INF = 0.15 ! single scatering albedo (PAR) for lower layer
!
REAL, PARAMETER                      :: XSSA_SUP_PIR = 0.80 ! single scatering albedo (PAR) for upper layer 
REAL, PARAMETER                      :: XSSA_INF_PIR = 0.80 ! single scatering albedo (PAR) for lower layer

! upper layer (calibration should depend on vegetation type...)
REAL, PARAMETER                      :: XK_SUP   = 2.  ! k_sup = 2*tan(zc); where zc is the crown angle
REAL, PARAMETER                      :: XXSI_SUP = 2.  ! hight/diameter ratio of the crown
!
! lower layer (calibration should depend on vegetation type...)
REAL, PARAMETER                      :: XK_INF   = 1.  ! k_sup = 2*tan(zc); where zc is the crown angle
REAL, PARAMETER                      :: XXSI_INF = 7.  ! height/diameter ratio of the crown
!
REAL, PARAMETER, DIMENSION(NVEGTYPE_OLD) :: XLAI_SHADE = &       ! above this value shading is activated
           (/ 0., 0., 0., 2., 3., 6., 2., 2., 2., 2., 2., 2., 2., 2., 3., 2., 3., 2., 2. /)
REAL, PARAMETER, DIMENSION(NVEGTYPE_OLD) :: XXB_SUP = &       ! b_sup = 1/omega_sup(zs=0) -1 calibration avec omega_sup(zs=0)=0.5
           (/ 1., 1., 1., 1., 2., 3., 2., 2., 2., 2., 2., 2., 1., 1., 2., 1., 2., 2., 1. /)
REAL, PARAMETER, DIMENSION(NVEGTYPE_OLD) :: XXB_INF = &       ! b_sup = 1/omega_sup(zs=0) -1 calibration avec omega_sup(zs=0)=0.5
           (/ 1., 1., 1., 4., 2., 4., 1., 1.5, 1.5, 1., 1., 1., 4., 4., 2., 4., 2., 1., 4. /)
!
! (Calvet et al. 2008) coefs for ratio of biomass to LAI with representation of nitrogen dilution fct of CO2
!REAL, PARAMETER, DIMENSION(NVEGTYPE_OLD) :: XPARAM = &       
!(/ 0., 0., 0., 2.56, 1.81, 1.81, 1.48, 1.48, 1.48, 1.81, 1.81, 1.81, 2.56, 2.56, 1.81, 2.56, &
!  1.81, 1.81, 2.56 /)
!
! Kronecker flag for nitrogen dilution hypothesis (based on Yin et al., 2002)
REAL, PARAMETER, DIMENSION(NVEGTYPE_OLD) :: XDILUDEC = &
(/ 0., 0., 0., 1., 0., 0., 1., 1., 1., 1., 1., 0., 1., 0., 0., 1., 1., 0., 1. /)
!
! Maximum Leaf photosynthetic capacity (kgCO2 m-2 s-1)
! Modified according to Kattge et al., 2009 median Vcmax at 25C values, except for TRBE : median-std
! For C3 PFTs : Ammax = Vcmax / 2.     (Jacobs, p 150)
! For C4 PFTs : Ammax = Vcmax
! Units : [Vcmax]=micromols_CO2 m-2 s-1, [Ammax]=kgCO2 m-2 s-1 --> [Ammax] = [Vcmax] * 44e-3 * 1e-6
REAL, PARAMETER, DIMENSION(NVEGTYPE_OLD) :: XAMAX   = &
(/ 1., 1., 1., 1.3E-6, 1.4E-6, 0.484E-6, 2.2E-6, 1.7E-6, 1.7E-6, 1.7E-6, 1.7E-6, 1.7E-6, &
   0.9E-6, 1.3E-6, 1.4E-6, 1.3E-6, 0.9E-6, 1.7E-6, 1.2E-6/)
!                           
INTEGER, PARAMETER, DIMENSION(NVEGTYPE_ECOSG) :: ITRANSFERT_ESG = &
        (/1,2,3,16,4,13,14,6,5,15,17,19,18,10,11,7,7,8,4,12/)

END MODULE MODD_CO2V_PAR






