!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
SUBROUTINE DX_AIR_COOLING_COIL_CV(PT_CANYON, PQ_CANYON, PPS, PRHOA,    &
                             PT_IN, PQ_IN, PCOP_RAT, PCAP_SYS_RAT,     &
                             PT_ADP, PF_WATER_COND,                    &
                             PM_SYS, PH_BLD_COOL, PH_WASTE, PLE_WASTE, &
                             PCOP, PCAP_SYS,  PT_OUT, PQ_OUT,          &
                             PDX_POWER,  PT_BLD_COOL      )
!
USE MODE_THERMOS
USE MODE_PSYCHRO
!
USE MODD_CSTS, ONLY : XCPD
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
! Extracted from EP_Engineer_ref p. 518
!
REAL, INTENT(IN) :: PT_CANYON        ! Canyon air temperature [K]
REAL, INTENT(IN) :: PQ_CANYON        ! Canyon air humidity ratio [kg kg-1]
REAL, INTENT(IN) :: PPS              ! Canyon air pressure [Pa]
REAL, INTENT(IN) :: PRHOA            ! Canyon air density [kg m-3]
REAL, INTENT(IN) :: PT_IN            ! Actual inlet air temperature [K]
REAL, INTENT(IN) :: PQ_IN            ! Actual inlet air humidity ratio [kg kg-1]
REAL, INTENT(IN) :: PCOP_RAT         ! Rated COP
REAL, INTENT(IN) :: PCAP_SYS_RAT     ! Rated capacity [W]
REAL, INTENT(IN) :: PT_ADP           ! Apparatus dewpoint [K]
REAL, INTENT(IN) :: PF_WATER_COND    ! fraction of evaporation of the condenser
REAL, INTENT(INOUT) :: PM_SYS        ! HVAC air mass flow rate [kg s-1]
REAL, INTENT(INOUT) :: PH_BLD_COOL   ! Sensible cooling load
REAL, INTENT(OUT) :: PH_WASTE        ! Sensible heat rejected by the condenser [W]
REAL, INTENT(OUT) :: PLE_WASTE       ! Latent heat rejected by the condenser [W]
REAL, INTENT(OUT) :: PCOP            ! Actual COP
REAL, INTENT(OUT) :: PCAP_SYS        ! Actual capacity [W]
REAL, INTENT(OUT) :: PT_OUT          ! Actual outlet temperature [K]
REAL, INTENT(OUT) :: PQ_OUT          ! Actual outlet humidity ratio [kg kg-1]
REAL, INTENT(OUT) :: PDX_POWER       ! Electrical power consumed by the DX unit [W]
REAL, INTENT(OUT) :: PT_BLD_COOL     ! Total energy supplied by the DX unit [W]
!
REAL :: ZTWB_CANYON   ! Canyon air wet-bulb temperature [ K]
REAL :: ZCAPTEMP      ! Total cooling capacity modifier curve function of temperature
REAL :: ZTW_IN        ! Wet-bulb temperature of the air entering the coil [ K]
REAL :: ZT_COND       ! Dry-bulb or wet-bulb air temperature entering the condenser [K]
REAL :: ZEIRTEMP      ! Energy input ratio modifier curve function of temperature
REAL :: ZPLR          ! Part load ratio
REAL :: ZPARTLOADF    ! Part load fraction correlation
REAL :: ZSHR          ! Actual coil sensible heat rate
REAL :: ZH_ADP        ! Enthalpy of air at ADP conditions [J/kg]
REAL :: ZH_OUT        ! Enthalpy of air leaving the cooling coil [J/kg]
REAL :: ZH_IN         ! Enthalpy of air entering the cooling coil [J/kg]
! Performance curves coefficients
REAL :: ZA1
REAL :: ZB1
REAL :: ZC1
REAL :: ZD1
REAL :: ZE1
REAL :: ZF1
! Total cooling capacity modifier curve function of flow fraction (desactivated)
!REAL :: ZCAPFLOW      
!REAL :: ZA2
!REAL :: ZB2
!REAL :: ZC2
REAL :: ZA3
REAL :: ZB3
REAL :: ZC3
REAL :: ZD3
REAL :: ZE3
REAL :: ZF3
! Energy input ratio modifier curve function of flow fraction (desactivated)
!REAL :: ZEIRFLOW      
!REAL :: ZA4
!REAL :: ZB4
!REAL :: ZC4
REAL :: ZA5
REAL :: ZB5
REAL :: ZC5
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('DX_AIR_COOLING_COIL_CV',0,ZHOOK_HANDLE)
!
!*      A.     Definitions
!              -----------
!
!
!*      B.     Performance curves
!              ------------------
! Total cooling capacity modifier curve function of temperature
! Obtained from default vaules of EnergyPlus SingleSpeedDX
! DOE-2.1E, COOL-CAP-FT for PTAC w/ SI temps
!  FurnaceDX DX Coil Cap-FT, Minimum Value of x, 12.77778;  
!                            Maximum Value of x, 23.88889;
!                            Minimum Value of y, 23.88889;
!                            Maximum Value of y, 46.11111
ZA1 =  0.942587793   !- Coefficient1 Constant
ZB1 =  0.009543347   !- Coefficient2 x
ZC1 =  0.00068377    !- Coefficient3 x**2
ZD1 =  -0.011042676  !- Coefficient4 y
ZE1 =  0.000005249   !- Coefficient5 y**2
ZF1 =  -0.00000972   !- Coefficient6 x*y
! Total cooling capacity modifier curve function of flow fraction
! DOE-2.1E, RATED-CCAP-FFLOW for PTAC
!  FurnaceDX DX Coil Cap-FF, Minimum Value of x, 0.5;  
!                            Maximum Value of x, 1.5;
!ZA2 =  0.8           !- Coefficient1 Constant
!ZB2 =  0.2           !- Coefficient2 x
!ZC2 =  0             !- Coefficient3 x**2
! Energy input ratio modifier curve function of temperature
! DOE-2.1E, COOL-EIR-FT for PTAC w/ SI temps
!  FurnaceDX DX Coil EIR-FT,  Minimum Value of x, 12.77778;  
!                            Maximum Value of x, 23.88889;
!                            Minimum Value of y, 23.88889;
!                            Maximum Value of y, 46.11111
ZA3 = 0.342414409    !- Coefficient1 Constant
ZB3 = 0.034885008    !- Coefficient2 x
ZC3 = -0.0006237     !- Coefficient3 x**2
ZD3 = 0.004977216    !- Coefficient4 y
ZE3 = 0.000437951    !- Coefficient5 y**2
ZF3 = -0.000728028   !- Coefficient6 x*y
! Energy input ratio modifier curve function of flow fraction
! DOE-2.1E, RATED-CEIR-FFLOW for PTAC
!  FurnaceDX DX Coil EIR-FF, Minimum Value of x, 0.5;  
!                            Maximum Value of x, 1.5;
!ZA4 = 1.1552         !- Coefficient1 Constant
!ZB4 = -0.1808        !- Coefficient2 x
!ZC4 = 0.0256         !- Coefficient3 x**2
! Part load fraction correlation
! PLF = l.- Cd(1.-PLR) where Cd = 0.15
!  FurnaceDX DX Coil PLF, Minimum Value of x, 0.0;  
!                         Maximum Value of x, 1.0;
ZA5 = 0.85           !- Coefficient1 Constant
ZB5 = 0.15           !- Coefficient2 x
ZC5 = 0              !- Coefficient3 x**2
!
!
!*      C.     Total cooling capacity 
!              ----------------------
! 
IF (PM_SYS/PRHOA/PCAP_SYS_RAT < 0.00004027) THEN
     PM_SYS = 0.00004027*PCAP_SYS_RAT*PRHOA
!     PRINT*,'ERROR: HVAC supply air flow rate must be greater than 0.00004027m3/s/W'
ELSE IF (PM_SYS/PRHOA/PCAP_SYS_RAT > 0.00006041) THEN
     PM_SYS = 0.00006041*PCAP_SYS_RAT*PRHOA
!     PRINT*,'ERROR: HVAC supply air flow rate must be lower than 0.00006041m3/s/W'
END IF
!
! Wet-bulb temperature entering the cooling coil
ZTW_IN      = TWB_FROM_TPQ(PT_IN, PPS, PQ_IN)
!
! Dry-bulb temperature of the air entering an air-cooled condenser 
! or wet-bulb temp entering a water-cooled condenser
ZTWB_CANYON = TWB_FROM_TPQ(PT_CANYON, PPS, PQ_CANYON)
ZT_COND     = ZTWB_CANYON + (PT_CANYON - ZTWB_CANYON)*(1. - PF_WATER_COND)
!
! Total cooling capacity modifier curve function of temperature
ZCAPTEMP = ZA1 + (ZTW_IN -273.15) * (ZB1 + ZC1*(ZTW_IN -273.15)) &
               + (ZT_COND-273.15) * (ZD1 + ZE1*(ZT_COND-273.15)) &
               + (ZTW_IN -273.15) * (ZT_COND-273.15) * ZF1
!
! Total cooling capacity
PCAP_SYS = PCAP_SYS_RAT * ZCAPTEMP
!
!*      D. Coil outlet conditions
!          ----------------------
!
! Inlet air enthalpy
ZH_IN  = ENTH_FN_T_Q(PT_IN,PQ_IN)
! Apparatus dewpoint enthalpy
ZH_ADP = ENTH_FN_T_Q(PT_ADP,QSAT(PT_ADP,PPS))
!
! Cooling coil sensible heat rate
IF (ZH_IN - ZH_ADP < 10.) THEN
  !
  ZSHR = 1.
  !
  PT_OUT = PT_ADP
  PQ_OUT = PQ_IN
  ZH_OUT = ENTH_FN_T_Q(PT_OUT,PQ_OUT)
  PT_BLD_COOL = 0.0
  !
ELSE
  !
  ZSHR  = MIN(XCPD*(PT_IN - PT_ADP)/(ZH_IN - ZH_ADP), 1.) !
  !
  ! Thermal load limited by the system capacity
  IF ( PH_BLD_COOL > PCAP_SYS * ZSHR ) PH_BLD_COOL = PCAP_SYS * ZSHR
  !
  ! Outlet air temperature
  PT_OUT = PT_IN - PH_BLD_COOL / PM_SYS / XCPD
  IF (PT_OUT < PT_ADP) PT_OUT = PT_ADP
  !
  ! Outlet air enthalpy
  ZH_OUT = ZH_IN - XCPD * (PT_IN - PT_OUT) / ZSHR
  !
  ! Outlet air humidity ratio
  PQ_OUT = Q_FN_T_ENTH(PT_OUT, ZH_OUT)
  !
  ! Total thermal energy supplied by the cooling coil
  PT_BLD_COOL = PM_SYS*(ZH_IN - ZH_OUT)
  !
END IF
!
!
!*      D. HVAC efficiency and electrical power consumed by the DX unit
!          ------------------------------------------------------------
! 
! Energy input ratio modifier curve function of temperature
ZEIRTEMP = ZA3 + (ZTW_IN -273.15) * (ZB3 + ZC3*(ZTW_IN -273.15))  &
               + (ZT_COND-273.15) * (ZD3 + ZE3*(ZT_COND-273.15))  &
               + (ZTW_IN -273.15) * (ZT_COND-273.15) * ZF3
!
! HVAC coefficient of performance
PCOP = PCOP_RAT / ZEIRTEMP
!
! Part load fraction correlation
ZPLR = PH_BLD_COOL / (PCAP_SYS * ZSHR)
!
IF (ZPLR > 1.0) ZPLR = 1.0
!
ZPARTLOADF = ZA5 + ZB5 * ZPLR + ZC5 * ZPLR**2
IF ( ZPARTLOADF < MIN(0.7, ZPLR) ) ZPARTLOADF = MIN(0.7, ZPLR)
!
! Electrical power consumed by the DX unit
! (compressors plus outdoor condenser fans)
PDX_POWER = PCAP_SYS / PCOP * ZPLR / ZPARTLOADF
!
!
!*      E. Waste heat emissions
!          --------------------
!
! Total heat rejected by the condenser
PLE_WASTE = (PT_BLD_COOL + PDX_POWER) * PF_WATER_COND
PH_WASTE  = (PT_BLD_COOL + PDX_POWER) * (1. - PF_WATER_COND)
!
IF (LHOOK) CALL DR_HOOK('DX_AIR_COOLING_COIL_CV',1,ZHOOK_HANDLE)
!
END SUBROUTINE DX_AIR_COOLING_COIL_CV
