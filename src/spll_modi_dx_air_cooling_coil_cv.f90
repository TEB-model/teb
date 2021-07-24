!auto_modi:spll_dx_air_cooling_coil_cv.D
MODULE MODI_DX_AIR_COOLING_COIL_CV
INTERFACE
SUBROUTINE DX_AIR_COOLING_COIL_CV(PT_CANYON, PQ_CANYON, PPS, PRHOA,    &
                             PT_IN, PQ_IN, PCOP_RAT, PCAP_SYS_RAT,     &
                             PT_ADP, PF_WATER_COND,                    &
                             PM_SYS, PH_BLD_COOL, PH_WASTE, PLE_WASTE, &
                             PCOP, PCAP_SYS,  PT_OUT, PQ_OUT,          &
                             PDX_POWER,  PT_BLD_COOL      )
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
END SUBROUTINE DX_AIR_COOLING_COIL_CV
END INTERFACE
END MODULE MODI_DX_AIR_COOLING_COIL_CV
