!auto_modi:spll_vegetation_evol.D
MODULE MODI_VEGETATION_EVOL
INTERFACE
    SUBROUTINE VEGETATION_EVOL(IO, DTI, PK, PEK, OAGRIP, PTSTEP, KMONTH, KDAY, PTIME, &
                               PLAT, PRHOA, P_CO2, ISSK, PRESP_BIOMASS_INST, PSWDIR)  
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_ISBA_n, ONLY : ISBA_P_t, ISBA_PE_t
USE MODD_SSO_n, ONLY : SSO_t
USE MODD_TYPE_DATE_SURF
IMPLICIT NONE
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(ISBA_P_t), INTENT(INOUT) :: PK
TYPE(ISBA_PE_t), INTENT(INOUT) :: PEK
LOGICAL,              INTENT(IN)    :: OAGRIP  ! agricultural practices
REAL,                 INTENT(IN)    :: PTSTEP  ! time step
INTEGER,              INTENT(IN)    :: KMONTH  ! current month
INTEGER,              INTENT(IN)    :: KDAY    ! current day
REAL,                 INTENT(IN)    :: PTIME   ! current time since midnight
REAL,   DIMENSION(:), INTENT(IN)    :: PLAT    ! latitude of each grid point
REAL,   DIMENSION(:), INTENT(IN)    :: PRHOA   ! air density
REAL,   DIMENSION(:), INTENT(IN)    :: P_CO2 ! CO2 concentration [ppmm]
TYPE(SSO_t), INTENT(INOUT) :: ISSK
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRESP_BIOMASS_INST ! instantaneous respiration of biomass (kgCO2/kgair m/s)
REAL, DIMENSION(:),   INTENT(IN),   OPTIONAL :: PSWDIR    ! Global incoming shortwave radiation (W m-2)
END SUBROUTINE VEGETATION_EVOL
END INTERFACE
END MODULE MODI_VEGETATION_EVOL
