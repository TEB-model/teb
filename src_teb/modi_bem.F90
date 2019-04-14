!auto_modi:spll_bem.D
MODULE MODI_BEM
INTERFACE
        SUBROUTINE BEM(BOP, T, B, DMT, PTSTEP, PSUNTIME, KDAY, PPS, PRHOA, PT_CAN,  &
                       PQ_CAN, PU_CAN, PHU_BLD, PT_RAD_IND, PFLX_BLD_FL, PFLX_BLD_MA,&
                       PRADHT_IN, PRAD_RF_MA, PRAD_RF_FL, PRAD_WL_MA, PRAD_WL_FL,&
                       PRAD_WIN_MA, PRAD_WIN_FL, PCONV_RF_BLD, PCONV_WL_BLD,     &
                       PCONV_WIN_BLD, PLOAD_IN_FL, PLOAD_IN_MA                )
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_t
IMPLICIT NONE
TYPE(BEM_OPTIONS_t), INTENT(INOUT) :: BOP
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(BEM_t), INTENT(INOUT) :: B
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: DMT
REAL,                INTENT(IN)   :: PTSTEP        ! Time step
REAL, DIMENSION(:),  INTENT(IN)   :: PSUNTIME       ! current solar time since midnight (solar time, s)
INTEGER,             INTENT(IN)   :: KDAY          ! Simulation day
REAL, DIMENSION(:),   INTENT(IN)  :: PPS          ! Canyon air pressure [Pa]
REAL, DIMENSION(:),   INTENT(IN)  :: PRHOA        ! Air density at the lowest level [kg m-3]
REAL, DIMENSION(:),   INTENT(IN)  :: PT_CAN    ! Canyon air temperature [K]
REAL, DIMENSION(:),   INTENT(IN)  :: PQ_CAN    ! Canyon air specific humidity [kg kg-1]
REAL, DIMENSION(:),   INTENT(IN)  :: PU_CAN    ! Canyon wind speed (m s-1)
REAL, DIMENSION(:),   INTENT(OUT)  :: PHU_BLD       ! Indoor relative humidity 0 < (-) < 1
REAL, DIMENSION(:),   INTENT(IN)  :: PT_RAD_IND    ! Indoor mean radiant temperature [K]
REAL, DIMENSION(:),   INTENT(OUT)  :: PFLX_BLD_FL! Heat flux from indoor air to floor 
REAL, DIMENSION(:),   INTENT(OUT)  :: PFLX_BLD_MA ! Heat flux from indoor air to mass 
REAL, DIMENSION(:),   INTENT(IN) :: PRADHT_IN     ! Indoor radiant heat transfer coefficient
REAL, DIMENSION(:)  , INTENT(IN)  :: PRAD_RF_MA  ! Rad. fluxes between roof and mass
REAL, DIMENSION(:)  , INTENT(IN)  :: PRAD_RF_FL ! Rad. fluxes between roof and floor
REAL, DIMENSION(:)  , INTENT(IN)  :: PRAD_WL_MA  ! Rad. fluxes between wall and mass
REAL, DIMENSION(:)  , INTENT(IN)  :: PRAD_WL_FL ! Rad. fluxes between wall and floor
REAL, DIMENSION(:)  , INTENT(IN)  :: PRAD_WIN_MA   ! Rad. fluxes between wind. and mass
REAL, DIMENSION(:)  , INTENT(IN)  :: PRAD_WIN_FL  ! Rad. fluxes between wind. and floor
REAL, DIMENSION(:)  , INTENT(IN)  :: PCONV_RF_BLD  ! Conv. fluxes between roof and indoor air
REAL, DIMENSION(:)  , INTENT(IN)  :: PCONV_WL_BLD  ! Conv. fluxes between wall and indoor air
REAL, DIMENSION(:)  , INTENT(IN)  :: PCONV_WIN_BLD   ! Conv. fluxes between wind. and indoor air
REAL, DIMENSION(:)  , INTENT(IN)  :: PLOAD_IN_FL  ! solar + int heat gain on floor W/m2 [floor]
REAL, DIMENSION(:)  , INTENT(IN)  :: PLOAD_IN_MA   ! solar + int heat gain on floor W/m2 [mass]
END SUBROUTINE BEM
END INTERFACE
END MODULE MODI_BEM
