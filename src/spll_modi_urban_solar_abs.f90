!auto_modi:spll_urban_solar_abs.D
MODULE MODI_URBAN_SOLAR_ABS
INTERFACE
    SUBROUTINE URBAN_SOLAR_ABS(TOP, T, B, DMT, PDIR_SW, PSCA_SW, PZENITH, PAZIM,   &
                               PFRAC_PANEL, PALB_PANEL, PALB_GD, PSVF_GD, PALB_GRF, &
                               PDN_RF, PDF_RF, PDN_RD, PDF_RD, PREC_SW_RD,          &
                               PREC_SW_SN_RD, PREC_SW_WL_A, PREC_SW_WL_B,           &
                               PREC_SW_GD, PREC_SW_RF, PDIR_ALB_TWN, PSCA_ALB_TWN,  &
                               PSW_RAD_GD, PREC_SW_WIN, PREF_SW_GRND, PREF_SW_FAC,  &
                               PE_SHADING, OSHAD_DAY, OSHADE, OALB_ONLY )  
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_t
IMPLICIT NONE
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(BEM_t), INTENT(INOUT) :: B
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: DMT
REAL, DIMENSION(:), INTENT(IN)    :: PDIR_SW           ! incoming direct solar radiation
REAL, DIMENSION(:), INTENT(IN)    :: PSCA_SW           ! scattered incoming solar rad.
REAL, DIMENSION(:), INTENT(IN)    :: PZENITH           ! solar zenithal angle
REAL, DIMENSION(:), INTENT(IN)    :: PAZIM             ! solar azimuthal angle
REAL, DIMENSION(:), INTENT(IN)    :: PFRAC_PANEL       ! Fraction of solar panel on roofs (-)
REAL, DIMENSION(:), INTENT(IN)    :: PALB_PANEL        ! Albedo     of solar panels (-)
REAL, DIMENSION(:), INTENT(IN)    :: PALB_GD       ! GD areas albedo
REAL, DIMENSION(:), INTENT(IN)    :: PSVF_GD       ! GD areas sky view factor
REAL, DIMENSION(:), INTENT(IN)    :: PALB_GRF      ! green roof albedo
REAL, DIMENSION(:), INTENT(IN)    :: PDN_RF          ! snow-covered roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDF_RF          ! snow-free    roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDN_RD          ! snow-covered road fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDF_RD          ! snow-free    road fraction
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_RD     ! solar radiation received
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_WL_A   ! solar radiation received
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_WL_B   ! solar radiation received
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_GD     ! solar radiation received
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_RF     ! solar radiation received
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_SN_RD  ! solar radiation received
REAL, DIMENSION(:), INTENT(OUT)    :: PDIR_ALB_TWN   ! town direct albedo
REAL, DIMENSION(:), INTENT(OUT)    :: PSCA_ALB_TWN   ! town diffuse albedo
REAL, DIMENSION(:), INTENT(OUT)    :: PSW_RAD_GD     ! solar radiation reaching GD areas
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_WIN      ! solar radiation received by windows
REAL, DIMENSION(:), INTENT(OUT)    :: PREF_SW_GRND     ! total solar radiation reflected by ground
REAL, DIMENSION(:), INTENT(OUT)    :: PREF_SW_FAC      ! total solar radiation reflected by wall
REAL, DIMENSION(:), INTENT(OUT)    :: PE_SHADING       ! Energy that is not reflected 
LOGICAL, DIMENSION(:),INTENT(INOUT):: OSHAD_DAY        ! has shading been necessary this day ?
LOGICAL, DIMENSION(:),INTENT(IN)   :: OSHADE           ! are building conditions favorable for 
LOGICAL, INTENT(IN), OPTIONAL :: OALB_ONLY
END SUBROUTINE URBAN_SOLAR_ABS
END INTERFACE
END MODULE MODI_URBAN_SOLAR_ABS
