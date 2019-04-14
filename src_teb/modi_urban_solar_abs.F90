!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_URBAN_SOLAR_ABS
INTERFACE
    SUBROUTINE URBAN_SOLAR_ABS(HBEM, HROAD_DIR, HWALL_OPT,                     &
                                 PDIR_SW, PSCA_SW, PZENITH, PAZIM,             &
                                 PBLD, PGARDEN, PROAD_DIR, PROAD, PFRAC_GR,    &
                                 PWALL_O_HOR, PCAN_HW_RATIO,                   &
                                 PALB_ROOF,                                    &
                                 PALB_ROAD, PSVF_ROAD, PALB_WALL, PSVF_WALL,   &
                                 PFRAC_PANEL, PALB_PANEL,                      &
                                 PALB_GARDEN, PSVF_GARDEN,                     &
                                 PALB_GREENROOF,                               &
                                 PASNOW_ROOF, PASNOW_ROAD,                     &
                                 PDN_ROOF, PDF_ROOF, PDN_ROAD, PDF_ROAD,       &
                                 PGR, PABS_WIN, PSHGC, PSHGC_SH, PALB_WIN,     &
                                 PABS_SW_ROOF, PABS_SW_ROAD,                   &
                                 PABS_SW_WALL_A, PABS_SW_WALL_B,               &
                                 PABS_SW_GARDEN, PABS_SW_GREENROOF,            &
                                 PABS_SW_SNOW_ROOF, PABS_SW_SNOW_ROAD,         &
                                 PABS_SW_SOLAR_PANEL,                          &
                                 PREC_SW_ROAD,  PREC_SW_SNOW_ROAD,             &
                                 PREC_SW_WALL_A, PREC_SW_WALL_B,               &
                                 PREC_SW_GARDEN, PREC_SW_ROOF,                 &
                                 PDIR_ALB_TOWN, PSCA_ALB_TOWN,                 &
                                 PSW_RAD_GARDEN, PABS_SW_WIN, PREC_SW_WIN,     &
                                 PTRAN_WIN,                                    &
                                 PREF_SW_GRND, PREF_SW_FAC, PTR_SW_WIN,        &
                                 PE_SHADING, OSHAD_DAY,                        &
                                 OSHADE                                        )  
 CHARACTER(LEN=3)  , INTENT(IN)    :: HBEM              ! Building Energy model 'DEF' or 'BEM'
 CHARACTER(LEN=4),   INTENT(IN)    :: HROAD_DIR         ! road direction option
 CHARACTER(LEN=4),   INTENT(IN)    :: HWALL_OPT         ! walls option
REAL, DIMENSION(:), INTENT(IN)    :: PDIR_SW           ! incoming direct solar radiation
REAL, DIMENSION(:), INTENT(IN)    :: PSCA_SW           ! scattered incoming solar rad.
REAL, DIMENSION(:), INTENT(IN)    :: PZENITH           ! solar zenithal angle
REAL, DIMENSION(:), INTENT(IN)    :: PAZIM             ! solar azimuthal angle
REAL, DIMENSION(:), INTENT(IN)    :: PBLD              ! buildings fraction
REAL, DIMENSION(:), INTENT(IN)    :: PGARDEN           ! GARDEN area fraction
REAL, DIMENSION(:), INTENT(IN)    :: PROAD_DIR         ! Road direction
REAL, DIMENSION(:), INTENT(IN)    :: PROAD             ! road fraction
REAL, DIMENSION(:), INTENT(IN)    :: PFRAC_GR          ! green roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_O_HOR       ! wall surf. / hor. surf
REAL, DIMENSION(:), INTENT(IN)    :: PCAN_HW_RATIO     ! canyon    h/W
REAL, DIMENSION(:), INTENT(IN)    :: PALB_ROOF         ! roof albedo
REAL, DIMENSION(:), INTENT(IN)    :: PALB_ROAD         ! road albedo
REAL, DIMENSION(:), INTENT(IN)    :: PSVF_ROAD         ! road sky view factor
REAL, DIMENSION(:), INTENT(IN)    :: PALB_WALL         ! wall albedo
REAL, DIMENSION(:), INTENT(IN)    :: PSVF_WALL         ! wall sky view factor
REAL, DIMENSION(:), INTENT(IN)    :: PFRAC_PANEL       ! Fraction of solar panel on roofs (-)
REAL, DIMENSION(:), INTENT(IN)    :: PALB_PANEL        ! Albedo     of solar panels (-)
REAL, DIMENSION(:), INTENT(IN)    :: PALB_GARDEN       ! GARDEN areas albedo
REAL, DIMENSION(:), INTENT(IN)    :: PSVF_GARDEN       ! GARDEN areas sky view factor
REAL, DIMENSION(:), INTENT(IN)    :: PALB_GREENROOF    ! green roof albedo
REAL, DIMENSION(:), INTENT(IN)    :: PASNOW_ROOF       ! roof snow albedo
REAL, DIMENSION(:), INTENT(IN)    :: PASNOW_ROAD       ! road snow albedo
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROOF          ! snow-covered roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDF_ROOF          ! snow-free    roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROAD          ! snow-covered road fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDF_ROAD          ! snow-free    road fraction
REAL, DIMENSION(:), INTENT(IN)    :: PGR               ! glazing ratio
REAL, DIMENSION(:), INTENT(IN)    :: PABS_WIN          ! solar absortance of windows
REAL, DIMENSION(:), INTENT(IN)     :: PSHGC            ! solar transmittance of windows
REAL, DIMENSION(:), INTENT(IN)     :: PSHGC_SH         ! window + shading solar heat gain coef.
REAL, DIMENSION(:), INTENT(OUT)    :: PALB_WIN         ! solar reflectance of windows
REAL, DIMENSION(:), INTENT(OUT)    :: PABS_SW_ROOF     ! solar radiation absorbed
REAL, DIMENSION(:), INTENT(OUT)    :: PABS_SW_ROAD     ! solar radiation absorbed
REAL, DIMENSION(:), INTENT(OUT)    :: PABS_SW_WALL_A   ! solar radiation absorbed
REAL, DIMENSION(:), INTENT(OUT)    :: PABS_SW_WALL_B   ! solar radiation absorbed
REAL, DIMENSION(:), INTENT(OUT)    :: PABS_SW_GARDEN   ! solar radiation absorbed
REAL, DIMENSION(:), INTENT(OUT)    :: PABS_SW_GREENROOF! solar radiation absorbed by green roofs
REAL, DIMENSION(:), INTENT(OUT)    :: PABS_SW_SNOW_ROOF! solar radiation absorbed
REAL, DIMENSION(:), INTENT(OUT)    :: PABS_SW_SNOW_ROAD! solar radiation absorbed
REAL, DIMENSION(:), INTENT(OUT)    :: PABS_SW_SOLAR_PANEL! solar radiation absorbed
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_ROAD     ! solar radiation received
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_WALL_A   ! solar radiation received
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_WALL_B   ! solar radiation received
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_GARDEN   ! solar radiation received
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_ROOF     ! solar radiation received
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_SNOW_ROAD! solar radiation received
REAL, DIMENSION(:), INTENT(OUT)    :: PDIR_ALB_TOWN    ! town direct albedo
REAL, DIMENSION(:), INTENT(OUT)    :: PSCA_ALB_TOWN    ! town diffuse albedo
REAL, DIMENSION(:), INTENT(OUT)    :: PSW_RAD_GARDEN   ! solar radiation reaching GARDEN areas
REAL, DIMENSION(:), INTENT(OUT)    :: PABS_SW_WIN      ! solar radiation absorbed by windows
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_WIN      ! solar radiation received by windows
REAL, DIMENSION(:), INTENT(OUT)    :: PREF_SW_GRND     ! total solar radiation reflected by ground
REAL, DIMENSION(:), INTENT(OUT)    :: PREF_SW_FAC      ! total solar radiation reflected by wall
REAL, DIMENSION(:), INTENT(OUT)    :: PTR_SW_WIN       ! Solar radiation transmitted
REAL, DIMENSION(:), INTENT(OUT)    :: PE_SHADING       ! Energy that is not reflected 
LOGICAL, DIMENSION(:),INTENT(INOUT):: OSHAD_DAY        ! has shading been necessary this day ?
LOGICAL, DIMENSION(:),INTENT(IN)   :: OSHADE           ! are building conditions favorable for 
REAL, DIMENSION(:), INTENT(IN)     :: PTRAN_WIN        ! window transmittivity
END SUBROUTINE URBAN_SOLAR_ABS
END INTERFACE
END MODULE MODI_URBAN_SOLAR_ABS
