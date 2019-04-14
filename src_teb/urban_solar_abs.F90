!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
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
!   ##########################################################################
!
!!****  *URBAN_SOLAR_ABS*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the solar radiation flux absorbed by roofs, roads and walls.
!     The absorption by roofs is trivial.
!         
!     
!!**  METHOD
!     ------
!
!
!        computation of input solar radiation on each surface
!        ****************************************************
!
!    direct fluxes:
!    -------------
!
!    dir_Rg_road (Wm-2) =   S * 2*theta0/pi
!                         - S *2/tan(zen) * h/W /pi * (1-cos(theta0))
!
!    dir_Rg_wall (Wm-2) =   S / tan(zen) /pi * (1-cos(theta0))
!                         + S * W/h * (1/2 -theta0/pi)
!
!   where zen      is the zenithal angle, from horizon
!         h/W      is the aspect ratio of the canyon
!         S        is the direct solar radiation flux on a horizontal surface
!
!         theta0 = arcsin(min(W/h * tan(zen),1))
!
!   The surfaces will keep (1-a) times these fluxes, and reflect the
!   remaining
!
!    scattered fluxes:
!    ----------------
!
!   sca_Rg_road = sca_Rg * SVF_road
!
!   sca_Rg_wall = sca_Rg * SVF_wall
!
!
!    solar flux and isotropic reflections :
!    ------------------------------------
!
!  after 0 reflection, the absorbed part of the flux is:
!
!      ARg_r(0) = (1-a_r) (sca_Rg_road + dir_Rg_road)
!
!      ARg_w(0) = (1-a_w) (sca_Rg_wall + dir_Rg_wall)
!  
!    and the reflected parts are
!
!      RRg_r(0) = a_r (sca_Rg_road + dir_Rg_road)
!
!      RRg_w(0) = a_w (sca_Rg_wall + dir_Rg_wall)
!
!  after n reflection:
!
!      ARg_r(n) = ARg_r(n-1) + RRg_w(n-1) * (1-  SVF_r)(1-a_r)
!
!      ARg_w(n) = ARg_w(n-1) + RRg_r(n-1) *      SVF_w (1-a_w)
!                            + RRg_w(n-1) * (1-2*SVF_w)(1-a_w)
!
!      RRg_r(n) = (1- SVF_r) a_r RRg_w(n-1)
!
!      RRg_w(n) =     SVF_w  a_w RRg_r(n-1)
!                +(1-2SVF_w) a_w RRg_w(n-1)
!
!
!   i.e.
!                                               n-1
!      ARg_r(n) = ARg_r(0) + (1-  SVF_r)(1-a_r) SUM RRg_w(k)
!                                               k=0
!
!                                               n-1
!      ARg_w(n) = ARg_w(0) +      SVF_w (1-a_w) SUM RRg_r(k)
!                                               k=0
!                                               n-1
!                          + (1-2*SVF_w)(1-a_w) SUM RRg_w(k)
!                                               k=0
!
! with
!
!     n                             n-1
!    SUM RRg_r(k) = (1-  SVF_r) a_r SUM RRg_w(k)      +  RRg_r(0)
!    k=0                            k=0
!
!     n                             n-1
!    SUM RRg_w(k) =      SVF_w  a_w SUM RRg_r(k) 
!    k=0                            k=0
!                                   n-1
!                  +(1-2*SVF_w) a_w SUM RRg_w(k)      +  RRg_w(0)
!                                   k=0
!
!
!   Then
!
!     n                                        n-1
!    SUM RRg_w(k) =  (1-2*SVF_w)       a_w     SUM RRg_w(k)
!    k=0                                       k=0
!                                              n-2
!                  + (1-  SVF_r) SVF_w a_w a_r SUM RRg_w(k) 
!                                              k=0
!
!                  + RRg_w(0) + SVF_w a_w RRg_r(0)
!
!
!
!
!  solving this system, lead after an infinity of reflections/absorptions:
!
!    inf                      RRg_w(0) + SVF_w a_w RRg_r(0)
!    SUM RRg_w(k) = ----------------------------------------------------
!    k=0             1 - (1-2*SVF_w) a_w - (1-  SVF_r) SVF_w a_w a_r
!
!
!    inf            (1-  SVF_r) a_r ( a_w SVF_w RRg_r(0) + RRg_w(0) )
!    SUM RRg_r(k) = ------------------------------------------------------------ + RRg_r(0)
!    k=0             1 - (1-2*SVF_w) a_w - (1-  SVF_r) SVF_w a_w a_r
!
!
! ARg_r(n) and ARg_w(n) follow
!
!
! If snow is present, the albedos in all these formulae (and only these,
! not the final net radiation budget) are modified by the albedo of the
! snow-covered surface.
! 
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    MODD_CST
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
!!      Original    23/01/98 
!!                  21/11/00 (V. Masson)  bug in reflections for roads
!!                     12/02 (A. Lemonsu) bug in diagnostic of albedo
!!                     12/11 (V. Masson ) adds road direction option
!!                     01/12 (V. Masson ) adds 2 different wall direct insulations
!!                     04/12 (G. Pigeon) add PTRAN_WIN 
!!                     09/12 (C. de Munck-A. Lemonsu) add green roofs
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,     ONLY : XPI
USE MODD_BEM_CST,  ONLY : XWIN_SW_MAX
USE MODD_SURF_PAR, ONLY : XUNDEF
!
USE MODI_WINDOW_SHADING
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
 CHARACTER(LEN=3)  , INTENT(IN)    :: HBEM              ! Building Energy model 'DEF' or 'BEM'
 CHARACTER(LEN=4),   INTENT(IN)    :: HROAD_DIR         ! road direction option
              ! 'UNIF' : classical TEB version, all roads uniformely present
              ! 'ORIE' : specified road ORIENTATION (in fact many roads  direction
              !          possible because there can be many TEB patches)
 CHARACTER(LEN=4),   INTENT(IN)    :: HWALL_OPT         ! walls option
              ! 'UNIF' : classical TEB version, all walls are identical
              ! 'TWO ' : the two opposite walls are different & receive different solar energy
REAL, DIMENSION(:), INTENT(IN)    :: PDIR_SW           ! incoming direct solar radiation
REAL, DIMENSION(:), INTENT(IN)    :: PSCA_SW           ! scattered incoming solar rad.
REAL, DIMENSION(:), INTENT(IN)    :: PZENITH           ! solar zenithal angle
REAL, DIMENSION(:), INTENT(IN)    :: PAZIM             ! solar azimuthal angle
!                                                      ! (radian from N, clockwise)
REAL, DIMENSION(:), INTENT(IN)    :: PBLD              ! buildings fraction
REAL, DIMENSION(:), INTENT(IN)    :: PGARDEN           ! GARDEN area fraction
REAL, DIMENSION(:), INTENT(IN)    :: PROAD_DIR         ! Road direction
                                                       ! (deg from N, clockwise)
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
!
REAL, DIMENSION(:), INTENT(IN)    :: PGR               ! glazing ratio
REAL, DIMENSION(:), INTENT(IN)    :: PABS_WIN          ! solar absortance of windows
!
!new arguments for shading
REAL, DIMENSION(:), INTENT(IN)     :: PSHGC            ! solar transmittance of windows
REAL, DIMENSION(:), INTENT(IN)     :: PSHGC_SH         ! window + shading solar heat gain coef.
!
REAL, DIMENSION(:), INTENT(OUT)    :: PALB_WIN         ! solar reflectance of windows
!
REAL, DIMENSION(:), INTENT(OUT)    :: PABS_SW_ROOF     ! solar radiation absorbed
!                                                      ! by snow-free roofs
REAL, DIMENSION(:), INTENT(OUT)    :: PABS_SW_ROAD     ! solar radiation absorbed
!                                                      ! by snow-free roads
REAL, DIMENSION(:), INTENT(OUT)    :: PABS_SW_WALL_A   ! solar radiation absorbed
REAL, DIMENSION(:), INTENT(OUT)    :: PABS_SW_WALL_B   ! solar radiation absorbed
!                                                      ! by walls
REAL, DIMENSION(:), INTENT(OUT)    :: PABS_SW_GARDEN   ! solar radiation absorbed
!                                                      ! by GARDEN areas
REAL, DIMENSION(:), INTENT(OUT)    :: PABS_SW_GREENROOF! solar radiation absorbed by green roofs
REAL, DIMENSION(:), INTENT(OUT)    :: PABS_SW_SNOW_ROOF! solar radiation absorbed
!                                                      ! by snow-covered roofs
REAL, DIMENSION(:), INTENT(OUT)    :: PABS_SW_SNOW_ROAD! solar radiation absorbed
!                                                      ! by snow-covered roads
REAL, DIMENSION(:), INTENT(OUT)    :: PABS_SW_SOLAR_PANEL! solar radiation absorbed
!                                                      ! by solar panels
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_ROAD     ! solar radiation received
!                                                      ! by snow-free roads
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_WALL_A   ! solar radiation received
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_WALL_B   ! solar radiation received
!                                                      ! by snow-free walls
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_GARDEN   ! solar radiation received
!                                                      ! by GARDEN areas
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_ROOF     ! solar radiation received
!                                                      ! by ROOF areas (below solar panels if any)
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_SNOW_ROAD! solar radiation received
!                                                      ! by snow-covered roads
REAL, DIMENSION(:), INTENT(OUT)    :: PDIR_ALB_TOWN    ! town direct albedo
REAL, DIMENSION(:), INTENT(OUT)    :: PSCA_ALB_TOWN    ! town diffuse albedo
!
REAL, DIMENSION(:), INTENT(OUT)    :: PSW_RAD_GARDEN   ! solar radiation reaching GARDEN areas
REAL, DIMENSION(:), INTENT(OUT)    :: PABS_SW_WIN      ! solar radiation absorbed by windows
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_WIN      ! solar radiation received by windows

REAL, DIMENSION(:), INTENT(OUT)    :: PREF_SW_GRND     ! total solar radiation reflected by ground
REAL, DIMENSION(:), INTENT(OUT)    :: PREF_SW_FAC      ! total solar radiation reflected by wall
REAL, DIMENSION(:), INTENT(OUT)    :: PTR_SW_WIN       ! Solar radiation transmitted
                                                       ! through glazing [W/m2(bld)]
!new arguments for shading
REAL, DIMENSION(:), INTENT(OUT)    :: PE_SHADING       ! Energy that is not reflected 
                                                       ! by the shading, nor transmitted through
                                                       ! the bld, nor absorbed by the
                                                       ! [W/m2(win)]
LOGICAL, DIMENSION(:),INTENT(INOUT):: OSHAD_DAY        ! has shading been necessary this day ?
LOGICAL, DIMENSION(:),INTENT(IN)   :: OSHADE           ! are building conditions favorable for 
!                                                      ! shading (independantly of solar irradiance) ?
REAL, DIMENSION(:), INTENT(IN)     :: PTRAN_WIN        ! window transmittivity
!
!*      0.2    declarations of local variables
!
!                                                           
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZDIR_SW               ! direct and diffuse incoming radiation
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZSCA_SW               ! with a minimum to compute albedo
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZTANZEN               ! tangente of solar zenithal angle
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZTHETA0               ! canyon angle for
!                                                       ! which solar
!                                                       ! radiation
!                                                       ! reaches the road
!
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZAALB_ROAD            ! averaged albedo
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZDIR_SW_ROAD          ! direct radiation reaching
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZDIR_SW_WALL_A        ! road, wall A,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZDIR_SW_WALL_B        ! wall B,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZDIR_SW_GARDEN        ! GARDEN areas,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZDIR_SW_WALL          ! and on average on 2 walls
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZSCA_SW_ROAD          ! diffuse radiation reaching
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZSCA_SW_WALL          ! road, wall,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZSCA_SW_GARDEN        ! and GARDEN areas
!
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_ROOF      ! solar radiation
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_ROAD      ! absorbed by roofs,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_WALL_A    ! road, wall A,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_WALL_B    ! wall B,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_WALL      ! both walls on average,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_GARDEN    ! GARDEN areas,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_GREENROOF ! green roof areas,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_PANEL     ! solar panels,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_WIN       ! window (abs+trans), and snow
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_SNOW_ROOF ! over roof, wall,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_SNOW_ROAD ! and GARDEN areas
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_ROOF      ! solar radiation
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_ROAD      ! absorbed by roofs,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_WALL      ! road, wall,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_GARDEN    ! GARDEN areas,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_GREENROOF ! green roof areas,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_PANEL     ! solar panels,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_WIN       ! window (abs+trans), and snow
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_SNOW_ROOF ! over roof and wall,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_SNOW_ROAD ! coming from diffuse rad.
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZDW                   ! difference of radiation
!                                                       ! absorbed by the 2 walls
!
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZROAD                 !
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZROAD_DIR             ! Road direction
!                                                       ! (radian from N, clockwise)
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZGARDEN               !
!
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZREC_DIR_SW_WIN
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZREC_SCA_SW_WIN
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZAALB_WALL
!
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZTRAN_WIN             ! solar transmittivity of windows
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_WIN              ! solar transmittivity of windows
LOGICAL, DIMENSION(SIZE(PDIR_SW)) :: G_EFF_SHAD         !TRUE if shading should be active considering actual rad
!
INTEGER                        :: JJ                    ! loop index
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('URBAN_SOLAR_ABS',0,ZHOOK_HANDLE)
PABS_SW_ROOF      (:) = 0.
PABS_SW_ROAD      (:) = 0.
PABS_SW_WALL_A    (:) = 0.
PABS_SW_WALL_B    (:) = 0.
PABS_SW_GARDEN    (:) = 0.
PABS_SW_GREENROOF (:) = 0.
PABS_SW_SNOW_ROOF (:) = 0.
PABS_SW_SNOW_ROAD (:) = 0.
PABS_SW_WIN       (:) = 0.
PREC_SW_WIN       (:) = 0.
!
PREC_SW_ROAD      (:) = 0.
PREC_SW_WALL_A    (:) = 0.
PREC_SW_WALL_B    (:) = 0.
PREC_SW_GARDEN    (:) = 0.
PREC_SW_SNOW_ROAD (:) = 0.
PREC_SW_ROOF      (:) = 0.
!
ZDIR_SW = MAX(PDIR_SW,0.)
ZSCA_SW = MAX(PSCA_SW,0.000001)
!
ZROAD_DIR = PROAD_DIR(:) * XPI/180.
!
DO JJ=1,SIZE(PROAD)
!
  IF (PROAD(JJ)+PGARDEN(JJ).NE.0.) THEN
    ZROAD(JJ)  = PROAD(JJ)  / (PROAD(JJ)+PGARDEN(JJ))
    ZGARDEN(JJ) =  PGARDEN(JJ) / (PROAD(JJ)+PGARDEN(JJ))
  ELSE
    ZROAD(JJ)=0.
    ZGARDEN(JJ)=0.
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*      1.     SOLAR RADIATIONS FOR ROOFS
!              --------------------------
!
!* One supposes that solar panels, if present, intercept all solar radiation
!
  ZABS_DIR_SW_PANEL    (JJ) = ZDIR_SW(JJ) * (1. - PALB_PANEL    (JJ))
  ZABS_SCA_SW_PANEL    (JJ) = ZSCA_SW(JJ) * (1. - PALB_PANEL    (JJ))
!
!* solar energy received by the surfaces below solar panels
  ZABS_DIR_SW_ROOF     (JJ) = ZDIR_SW(JJ) * (1. - PALB_ROOF     (JJ)) * (1.-PFRAC_PANEL(JJ))
  ZABS_DIR_SW_SNOW_ROOF(JJ) = ZDIR_SW(JJ) * (1. - PASNOW_ROOF   (JJ)) * (1.-PFRAC_PANEL(JJ))
  ZABS_DIR_SW_GREENROOF(JJ) = ZDIR_SW(JJ) * (1. - PALB_GREENROOF(JJ)) * (1.-PFRAC_PANEL(JJ))
  ZABS_SCA_SW_ROOF     (JJ) = ZSCA_SW(JJ) * (1. - PALB_ROOF     (JJ)) * (1.-PFRAC_PANEL(JJ))
  ZABS_SCA_SW_SNOW_ROOF(JJ) = ZSCA_SW(JJ) * (1. - PASNOW_ROOF   (JJ)) * (1.-PFRAC_PANEL(JJ))
  ZABS_SCA_SW_GREENROOF(JJ) = ZSCA_SW(JJ) * (1. - PALB_GREENROOF(JJ)) * (1.-PFRAC_PANEL(JJ))
!
!-------------------------------------------------------------------------------
!
!*      2.     SOLAR RADIATIONS FOR ROADS AND WALLS
!              ------------------------------------
!
  IF (ABS(0.5*XPI-PZENITH(JJ)) <  1.E-6) THEN
    IF(0.5*XPI-PZENITH(JJ) > 0.)  ZTANZEN(JJ)=TAN(0.5*XPI-1.E-6)
    IF(0.5*XPI-PZENITH(JJ) <= 0.) ZTANZEN(JJ)=TAN(0.5*XPI+1.E-6)
  ELSEIF (ABS(PZENITH(JJ)) <  1.E-6) THEN
    ZTANZEN(JJ)=SIGN(1.,PZENITH(JJ))*TAN(1.E-6)
  ELSE
    ZTANZEN(JJ) = TAN(PZENITH(JJ))
  ENDIF
!
!
  IF (PBLD(JJ) .GT. 0.) THEN
!
!*      2.1    radiation coefficients
!              ----------------------
!
IF (HROAD_DIR=='UNIF') THEN
  ZTHETA0(JJ) = ASIN( MIN(ABS( 1./ZTANZEN(JJ))/PCAN_HW_RATIO(JJ), 1. ) )
!
!*      2.2    direct solar radiation received by roads and GARDEN areas
!               -------------------------------------------------------
!
  ZDIR_SW_ROAD (JJ) = (  ZDIR_SW(JJ) * 2. * ZTHETA0(JJ) / XPI                 &
                       - ZDIR_SW(JJ) * 2. * ZTANZEN(JJ) / XPI                 &
                                 * PCAN_HW_RATIO(JJ) * (1.-COS(ZTHETA0(JJ)))  )  
!
  ZDIR_SW_GARDEN(JJ) = ( ZDIR_SW(JJ) * 2. * ZTHETA0(JJ) / XPI                 &
                       - ZDIR_SW(JJ) * 2. * ZTANZEN(JJ) / XPI                 &
                                 * PCAN_HW_RATIO(JJ) * (1.-COS(ZTHETA0(JJ)))  )  
ELSE
  ZDIR_SW_ROAD (JJ) = ZDIR_SW(JJ) * MAX(0.,                           &
   1.-PCAN_HW_RATIO(JJ)*ZTANZEN(JJ)*ABS(SIN(PAZIM(JJ)-ZROAD_DIR(JJ))) &
                                       )
  ZDIR_SW_GARDEN(JJ) = ZDIR_SW_ROAD(JJ)

END IF
!
!*      2.3    direct solar radiation received by walls
!              ----------------------------------------
!
  ZDIR_SW_WALL(JJ) = (ZDIR_SW(JJ) - (ZDIR_SW_ROAD(JJ)*ZROAD(JJ)+ZDIR_SW_GARDEN(JJ)*ZGARDEN(JJ))) &
                      * 0.5 / PCAN_HW_RATIO(JJ)  
!

  ELSE
!
   ZDIR_SW_ROAD (JJ) = ZDIR_SW(JJ)
   ZDIR_SW_GARDEN(JJ) = ZDIR_SW(JJ)
   ZDIR_SW_WALL (JJ) = 0.
!
  ENDIF
!
IF (HROAD_DIR=='UNIF' .OR. HWALL_OPT=='UNIF') THEN
!* if walls are averaged, then
  ZDIR_SW_WALL_A(JJ) = ZDIR_SW_WALL(JJ)
  ZDIR_SW_WALL_B(JJ) = ZDIR_SW_WALL(JJ)
ELSE
!* if walls are separated, then radiation reaches the wall facing sun
! Note that wall A is the one facing mostly to the South (depending to
! road orientation), and wall B in the one facing mostly to the North
!
! In case of N-S road, wall A is the West  wall (= East-facing  wall), 
!                  and wall B is the East  wall (= West-facing  wall)
! In case of E-W road, wall A is the North wall (= South-facing wall), 
!                  and wall B is the South wall (= North-facing wall)
  IF (SIN(PAZIM(JJ)-ZROAD_DIR(JJ))>0.) THEN
    ZDIR_SW_WALL_A(JJ) = 2.* ZDIR_SW_WALL(JJ)
    ZDIR_SW_WALL_B(JJ) = 0.
  ELSE
    ZDIR_SW_WALL_A(JJ) = 0.
    ZDIR_SW_WALL_B(JJ) = 2.* ZDIR_SW_WALL(JJ)
  END IF
END IF
!
!
!
!*      2.4    diffuse solar radiation received by roads and GARDEN areas
!              ---------------------------------------------------------
!
  ZSCA_SW_ROAD (JJ) = ZSCA_SW(JJ) * PSVF_ROAD(JJ)
!
  ZSCA_SW_GARDEN(JJ) = ZSCA_SW(JJ) * PSVF_GARDEN(JJ)
!
!*      2.5    diffuse solar radiation received by walls
!              -----------------------------------------
!
  ZSCA_SW_WALL(JJ) = ZSCA_SW(JJ) * PSVF_WALL(JJ)
!
!*      2.6    total solar radiation received by GARDEN areas
!              ---------------------------------------------
!
  PSW_RAD_GARDEN(JJ) = ZDIR_SW_GARDEN(JJ) + ZSCA_SW_GARDEN(JJ)
!
!*      2.7    averaged albedos when snow is present
!              -------------------------------------
!
  ZAALB_ROAD (JJ) =  PDF_ROAD (JJ) * PALB_ROAD   (JJ) &
                  + PDN_ROAD (JJ) * PASNOW_ROAD (JJ)  
!
!
ENDDO
!
!*      2.7b    averaged facade albedo
!              -------------------------------------
!
IF (HBEM=='BEM') THEN
    !
    ZTRAN_WIN         (:) = PTRAN_WIN(:)
    !
    G_EFF_SHAD(:) = OSHADE(:)                              .AND. &
          (ZDIR_SW_WALL(:) + ZSCA_SW_WALL(:) > XWIN_SW_MAX) 
    !
    OSHAD_DAY(:)  = G_EFF_SHAD(:) .OR. OSHAD_DAY(:)
    !
    CALL WINDOW_SHADING(PSHGC, PSHGC_SH, OSHAD_DAY, PALB_WALL, &
                  PABS_WIN, ZABS_WIN, PALB_WIN, ZTRAN_WIN      )
    !
ELSE
  !
  ZABS_WIN (:) = 0.
  PALB_WIN (:) = 0.
  ZTRAN_WIN(:) = 0.
  !
ENDIF
!
ZAALB_WALL(:) =  PGR(:) * PALB_WIN(:) + (1.-PGR(:)) * PALB_WALL(:)
!
!*      2.8    absorption of direct incoming solar radiation
!              ---------------------------------------------
!
!
 CALL SOLAR_REFLECTIONS(ZDIR_SW_ROAD,ZDIR_SW_WALL,                                 &
                         ZDIR_SW_GARDEN,                                          &
                         ZABS_DIR_SW_ROAD, ZABS_DIR_SW_SNOW_ROAD,                 &
                         ZABS_DIR_SW_WALL, ZABS_DIR_SW_GARDEN, ZABS_DIR_SW_WIN    )  
!
IF (HROAD_DIR=='UNIF' .OR. HWALL_OPT=='UNIF') THEN
!* if walls are averaged, then
  ZABS_DIR_SW_WALL_A = ZABS_DIR_SW_WALL
  ZABS_DIR_SW_WALL_B = ZABS_DIR_SW_WALL
ELSE
!* if walls are separated, then radiation reaches the wall facing sun
! Note that wall A is the one facing mostly to the North (depending to
! road orientation), and wall B in the one facing mostly to the South.
  ZDW = (1.-PALB_WALL(:)) * ZAALB_WALL(:) * (1.-2.*PSVF_WALL(:)) &
       / (1.+ZAALB_WALL(:)*(1.-2.*PSVF_WALL(:)))                 &
       * 0.5 * (ZDIR_SW_WALL_A(:)-ZDIR_SW_WALL_B(:))            &
       + 0.5 * (1.-PALB_WALL(:)) * (ZDIR_SW_WALL_A-ZDIR_SW_WALL_B)
  ZABS_DIR_SW_WALL_A = ZABS_DIR_SW_WALL + ZDW
  ZABS_DIR_SW_WALL_B = ZABS_DIR_SW_WALL - ZDW
END IF
!
!*      2.9    absorption of diffuse incoming solar radiation
!              ----------------------------------------------
!
 CALL SOLAR_REFLECTIONS(ZSCA_SW_ROAD,ZSCA_SW_WALL, ZSCA_SW_GARDEN,                 &
                         ZABS_SCA_SW_ROAD, ZABS_SCA_SW_SNOW_ROAD,                 &
                         ZABS_SCA_SW_WALL, ZABS_SCA_SW_GARDEN, ZABS_SCA_SW_WIN    )  
!
! solar flux reflected for wall and road
!
PREF_SW_GRND = ZROAD   * PALB_ROAD   / (1. - PALB_ROAD  )                   * &
                                   (ZABS_DIR_SW_ROAD + ZABS_SCA_SW_ROAD)      &
             + ZGARDEN * PALB_GARDEN / (1. - PALB_GARDEN)                   * &
                                   (ZABS_DIR_SW_GARDEN + ZABS_SCA_SW_GARDEN)
!
PREF_SW_FAC = (1 - PGR) * PALB_WALL  / (1. - PALB_WALL)                     * &
                                     (ZABS_DIR_SW_WALL + ZABS_SCA_SW_WALL)    &
            +      PGR  * PALB_WIN   / (1 - PALB_WIN  )                     * &
                                   (ZABS_DIR_SW_WIN + ZABS_SCA_SW_WIN)
!
!-------------------------------------------------------------------------------
!
!*      3.     Town albedo
!              -----------
!
!*      3.1    direct albedo
!              -------------
!
 CALL TOWN_ALBEDO(ZDIR_SW,ZABS_DIR_SW_ROOF,ZABS_DIR_SW_SNOW_ROOF,              &
                   ZABS_DIR_SW_ROAD, ZABS_DIR_SW_SNOW_ROAD,ZABS_DIR_SW_WALL,  &
                   ZABS_DIR_SW_GARDEN, ZABS_DIR_SW_GREENROOF, ZABS_DIR_SW_WIN,&
                   ZABS_DIR_SW_PANEL, PDIR_ALB_TOWN                           )  
!
!*      3.2    direct albedo
!              -------------
!
 CALL TOWN_ALBEDO(ZSCA_SW,ZABS_SCA_SW_ROOF,ZABS_SCA_SW_SNOW_ROOF,              &
                   ZABS_SCA_SW_ROAD, ZABS_SCA_SW_SNOW_ROAD,ZABS_SCA_SW_WALL,  &
                   ZABS_SCA_SW_GARDEN, ZABS_SCA_SW_GREENROOF, ZABS_SCA_SW_WIN,&
                   ZABS_SCA_SW_PANEL, PSCA_ALB_TOWN                           )  
!
!-------------------------------------------------------------------------------
!
!*      4.     Trivial cases
!              -------------
!
WHERE(PDIR_SW(:)==0.)
  ZABS_DIR_SW_ROOF      (:) = 0.
  ZABS_DIR_SW_ROAD      (:) = 0.
  ZABS_DIR_SW_WALL_A    (:) = 0.
  ZABS_DIR_SW_WALL_B    (:) = 0.
  ZABS_DIR_SW_GARDEN    (:) = 0.
  ZABS_DIR_SW_GREENROOF (:) = 0.
  ZABS_DIR_SW_PANEL     (:) = 0.
  ZABS_DIR_SW_WIN       (:) = 0.
  ZABS_DIR_SW_SNOW_ROOF (:) = 0.
  ZABS_DIR_SW_SNOW_ROAD (:) = 0.
END WHERE
!
WHERE(PSCA_SW(:)==0.)
  ZABS_SCA_SW_ROOF      (:) = 0.
  ZABS_SCA_SW_ROAD      (:) = 0.
  ZABS_SCA_SW_WALL      (:) = 0.
  ZABS_SCA_SW_GARDEN    (:) = 0.
  ZABS_SCA_SW_GREENROOF (:) = 0.
  ZABS_SCA_SW_PANEL     (:) = 0.
  ZABS_SCA_SW_WIN       (:) = 0.
  ZABS_SCA_SW_SNOW_ROOF (:) = 0.
  ZABS_SCA_SW_SNOW_ROAD (:) = 0.
END WHERE
!
WHERE (PDIR_ALB_TOWN==XUNDEF) PDIR_ALB_TOWN = PSCA_ALB_TOWN
!
!-------------------------------------------------------------------------------
!
DO JJ=1,SIZE(PROAD)
!
!*      5.     Total solar radiation absorbed by each surface
!              ----------------------------------------------
!
! solar radiation absorbed by roofs
!
  PABS_SW_ROOF     (JJ) = ZABS_DIR_SW_ROOF     (JJ) + ZABS_SCA_SW_ROOF     (JJ)
!
! solar radiation absorbed by roads
!
  PABS_SW_ROAD     (JJ) = ZABS_DIR_SW_ROAD     (JJ) + ZABS_SCA_SW_ROAD     (JJ)
!
! solar radiation absorbed by GARDEN areas
!
  PABS_SW_GARDEN   (JJ) = ZABS_DIR_SW_GARDEN   (JJ) + ZABS_SCA_SW_GARDEN   (JJ)
!
! solar radiation absorbed by GREENROOF areas
!
  PABS_SW_GREENROOF(JJ) = ZABS_DIR_SW_GREENROOF(JJ) + ZABS_SCA_SW_GREENROOF(JJ)
!
! solar radiation absorbed by solar panels
!
  PABS_SW_SOLAR_PANEL(JJ)= ZABS_DIR_SW_PANEL    (JJ) + ZABS_SCA_SW_PANEL    (JJ)
!
! solar radiation absorbed by walls
!
  PABS_SW_WALL_A   (JJ) = ZABS_DIR_SW_WALL_A   (JJ) + ZABS_SCA_SW_WALL     (JJ)
  PABS_SW_WALL_B   (JJ) = ZABS_DIR_SW_WALL_B   (JJ) + ZABS_SCA_SW_WALL     (JJ)
!
! solar radiation absorbed (but not transmitted) by windows
!
  ZREC_DIR_SW_WIN(JJ) = ZABS_DIR_SW_WIN(JJ) / (1.-PALB_WIN(JJ))
  ZREC_SCA_SW_WIN(JJ) = ZABS_SCA_SW_WIN(JJ) / (1.-PALB_WIN(JJ))
!
  PABS_SW_WIN      (JJ) = (ZREC_DIR_SW_WIN  (JJ) + ZREC_SCA_SW_WIN   (JJ)) * ZABS_WIN(JJ)
!
!
! solar radiation absorbed by snow on roofs
!
  PABS_SW_SNOW_ROOF (JJ) = ZABS_DIR_SW_SNOW_ROOF (JJ) + ZABS_SCA_SW_SNOW_ROOF (JJ)
!
! solar radiation absorbed by snow on roads
!
  PABS_SW_SNOW_ROAD (JJ) = ZABS_DIR_SW_SNOW_ROAD (JJ) + ZABS_SCA_SW_SNOW_ROAD (JJ)
!
!-------------------------------------------------------------------------------
!
!*      6.     total solar radiation received by roads and GARDEN areas
!              -------------------------------------------------------
!
  PREC_SW_ROAD      (JJ) = PABS_SW_ROAD      (JJ)/(1.-PALB_ROAD   (JJ))
!
  PREC_SW_SNOW_ROAD (JJ) = PABS_SW_SNOW_ROAD (JJ)/(1.-PASNOW_ROAD (JJ))
!
  PREC_SW_WALL_A    (JJ) = PABS_SW_WALL_A    (JJ)/(1.-PALB_WALL   (JJ))
  PREC_SW_WALL_B    (JJ) = PABS_SW_WALL_B    (JJ)/(1.-PALB_WALL   (JJ))
!
  PREC_SW_WIN       (JJ) = ZREC_DIR_SW_WIN(JJ) + ZREC_SCA_SW_WIN(JJ)
!
  PREC_SW_GARDEN    (JJ) = PABS_SW_GARDEN    (JJ)/(1.-PALB_GARDEN (JJ))
!
!*      6.2    total solar radiation received by roof surfaces below solar panels
!
  PREC_SW_ROOF      (JJ) = (PDIR_SW(JJ) + PSCA_SW(JJ)) * (1.-PFRAC_PANEL(JJ))
!
!-------------------------------------------------------------------------------
!
!*      7.     total solar radiation transmitted inside building
!*             and energy not ref., nor absorbed, nor transmitted
!              --------------------------------------------------
!
    PTR_SW_WIN(JJ) = PREC_SW_WIN(JJ) * ZTRAN_WIN(JJ)
!               [W/m2(bld)]
    PE_SHADING(JJ) = PREC_SW_WIN(JJ) * (1. - PALB_WIN(JJ) - ZABS_WIN(JJ) - ZTRAN_WIN(JJ))
!               [W/m2(win)]

ENDDO
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('URBAN_SOLAR_ABS',1,ZHOOK_HANDLE)
 CONTAINS
!
!-------------------------------------------------------------------------------
SUBROUTINE SOLAR_REFLECTIONS(ZSW_ROAD,ZSW_WALL,ZSW_GARDEN,             &
                               ZABS_SW_ROAD,ZABS_SW_SNOW_ROAD,         &
                               ZABS_SW_WALL,                           &
                               ZABS_SW_GARDEN, ZABS_SW_WIN             )  
!
REAL, DIMENSION(:), INTENT(IN) :: ZSW_ROAD          ! solar radiation received by road,
REAL, DIMENSION(:), INTENT(IN) :: ZSW_WALL          ! wall, and GARDEN areas 
REAL, DIMENSION(:), INTENT(IN) :: ZSW_GARDEN        ! before reflection
REAL, DIMENSION(:), INTENT(OUT):: ZABS_SW_ROAD      ! solar radiation absorbed by
REAL, DIMENSION(:), INTENT(OUT):: ZABS_SW_SNOW_ROAD ! solar radiation absorbed by
REAL, DIMENSION(:), INTENT(OUT):: ZABS_SW_WALL      ! road, snow over road, and wall 
REAL, DIMENSION(:), INTENT(OUT):: ZABS_SW_GARDEN    ! solar radiation absorbed by garden
REAL, DIMENSION(:), INTENT(OUT):: ZABS_SW_WIN       ! solar radiation absorbed by window
!
REAL, DIMENSION(SIZE(ZSW_ROAD)) :: ZREF0_SW_ROAD    ! first solar reflection
REAL, DIMENSION(SIZE(ZSW_ROAD)) :: ZREF0_SW_WALL    ! against road, wall
REAL, DIMENSION(SIZE(ZSW_ROAD)) :: ZREF0_SW_GARDEN  ! and GARDEN areas
REAL, DIMENSION(SIZE(ZSW_ROAD)) :: ZSREF_SW_ROAD    ! sum of all reflections
REAL, DIMENSION(SIZE(ZSW_ROAD)) :: ZSREF_SW_WALL    ! against road, wall,
REAL, DIMENSION(SIZE(ZSW_ROAD)) :: ZSREF_SW_GARDEN  ! and GARDEN areas
!
REAL, DIMENSION(SIZE(ZSW_ROAD)) :: ZWORK1, ZWORK2, ZDENOM
INTEGER :: JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SOLAR_REFLECTIONS',0,ZHOOK_HANDLE)
!
DO JJ=1,SIZE(ZSW_ROAD)
!*      A.     first solar radiation reflection
!              --------------------------------
!
  ZREF0_SW_ROAD(JJ)  = ZAALB_ROAD(JJ)   * ZSW_ROAD(JJ) 
!
  ZREF0_SW_GARDEN(JJ) = PALB_GARDEN(JJ) * ZSW_GARDEN(JJ)
!
  ZREF0_SW_WALL(JJ)  = ZAALB_WALL(JJ)   * ZSW_WALL(JJ)
!
!*      B.     sum of solar radiation reflected
!              --------------------------------
!

  ZDENOM(JJ) = 1. - (1.-2.*PSVF_WALL(JJ)) * PALB_WALL(JJ) &
                 - (1. - PSVF_ROAD(JJ))* &
      PSVF_WALL(JJ)*PALB_WALL(JJ)*ZAALB_ROAD(JJ)*ZROAD(JJ) &
                 - (1. - PSVF_GARDEN(JJ))* &
      PSVF_WALL(JJ)*PALB_WALL(JJ)*PALB_GARDEN(JJ)*ZGARDEN(JJ)  

 ZWORK1(JJ) = PSVF_WALL(JJ) * PALB_WALL(JJ) * ZROAD(JJ)
 ZWORK2(JJ) = PSVF_WALL(JJ) * PALB_WALL(JJ) * ZGARDEN(JJ)
!
!
  ZSREF_SW_WALL(JJ)  = (  ZREF0_SW_WALL(JJ)   &
                           + ZWORK1(JJ) *ZREF0_SW_ROAD(JJ)   &
                           + ZWORK2(JJ) *ZREF0_SW_GARDEN(JJ)) &
                        / ZDENOM(JJ)  

  ZSREF_SW_ROAD(JJ)  = ( (1.-   PSVF_ROAD(JJ) ) * ZAALB_ROAD(JJ)                     &
                                                * ZREF0_SW_WALL(JJ)                  &
                         +(1.-   PSVF_ROAD(JJ) ) * ZAALB_ROAD(JJ)                     &
                                                * ZWORK1(JJ)                     &
                                                * ZREF0_SW_ROAD(JJ)                  &
                         +(1.-   PSVF_ROAD(JJ) ) * ZAALB_ROAD(JJ)                     &
                                                * ZWORK2(JJ)                        &
                                                * ZREF0_SW_GARDEN(JJ)               )&
                        / ZDENOM(JJ) &
                        + ZREF0_SW_ROAD(JJ)   

  ZSREF_SW_GARDEN(JJ) = ( (1.- PSVF_GARDEN(JJ)) * PALB_GARDEN(JJ)                    &
                                                * ZREF0_SW_WALL(JJ)                  &
                         +(1.-  PSVF_GARDEN(JJ)) * PALB_GARDEN(JJ)                    &
                                                * ZWORK1(JJ)                        &
                                                * ZREF0_SW_ROAD(JJ)                  &
                         +(1.-  PSVF_GARDEN(JJ)) * PALB_GARDEN(JJ)                    &
                                                * ZWORK2(JJ)                    &
                                                * ZREF0_SW_GARDEN(JJ)               )&
                        / ZDENOM(JJ)  &
                        + ZREF0_SW_GARDEN(JJ)  
!
!
!*      C.     total solar radiation received by roads and GARDEN areas
!              -------------------------------------------------------
!
  ZABS_SW_ROAD(JJ)       = (1.-PALB_ROAD(JJ)   )                              &
                          * (   ZSW_ROAD(JJ)                                   &
                            + ZSREF_SW_WALL(JJ)  * (1.- PSVF_ROAD(JJ) ) )  
!
  ZABS_SW_SNOW_ROAD(JJ)  = (1.-PASNOW_ROAD(JJ) )                              &
                          * (   ZSW_ROAD(JJ)                                   &
                            + ZSREF_SW_WALL(JJ)  * (1.- PSVF_ROAD(JJ) ) )  
!
  ZABS_SW_GARDEN(JJ)      = (1.-PALB_GARDEN(JJ))                               &
                          * (   ZSW_GARDEN(JJ)                                  &
                            + ZSREF_SW_WALL(JJ)  * (1.- PSVF_GARDEN(JJ)) )  
!
!
!*      D.     total solar radiation received by walls
!              ---------------------------------------
!
  ZABS_SW_WALL(JJ)      = (1.-PALB_WALL(JJ))                                    &
                         * (   ZSW_WALL(JJ)                                      &
                           + ZSREF_SW_ROAD(JJ)   *        PSVF_WALL(JJ)*ZROAD(JJ)  &
                           + ZSREF_SW_GARDEN(JJ) *        PSVF_WALL(JJ)*ZGARDEN(JJ) &
                           + ZSREF_SW_WALL(JJ)   * (1.-2.*PSVF_WALL(JJ))          )  
!
  ZABS_SW_WIN (JJ)      = (1.-PALB_WIN (JJ))                                    &
                         * (   ZSW_WALL(JJ)                                      &
                           + ZSREF_SW_ROAD(JJ)   *        PSVF_WALL(JJ)*ZROAD(JJ)  &
                           + ZSREF_SW_GARDEN(JJ) *        PSVF_WALL(JJ)*ZGARDEN(JJ) &
                           + ZSREF_SW_WALL(JJ)   * (1.-2.*PSVF_WALL(JJ))          )  
!
ENDDO
!
IF (LHOOK) CALL DR_HOOK('SOLAR_REFLECTIONS',1,ZHOOK_HANDLE)
!
END SUBROUTINE SOLAR_REFLECTIONS
!
!-------------------------------------------------------------------------------
!
SUBROUTINE TOWN_ALBEDO(ZSW,ZABS_SW_ROOF,ZABS_SW_SNOW_ROOF,              &
                         ZABS_SW_ROAD, ZABS_SW_SNOW_ROAD,ZABS_SW_WALL,  &
                         ZABS_SW_GARDEN, ZABS_SW_GREENROOF, ZABS_SW_WIN,&
                         ZABS_SW_PANEL, ZALBEDO                         )  
!
REAL, DIMENSION(:), INTENT(IN) :: ZSW               ! incoming solar radiation
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_ROOF      ! solar radiation absorbed by roofs
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_ROAD      ! solar radiation absorbed by roads
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_WALL      ! solar radiation absorbed by walls
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_WIN       ! solar radiation absorbed & transmitted by windows
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_GARDEN    ! solar radiation absorbed by GARDEN areas
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_GREENROOF ! solar radiation absorbed by green roof areas
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_SNOW_ROOF ! solar radiation absorbed by roof snow
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_SNOW_ROAD ! solar radiation absorbed by road snow
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_PANEL     ! solar radiation absorbed by solar panels
REAL, DIMENSION(:), INTENT(OUT):: ZALBEDO           ! town averaged albedo

REAL, DIMENSION(SIZE(ZSW))     :: ZSW_UP            ! outgoing solar radiation
INTEGER :: JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('TOWN_ALBEDO',0,ZHOOK_HANDLE)
DO JJ=1,SIZE(ZSW)

  ZSW_UP(JJ) = ZSW(JJ)                                                           &
            - ( PBLD(JJ)   *(1.-PFRAC_GR(JJ))*PDF_ROOF(JJ) *ZABS_SW_ROOF(JJ)     &
               +PBLD(JJ)   *(1.-PFRAC_GR(JJ))*PDN_ROOF(JJ) *ZABS_SW_SNOW_ROOF(JJ)&
               +PBLD(JJ)   *    PFRAC_GR(JJ)               *ZABS_SW_GREENROOF(JJ)&
               +PBLD(JJ)   *    PFRAC_PANEL(JJ)            *ZABS_SW_PANEL(JJ)    &
               +PROAD(JJ)                    *PDF_ROAD(JJ) *ZABS_SW_ROAD (JJ)    &
               +PROAD(JJ)                    *PDN_ROAD(JJ) *ZABS_SW_SNOW_ROAD(JJ)&
               +PGARDEN(JJ)                                *ZABS_SW_GARDEN(JJ)   &
               +PWALL_O_HOR(JJ)              *(1.-PGR(JJ)) *ZABS_SW_WALL(JJ)     &
               +PWALL_O_HOR(JJ)              *    PGR(JJ)  *ZABS_SW_WIN (JJ)     )  
!
  IF (ZSW(JJ)>0.) THEN
    ZALBEDO(JJ)  = ZSW_UP(JJ) / ZSW(JJ)
  ELSE
    ZALBEDO(JJ)  = XUNDEF
  END IF
!
ENDDO
IF (LHOOK) CALL DR_HOOK('TOWN_ALBEDO',1,ZHOOK_HANDLE)
!
END SUBROUTINE TOWN_ALBEDO
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE URBAN_SOLAR_ABS
