!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE URBAN_SOLAR_ABS(TOP, T, B, DMT, PDIR_SW, PSCA_SW, PZENITH, PAZIM,   &
                               PFRAC_PANEL, PALB_PANEL, PALB_GD, PSVF_GD, PALB_GRF, &
                               PDN_RF, PDF_RF, PDN_RD, PDF_RD, PREC_SW_RD,          &
                               PREC_SW_SN_RD, PREC_SW_WL_A, PREC_SW_WL_B,           &
                               PREC_SW_GD, PREC_SW_RF, PDIR_ALB_TWN, PSCA_ALB_TWN,  &
                               PSW_RAD_GD, PREC_SW_WIN, PREF_SW_GRND, PREF_SW_FAC,  &
                               PE_SHADING, OSHAD_DAY, OSHADE, OALB_ONLY )  
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
!!                     04/12 (G. Pigeon) add B%XTRAN_WIN 
!!                     09/12 (C. de Munck-A. Lemonsu) add green roofs
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_t
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
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
              ! 'UNIF' : classical TEB version, all walls are identical
              ! 'TWO ' : the two opposite walls are different & receive different solar energy
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(BEM_t), INTENT(INOUT) :: B
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: DMT
!
REAL, DIMENSION(:), INTENT(IN)    :: PDIR_SW           ! incoming direct solar radiation
REAL, DIMENSION(:), INTENT(IN)    :: PSCA_SW           ! scattered incoming solar rad.
REAL, DIMENSION(:), INTENT(IN)    :: PZENITH           ! solar zenithal angle
REAL, DIMENSION(:), INTENT(IN)    :: PAZIM             ! solar azimuthal angle
!                                                      ! (radian from N, clockwise)
REAL, DIMENSION(:), INTENT(IN)    :: PFRAC_PANEL       ! Fraction of solar panel on roofs (-)
REAL, DIMENSION(:), INTENT(IN)    :: PALB_PANEL        ! Albedo     of solar panels (-)
REAL, DIMENSION(:), INTENT(IN)    :: PALB_GD       ! GD areas albedo
REAL, DIMENSION(:), INTENT(IN)    :: PSVF_GD       ! GD areas sky view factor
REAL, DIMENSION(:), INTENT(IN)    :: PALB_GRF      ! green roof albedo
REAL, DIMENSION(:), INTENT(IN)    :: PDN_RF          ! snow-covered roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDF_RF          ! snow-free    roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDN_RD          ! snow-covered road fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDF_RD          ! snow-free    road fraction
!
!new arguments for shading
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_RD     ! solar radiation received
!                                                      ! by snow-free roads
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_WL_A   ! solar radiation received
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_WL_B   ! solar radiation received
!                                                    ! by snow-free walls
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_GD     ! solar radiation received
!                                                    ! by GD areas
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_RF     ! solar radiation received
!                                                    ! by RF areas (below solar panels if any)
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_SN_RD  ! solar radiation received
!                                                    ! by snow-covered roads
REAL, DIMENSION(:), INTENT(OUT)    :: PDIR_ALB_TWN   ! town direct albedo
REAL, DIMENSION(:), INTENT(OUT)    :: PSCA_ALB_TWN   ! town diffuse albedo
!
REAL, DIMENSION(:), INTENT(OUT)    :: PSW_RAD_GD     ! solar radiation reaching GD areas
REAL, DIMENSION(:), INTENT(OUT)    :: PREC_SW_WIN      ! solar radiation received by windows

REAL, DIMENSION(:), INTENT(OUT)    :: PREF_SW_GRND     ! total solar radiation reflected by ground
REAL, DIMENSION(:), INTENT(OUT)    :: PREF_SW_FAC      ! total solar radiation reflected by wall
!new arguments for shading
REAL, DIMENSION(:), INTENT(OUT)    :: PE_SHADING       ! Energy that is not reflected 
                                                       ! by the shading, nor transmitted through
                                                       ! the bld, nor absorbed by the
                                                       ! [W/m2(win)]
LOGICAL, DIMENSION(:),INTENT(INOUT):: OSHAD_DAY        ! has shading been necessary this day ?
LOGICAL, DIMENSION(:),INTENT(IN)   :: OSHADE           ! are building conditions favorable for 
!                                                      ! shading (independantly of solar irradiance) ?
!
LOGICAL, INTENT(IN), OPTIONAL :: OALB_ONLY
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
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZAALB_RD            ! averaged albedo
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZDIR_SW_RD          ! direct radiation reaching
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZDIR_SW_WL_A        ! road, wall A,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZDIR_SW_WL_B        ! wall B,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZDIR_SW_GD          ! GD areas,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZDIR_SW_WL          ! and on average on 2 walls
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZSCA_SW_RD          ! diffuse radiation reaching
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZSCA_SW_WL          ! road, wall,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZSCA_SW_GD          ! and GD areas
!
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_RF      ! solar radiation
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_RD      ! absorbed by roofs,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_WL_A    ! road, wall A,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_WL_B    ! wall B,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_WL      ! both walls on average,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_GD      ! GD areas,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_GRF     ! green roof areas,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_PANEL     ! solar panels,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_WIN       ! window (abs+trans), and snow
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_SN_RF ! over roof, wall,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_SN_RD ! and GD areas
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_RF      ! solar radiation
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_RD      ! absorbed by roofs,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_WL      ! road, wall,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_GD    ! GD areas,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_GRF   ! green roof areas,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_PANEL     ! solar panels,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_WIN       ! window (abs+trans), and snow
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_SN_RF ! over roof and wall,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_SN_RD ! coming from diffuse rad.
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZDW                   ! difference of radiation
!                                                       ! absorbed by the 2 walls
!
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZRD                 !
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZRD_DIR             ! Road direction
!                                                     ! (radian from N, clockwise)
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZGD               !
!
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZREC_DIR_SW_WIN
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZREC_SCA_SW_WIN
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZAALB_WL
!
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZTRAN_WIN             ! solar transmittivity of windows
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_WIN              ! solar transmittivity of windows
LOGICAL, DIMENSION(SIZE(PDIR_SW)) :: G_EFF_SHAD         !TRUE if shading should be active considering actual rad
!
LOGICAL :: GALB_ONLY
INTEGER                        :: JJ                    ! loop index
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('URBAN_SOLAR_ABS',0,ZHOOK_HANDLE)
!
GALB_ONLY = .FALSE.
IF (PRESENT(OALB_ONLY)) GALB_ONLY = OALB_ONLY
!
ZDIR_SW = MAX(PDIR_SW,0.)
ZSCA_SW = MAX(PSCA_SW,0.000001)
!
ZRD_DIR = T%XROAD_DIR(:) * XPI/180.
!
DO JJ=1,SIZE(T%XROAD)
!
  IF (T%XROAD(JJ)+T%XGARDEN(JJ).NE.0.) THEN
    ZRD(JJ)   = T%XROAD(JJ) / (T%XROAD(JJ)+T%XGARDEN(JJ))
    ZGD(JJ) = T%XGARDEN(JJ) / (T%XROAD(JJ)+T%XGARDEN(JJ))
  ELSE
    ZRD(JJ)=0.
    ZGD(JJ)=0.
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
  ZABS_DIR_SW_RF   (JJ) = ZDIR_SW(JJ) * (1. - T%XALB_ROOF     (JJ)) * (1.-PFRAC_PANEL(JJ))
  ZABS_DIR_SW_SN_RF(JJ) = ZDIR_SW(JJ) * (1. - T%TSNOW_ROOF%ALB(JJ)) * (1.-PFRAC_PANEL(JJ))
  ZABS_DIR_SW_GRF  (JJ) = ZDIR_SW(JJ) * (1. - PALB_GRF        (JJ)) * (1.-PFRAC_PANEL(JJ))
  ZABS_SCA_SW_RF   (JJ) = ZSCA_SW(JJ) * (1. - T%XALB_ROOF     (JJ)) * (1.-PFRAC_PANEL(JJ))
  ZABS_SCA_SW_SN_RF(JJ) = ZSCA_SW(JJ) * (1. - T%TSNOW_ROOF%ALB(JJ)) * (1.-PFRAC_PANEL(JJ))
  ZABS_SCA_SW_GRF  (JJ) = ZSCA_SW(JJ) * (1. - PALB_GRF        (JJ)) * (1.-PFRAC_PANEL(JJ))
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
  IF (T%XBLD(JJ) .GT. 0.) THEN
!
!*      2.1    radiation coefficients
!              ----------------------
!
IF (TOP%CROAD_DIR=='UNIF') THEN
  ZTHETA0(JJ) = ASIN( MIN(ABS( 1./ZTANZEN(JJ))/T%XCAN_HW_RATIO(JJ), 1. ) )
!
!*      2.2    direct solar radiation received by roads and GARDEN areas
!               -------------------------------------------------------
!
  ZDIR_SW_RD (JJ) = (  ZDIR_SW(JJ) * 2. * ZTHETA0(JJ) / XPI                 &
                     - ZDIR_SW(JJ) * 2. * ZTANZEN(JJ) / XPI                 &
                        * T%XCAN_HW_RATIO(JJ) * (1.-COS(ZTHETA0(JJ)))  )  
!
  ZDIR_SW_GD(JJ) = ( ZDIR_SW(JJ) * 2. * ZTHETA0(JJ) / XPI                 &
                   - ZDIR_SW(JJ) * 2. * ZTANZEN(JJ) / XPI                 &
                       * T%XCAN_HW_RATIO(JJ) * (1.-COS(ZTHETA0(JJ)))  )  
ELSE
  ZDIR_SW_RD (JJ) = ZDIR_SW(JJ) * &
          MAX(0.,1.-T%XCAN_HW_RATIO(JJ)*ZTANZEN(JJ)*ABS(SIN(PAZIM(JJ)-ZRD_DIR(JJ))))
  ZDIR_SW_GD(JJ) = ZDIR_SW_RD(JJ)

END IF
!
!*      2.3    direct solar radiation received by walls
!              ----------------------------------------
!
  ZDIR_SW_WL(JJ) = (ZDIR_SW(JJ) - (ZDIR_SW_RD(JJ)*ZRD(JJ)+ZDIR_SW_GD(JJ)*ZGD(JJ))) &
                      * 0.5 / T%XCAN_HW_RATIO(JJ)  
!

  ELSE
!
   ZDIR_SW_RD (JJ) = ZDIR_SW(JJ)
   ZDIR_SW_GD (JJ) = ZDIR_SW(JJ)
   ZDIR_SW_WL (JJ) = 0.
!
  ENDIF
!
IF (TOP%CROAD_DIR=='UNIF' .OR. TOP%CWALL_OPT=='UNIF') THEN
!* if walls are averaged, then
  ZDIR_SW_WL_A(JJ) = ZDIR_SW_WL(JJ)
  ZDIR_SW_WL_B(JJ) = ZDIR_SW_WL(JJ)
ELSE
!* if walls are separated, then radiation reaches the wall facing sun
! Note that wall A is the one facing mostly to the South (depending to
! road orientation), and wall B in the one facing mostly to the North
!
! In case of N-S road, wall A is the West  wall (= East-facing  wall), 
!                  and wall B is the East  wall (= West-facing  wall)
! In case of E-W road, wall A is the North wall (= South-facing wall), 
!                  and wall B is the South wall (= North-facing wall)
  IF (SIN(PAZIM(JJ)-ZRD_DIR(JJ))>0.) THEN
    ZDIR_SW_WL_A(JJ) = 2.* ZDIR_SW_WL(JJ)
    ZDIR_SW_WL_B(JJ) = 0.
  ELSE
    ZDIR_SW_WL_A(JJ) = 0.
    ZDIR_SW_WL_B(JJ) = 2.* ZDIR_SW_WL(JJ)
  END IF
END IF
!
!
!
!*      2.4    diffuse solar radiation received by roads and GARDEN areas
!              ---------------------------------------------------------
!
  ZSCA_SW_RD (JJ) = ZSCA_SW(JJ) * T%XSVF_ROAD(JJ)
!
  ZSCA_SW_GD (JJ) = ZSCA_SW(JJ) * PSVF_GD(JJ)
!
!*      2.5    diffuse solar radiation received by walls
!              -----------------------------------------
!
  ZSCA_SW_WL (JJ) = ZSCA_SW(JJ) * T%XSVF_WALL(JJ)
!
!*      2.6    total solar radiation received by GARDEN areas
!              ---------------------------------------------
!
  PSW_RAD_GD (JJ) = ZDIR_SW_GD(JJ) + ZSCA_SW_GD(JJ)
!
!*      2.7    averaged albedos when snow is present
!              -------------------------------------
!
  ZAALB_RD   (JJ) =  PDF_RD (JJ) * T%XALB_ROAD   (JJ) + PDN_RD (JJ) * T%TSNOW_ROAD%ALB (JJ)  
!
!
ENDDO
!
!*      2.7b    averaged facade albedo
!              -------------------------------------
!
IF (TOP%CBEM=='BEM') THEN
    !
    ZTRAN_WIN(:) = B%XTRAN_WIN(:)
    !
    G_EFF_SHAD(:) = OSHADE(:).AND.(ZDIR_SW_WL(:) + ZSCA_SW_WL(:) > XWIN_SW_MAX) 
    !
    OSHAD_DAY(:)  = G_EFF_SHAD(:) .OR. OSHAD_DAY(:)
    !
    CALL WINDOW_SHADING(B%XSHGC, B%XSHGC_SH, OSHAD_DAY, T%XALB_WALL, &
                        B%XABS_WIN, ZABS_WIN, B%XALB_WIN, ZTRAN_WIN      )
    !
ELSE
  !
  ZABS_WIN (:) = 0.
  B%XALB_WIN (:) = 0.
  ZTRAN_WIN(:) = 0.
  !
ENDIF
!
ZAALB_WL(:) =  B%XGR(:) * B%XALB_WIN(:) + (1.-B%XGR(:)) * T%XALB_WALL(:)
!
!*      2.8    absorption of direct incoming solar radiation
!              ---------------------------------------------
!
!
 CALL SOLAR_REFLECTIONS(ZDIR_SW_RD, ZDIR_SW_WL, ZDIR_SW_GD, ZABS_DIR_SW_RD, &
                        ZABS_DIR_SW_SN_RD, ZABS_DIR_SW_WL, ZABS_DIR_SW_GD,  &
                        ZABS_DIR_SW_WIN    )  
!
IF (TOP%CROAD_DIR=='UNIF' .OR. TOP%CWALL_OPT=='UNIF') THEN
!* if walls are averaged, then
  ZABS_DIR_SW_WL_A = ZABS_DIR_SW_WL
  ZABS_DIR_SW_WL_B = ZABS_DIR_SW_WL
ELSE
!* if walls are separated, then radiation reaches the wall facing sun
! Note that wall A is the one facing mostly to the North (depending to
! road orientation), and wall B in the one facing mostly to the South.
  ZDW = (1.-T%XALB_WALL(:)) * ZAALB_WL(:) * (1.-2.*T%XSVF_WALL(:)) &
       / (1.+ZAALB_WL(:)*(1.-2.*T%XSVF_WALL(:)))                 &
       * 0.5 * (ZDIR_SW_WL_A(:)-ZDIR_SW_WL_B(:))            &
       + 0.5 * (1.-T%XALB_WALL(:)) * (ZDIR_SW_WL_A-ZDIR_SW_WL_B)
  ZABS_DIR_SW_WL_A = ZABS_DIR_SW_WL + ZDW
  ZABS_DIR_SW_WL_B = ZABS_DIR_SW_WL - ZDW
END IF
!
!*      2.9    absorption of diffuse incoming solar radiation
!              ----------------------------------------------
!
 CALL SOLAR_REFLECTIONS(ZSCA_SW_RD,ZSCA_SW_WL, ZSCA_SW_GD,              &
                        ZABS_SCA_SW_RD, ZABS_SCA_SW_SN_RD,              &
                        ZABS_SCA_SW_WL, ZABS_SCA_SW_GD, ZABS_SCA_SW_WIN    )  
!
! solar flux reflected for wall and road
!
PREF_SW_GRND = ZRD * T%XALB_ROAD / (1.-T%XALB_ROAD ) * (ZABS_DIR_SW_RD + ZABS_SCA_SW_RD) &
             + ZGD * PALB_GD     / (1.-PALB_GD     ) * (ZABS_DIR_SW_GD + ZABS_SCA_SW_GD)
!
PREF_SW_FAC = (1 - B%XGR) * T%XALB_WALL / (1.-T%XALB_WALL) * (ZABS_DIR_SW_WL + ZABS_SCA_SW_WL)    &
            +      B%XGR  * B%XALB_WIN  / (1 - B%XALB_WIN) * (ZABS_DIR_SW_WIN + ZABS_SCA_SW_WIN)
!
!-------------------------------------------------------------------------------
!
!*      3.     Town albedo
!              -----------
!
!*      3.1    direct albedo
!              -------------
!
 CALL TOWN_ALBEDO(ZDIR_SW,ZABS_DIR_SW_RF,ZABS_DIR_SW_SN_RF,           &
                   ZABS_DIR_SW_RD, ZABS_DIR_SW_SN_RD,ZABS_DIR_SW_WL,  &
                   ZABS_DIR_SW_GD, ZABS_DIR_SW_GRF, ZABS_DIR_SW_WIN,  &
                   ZABS_DIR_SW_PANEL, PDIR_ALB_TWN                    )  
!
!*      3.2    direct albedo
!              -------------
!
 CALL TOWN_ALBEDO(ZSCA_SW,ZABS_SCA_SW_RF,ZABS_SCA_SW_SN_RF,           &
                   ZABS_SCA_SW_RD, ZABS_SCA_SW_SN_RD,ZABS_SCA_SW_WL,  &
                   ZABS_SCA_SW_GD, ZABS_SCA_SW_GRF, ZABS_SCA_SW_WIN,  &
                   ZABS_SCA_SW_PANEL, PSCA_ALB_TWN                    )  
!
WHERE (PDIR_ALB_TWN==XUNDEF) PDIR_ALB_TWN = PSCA_ALB_TWN
!
IF (GALB_ONLY) THEN
  IF (LHOOK) CALL DR_HOOK('URBAN_SOLAR_ABS',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!-------------------------------------------------------------------------------
!
!*      4.     Trivial cases
!              -------------
!
WHERE(PDIR_SW(:)==0.)
  ZABS_DIR_SW_RF   (:) = 0.
  ZABS_DIR_SW_RD   (:) = 0.
  ZABS_DIR_SW_WL_A (:) = 0.
  ZABS_DIR_SW_WL_B (:) = 0.
  ZABS_DIR_SW_GD   (:) = 0.
  ZABS_DIR_SW_GRF  (:) = 0.
  ZABS_DIR_SW_PANEL(:) = 0.
  ZABS_DIR_SW_WIN  (:) = 0.
  ZABS_DIR_SW_SN_RF(:) = 0.
  ZABS_DIR_SW_SN_RD(:) = 0.
END WHERE
!
WHERE(PSCA_SW(:)==0.)
  ZABS_SCA_SW_RF   (:) = 0.
  ZABS_SCA_SW_RD   (:) = 0.
  ZABS_SCA_SW_WL   (:) = 0.
  ZABS_SCA_SW_GD   (:) = 0.
  ZABS_SCA_SW_GRF  (:) = 0.
  ZABS_SCA_SW_PANEL(:) = 0.
  ZABS_SCA_SW_WIN  (:) = 0.
  ZABS_SCA_SW_SN_RF(:) = 0.
  ZABS_SCA_SW_SN_RD(:) = 0.
END WHERE
!
DMT%XABS_SW_ROOF      (:) = 0.
DMT%XABS_SW_ROAD      (:) = 0.
DMT%XABS_SW_WALL_A    (:) = 0.
DMT%XABS_SW_WALL_B    (:) = 0.
DMT%XABS_SW_GARDEN    (:) = 0.
DMT%XABS_SW_GREENROOF (:) = 0.
DMT%XABS_SW_SNOW_ROOF (:) = 0.
DMT%XABS_SW_SNOW_ROAD (:) = 0.
IF (TOP%CBEM=="BEM") THEN
  DMT%XABS_SW_WIN  (:) = 0.
ENDIF
!
PREC_SW_WIN   (:) = 0.
PREC_SW_RD    (:) = 0.
PREC_SW_WL_A  (:) = 0.
PREC_SW_WL_B  (:) = 0.
PREC_SW_GD    (:) = 0.
PREC_SW_SN_RD (:) = 0.
PREC_SW_RF    (:) = 0.
!
!-------------------------------------------------------------------------------
!
!
IF (TOP%CBEM=='BEM') THEN
  !
  DO JJ=1,SIZE(T%XROAD)
    !
    ! solar radiation absorbed (but not transmitted) by windows
    !
    ZREC_DIR_SW_WIN (JJ) = ZABS_DIR_SW_WIN(JJ) / (1.-B%XALB_WIN(JJ))
    ZREC_SCA_SW_WIN (JJ) = ZABS_SCA_SW_WIN(JJ) / (1.-B%XALB_WIN(JJ))
    !
    PREC_SW_WIN     (JJ) = ZREC_DIR_SW_WIN(JJ) + ZREC_SCA_SW_WIN(JJ)
    !
    DMT%XABS_SW_WIN(JJ) = (ZREC_DIR_SW_WIN  (JJ) + ZREC_SCA_SW_WIN   (JJ)) * ZABS_WIN(JJ)
    !
    DMT%XTR_SW_WIN (JJ) = PREC_SW_WIN(JJ) * ZTRAN_WIN(JJ)
    !
  ENDDO
  !
ENDIF
!
DO JJ=1,SIZE(T%XROAD)
!
!*      5.     Total solar radiation absorbed by each surface
!              ----------------------------------------------
!
! solar radiation absorbed by roofs
!
  DMT%XABS_SW_ROOF     (JJ) = ZABS_DIR_SW_RF     (JJ) + ZABS_SCA_SW_RF     (JJ)
!
! solar radiation absorbed by roads
!
  DMT%XABS_SW_ROAD     (JJ) = ZABS_DIR_SW_RD     (JJ) + ZABS_SCA_SW_RD     (JJ)
!
! solar radiation absorbed by GARDEN areas
!
  DMT%XABS_SW_GARDEN   (JJ) = ZABS_DIR_SW_GD   (JJ) + ZABS_SCA_SW_GD   (JJ)
!
! solar radiation absorbed by GRF areas
!
  DMT%XABS_SW_GREENROOF(JJ) = ZABS_DIR_SW_GRF(JJ) + ZABS_SCA_SW_GRF(JJ)
!
! solar radiation absorbed by walls
!
  DMT%XABS_SW_WALL_A (JJ) = ZABS_DIR_SW_WL_A   (JJ) + ZABS_SCA_SW_WL     (JJ)
  DMT%XABS_SW_WALL_B (JJ) = ZABS_DIR_SW_WL_B   (JJ) + ZABS_SCA_SW_WL     (JJ)
!
!
! solar radiation absorbed by snow on roofs
!
  DMT%XABS_SW_SNOW_ROOF (JJ) = ZABS_DIR_SW_SN_RF (JJ) + ZABS_SCA_SW_SN_RF (JJ)
!
! solar radiation absorbed by snow on roads
!
  DMT%XABS_SW_SNOW_ROAD (JJ) = ZABS_DIR_SW_SN_RD (JJ) + ZABS_SCA_SW_SN_RD (JJ)
!
!-------------------------------------------------------------------------------
!
!*      6.     total solar radiation received by roads and GARDEN areas
!              -------------------------------------------------------
!
  PREC_SW_RD      (JJ) = DMT%XABS_SW_ROAD      (JJ)/(1.-T%XALB_ROAD   (JJ))
!
  PREC_SW_SN_RD   (JJ) = DMT%XABS_SW_SNOW_ROAD (JJ)/(1.-T%TSNOW_ROAD%ALB(JJ))
!
  PREC_SW_WL_A    (JJ) = DMT%XABS_SW_WALL_A    (JJ)/(1.-T%XALB_WALL   (JJ))
  PREC_SW_WL_B    (JJ) = DMT%XABS_SW_WALL_B    (JJ)/(1.-T%XALB_WALL   (JJ))
!
  PREC_SW_GD      (JJ) = DMT%XABS_SW_GARDEN    (JJ)/(1.-PALB_GD (JJ))
!
!*      6.2    total solar radiation received by roof surfaces below solar panels
!
  PREC_SW_RF      (JJ) = (PDIR_SW(JJ) + PSCA_SW(JJ)) * (1.-PFRAC_PANEL(JJ))
!
!-------------------------------------------------------------------------------
!
!*      7.     total solar radiation transmitted inside building
!*             and energy not ref., nor absorbed, nor transmitted
!              --------------------------------------------------
!
!               [W/m2(bld)]
  PE_SHADING(JJ) = PREC_SW_WIN(JJ) * (1. - B%XALB_WIN(JJ) - ZABS_WIN(JJ) - ZTRAN_WIN(JJ))
!               [W/m2(win)]

ENDDO
!
IF (TOP%LSOLAR_PANEL) THEN
  !
  DO JJ=1,SIZE(T%XROAD)
    !
    ! solar radiation absorbed by solar panels
    !
    DMT%XABS_SW_PANEL(JJ) = ZABS_DIR_SW_PANEL    (JJ) + ZABS_SCA_SW_PANEL    (JJ)
    !
  ENDDO
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('URBAN_SOLAR_ABS',1,ZHOOK_HANDLE)
CONTAINS
!
!-------------------------------------------------------------------------------
SUBROUTINE SOLAR_REFLECTIONS(ZSW_RD,ZSW_WL, ZSW_GD, ZABS_SW_RD,ZABS_SW_SN_RD, &
                             ZABS_SW_WL, ZABS_SW_GD, ZABS_SW_WIN             )  
!
REAL, DIMENSION(:), INTENT(IN) :: ZSW_RD          ! solar radiation received by road,
REAL, DIMENSION(:), INTENT(IN) :: ZSW_WL          ! wall, and GD areas 
REAL, DIMENSION(:), INTENT(IN) :: ZSW_GD        ! before reflection
REAL, DIMENSION(:), INTENT(OUT):: ZABS_SW_RD      ! solar radiation absorbed by
REAL, DIMENSION(:), INTENT(OUT):: ZABS_SW_SN_RD ! solar radiation absorbed by
REAL, DIMENSION(:), INTENT(OUT):: ZABS_SW_WL      ! road, snow over road, and wall 
REAL, DIMENSION(:), INTENT(OUT):: ZABS_SW_GD    ! solar radiation absorbed by garden
REAL, DIMENSION(:), INTENT(OUT):: ZABS_SW_WIN       ! solar radiation absorbed by window
!
REAL, DIMENSION(SIZE(ZSW_RD)) :: ZREF0_SW_RD    ! first solar reflection
REAL, DIMENSION(SIZE(ZSW_RD)) :: ZREF0_SW_WL    ! against road, wall
REAL, DIMENSION(SIZE(ZSW_RD)) :: ZREF0_SW_GD  ! and GD areas
REAL, DIMENSION(SIZE(ZSW_RD)) :: ZSREF_SW_RD    ! sum of all reflections
REAL, DIMENSION(SIZE(ZSW_RD)) :: ZSREF_SW_WL    ! against road, wall,
REAL, DIMENSION(SIZE(ZSW_RD)) :: ZSREF_SW_GD  ! and GD areas
!
REAL, DIMENSION(SIZE(ZSW_RD)) :: ZWORK1, ZWORK2, ZDENOM
INTEGER :: JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SOLAR_REFLECTIONS',0,ZHOOK_HANDLE)
!
DO JJ=1,SIZE(ZSW_RD)
!*      A.     first solar radiation reflection
!              --------------------------------
!
  ZREF0_SW_RD(JJ) = ZAALB_RD(JJ) * ZSW_RD(JJ) 
!
  ZREF0_SW_GD(JJ) = PALB_GD(JJ)  * ZSW_GD(JJ)
!
  ZREF0_SW_WL(JJ) = ZAALB_WL(JJ) * ZSW_WL(JJ)
!
!*      B.     sum of solar radiation reflected
!              --------------------------------
!

  ZDENOM(JJ) = 1. - (1.-2.*T%XSVF_WALL(JJ)) * T%XALB_WALL(JJ) - (1. - T%XSVF_ROAD(JJ))* &
               T%XSVF_WALL(JJ)*T%XALB_WALL(JJ)*ZAALB_RD(JJ)*ZRD(JJ) &
                  - (1. - PSVF_GD(JJ))* &
                T%XSVF_WALL(JJ)*T%XALB_WALL(JJ)*PALB_GD(JJ)*ZGD(JJ)  

  ZWORK1(JJ) = T%XSVF_WALL(JJ) * T%XALB_WALL(JJ) * ZRD(JJ)
  ZWORK2(JJ) = T%XSVF_WALL(JJ) * T%XALB_WALL(JJ) * ZGD(JJ)
!
!
  ZSREF_SW_WL(JJ) = (  ZREF0_SW_WL(JJ) + ZWORK1(JJ) *ZREF0_SW_RD(JJ)   &
                        + ZWORK2(JJ) *ZREF0_SW_GD(JJ)) / ZDENOM(JJ)  

  ZSREF_SW_RD(JJ) = ((1.- T%XSVF_ROAD(JJ)) * ZAALB_RD(JJ) * ZREF0_SW_WL(JJ) &
                    +(1.- T%XSVF_ROAD(JJ)) * ZAALB_RD(JJ) * ZWORK1(JJ) * ZREF0_SW_RD(JJ)  &
                    +(1.- T%XSVF_ROAD(JJ)) * ZAALB_RD(JJ) * ZWORK2(JJ) * ZREF0_SW_GD(JJ)) &
                      / ZDENOM(JJ) + ZREF0_SW_RD(JJ)   

  ZSREF_SW_GD(JJ) = ((1.- PSVF_GD(JJ)) * PALB_GD(JJ) * ZREF0_SW_WL(JJ)   &
                    +(1.- PSVF_GD(JJ)) * PALB_GD(JJ) * ZWORK1(JJ) * ZREF0_SW_RD(JJ) &
                    +(1.- PSVF_GD(JJ)) * PALB_GD(JJ) * ZWORK2(JJ) * ZREF0_SW_GD(JJ) )&
                      / ZDENOM(JJ) + ZREF0_SW_GD(JJ)  
!
!*      C.     total solar radiation received by roads and GD areas
!              -------------------------------------------------------
!
  ZABS_SW_RD(JJ)    = (1.-T%XALB_ROAD(JJ)) * (ZSW_RD(JJ) + ZSREF_SW_WL(JJ) * (1.- T%XSVF_ROAD(JJ)))  
!
  ZABS_SW_SN_RD(JJ) = (1.-T%TSNOW_ROAD%ALB(JJ)) * &
                        (ZSW_RD(JJ) + ZSREF_SW_WL(JJ) * (1.- T%XSVF_ROAD(JJ)))  
!
  ZABS_SW_GD(JJ)    = (1.-PALB_GD(JJ)) * (ZSW_GD(JJ) + ZSREF_SW_WL(JJ)  * (1.- PSVF_GD(JJ)))  
!
!
!*      D.     total solar radiation received by walls
!              ---------------------------------------
!
  ZABS_SW_WL(JJ)      = (1.-T%XALB_WALL(JJ))                                  &
                         * (ZSW_WL(JJ)                                        &
                           + ZSREF_SW_RD(JJ) *        T%XSVF_WALL(JJ)*ZRD(JJ) &
                           + ZSREF_SW_GD(JJ) *        T%XSVF_WALL(JJ)*ZGD(JJ) &
                           + ZSREF_SW_WL(JJ) * (1.-2.*T%XSVF_WALL(JJ))          )  
!
  ZABS_SW_WIN (JJ)      = (1.-B%XALB_WIN (JJ))                                  &
                         * (ZSW_WL(JJ)                                          &
                           + ZSREF_SW_RD(JJ)   *        T%XSVF_WALL(JJ)*ZRD(JJ) &
                           + ZSREF_SW_GD(JJ) *          T%XSVF_WALL(JJ)*ZGD(JJ) &
                           + ZSREF_SW_WL(JJ)   * (1.-2.*T%XSVF_WALL(JJ))          )  
!
ENDDO
!
IF (LHOOK) CALL DR_HOOK('SOLAR_REFLECTIONS',1,ZHOOK_HANDLE)
!
END SUBROUTINE SOLAR_REFLECTIONS
!
!-------------------------------------------------------------------------------
!
SUBROUTINE TOWN_ALBEDO(ZSW,ZABS_SW_RF,ZABS_SW_SN_RF,ZABS_SW_RD,ZABS_SW_SN_RD,&
                       ZABS_SW_WL,ZABS_SW_GD, ZABS_SW_GRF, ZABS_SW_WIN,      &
                         ZABS_SW_PANEL, ZALBEDO                         )  
!
REAL, DIMENSION(:), INTENT(IN) :: ZSW               ! incoming solar radiation
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_RF      ! solar radiation absorbed by roofs
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_RD      ! solar radiation absorbed by roads
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_WL      ! solar radiation absorbed by walls
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_WIN       ! solar radiation absorbed & transmitted by windows
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_GD    ! solar radiation absorbed by GARDEN areas
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_GRF ! solar radiation absorbed by green roof areas
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_SN_RF ! solar radiation absorbed by roof snow
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_SN_RD ! solar radiation absorbed by road snow
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_PANEL     ! solar radiation absorbed by solar panels
REAL, DIMENSION(:), INTENT(OUT):: ZALBEDO           ! town averaged albedo

REAL, DIMENSION(SIZE(ZSW))     :: ZSW_UP            ! outgoing solar radiation
INTEGER :: JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('TOWN_ALBEDO',0,ZHOOK_HANDLE)
DO JJ=1,SIZE(ZSW)

  ZSW_UP(JJ) = ZSW(JJ)                                                           &
            - ( T%XBLD(JJ)   *(1.-T%XGREENROOF(JJ))*PDF_RF(JJ) *ZABS_SW_RF   (JJ)&
               +T%XBLD(JJ)   *(1.-T%XGREENROOF(JJ))*PDN_RF(JJ) *ZABS_SW_SN_RF(JJ)&
               +T%XBLD(JJ)   *    T%XGREENROOF(JJ)             *ZABS_SW_GRF   (JJ)&
               +T%XBLD(JJ)   *    PFRAC_PANEL(JJ)              *ZABS_SW_PANEL (JJ)&
               +T%XROAD(JJ)                        *PDF_RD(JJ) *ZABS_SW_RD (JJ)    &
               +T%XROAD(JJ)                        *PDN_RD(JJ) *ZABS_SW_SN_RD(JJ)&
               +T%XGARDEN(JJ)                                  *ZABS_SW_GD(JJ)   &
               +T%XWALL_O_HOR(JJ)              *(1.-B%XGR(JJ)) *ZABS_SW_WL(JJ)     &
               +T%XWALL_O_HOR(JJ)              *    B%XGR(JJ)  *ZABS_SW_WIN (JJ)     )  
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
