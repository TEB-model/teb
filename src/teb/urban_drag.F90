!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE URBAN_DRAG(TOP, T, B, HIMPLICIT_WIND, PTSTEP, PT_CANYON, PQ_CANYON, &
                          PU_CANYON, PT_LOWCAN, PQ_LOWCAN, PU_LOWCAN, PZ_LOWCAN,   &
                          PTS_ROOF, PTS_ROAD, PTS_WALL, PTS_GARDEN,                &
                          PDELT_SNOW_ROOF, PDELT_SNOW_ROAD,  PEXNS, PEXNA, PTA,    &
                          PQA, PPS, PRHOA,PZREF, PUREF, PVMOD, PWS_ROOF_MAX,       &
                          PWS_ROAD_MAX, PPEW_A_COEF, PPEW_B_COEF,                  &
                          PPEW_A_COEF_LOWCAN, PPEW_B_COEF_LOWCAN, PQSAT_ROOF,      &
                          PQSAT_ROAD, PDELT_ROOF, PDELT_ROAD, PCD, PCDN, PAC_ROOF, &
                          PAC_ROOF_WAT, PAC_WALL, PAC_ROAD, PAC_ROAD_WAT, PAC_TOP, &
                          PAC_GARDEN, PRI, PUW_ROAD, PUW_ROOF, PDUWDU_ROAD,        &
                          PDUWDU_ROOF, PUSTAR_TOWN, PAC_WIN ) 
!   ##########################################################################
!
!!****  *URBAN_DRAG*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the surface drag over artificial surfaces as towns, 
!     taking into account the canyon like geometry of urbanized areas.
!         
!     
!!**  METHOD
!!    ------
!
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
!!      Original    20/01/98 
!!          01/00 (V. Masson)  separation of skimming, wake and isolated flows
!!          09/00 (V. Masson)  use of Z0 for roads
!!          12/02 (A. Lemonsu) convective speed w* in canyon
!             04 (A. Lemonsu) z0h=z0m for resistance canyon-atmosphere
!          03/08 (S. Leroyer) debug PU_CAN (1. * H/3)
!          12/08 (S. Leroyer) option (TOP%CZ0H) for z0h applied on roof, road and town
!!         09/12 B. Decharme new wind implicitation
!          11/11 (G. Pigeon) apply only urban_exch_coef when necessary if
!                            canopy/no canopy
!          09/12 (G. Pigeon) add new formulation for outdoor conv. coef for
!                            wall/roof/window
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_BEM_n, ONLY : BEM_t
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_CSTS,ONLY : XLVTT, XPI, XCPD, XG, XKARMAN
!
!USE MODE_SBLS
USE MODE_THERMOS
USE MODI_URBAN_EXCH_COEF
USE MODE_CONV_DOE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(BEM_t), INTENT(INOUT) :: B
!
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
!                                                     ! 'OLD' = direct
!                                                     ! 'NEW' = Taylor serie, order 1
!
REAL,               INTENT(IN)    :: PTSTEP         ! time-step
REAL, DIMENSION(:), INTENT(IN)    :: PT_CANYON      ! canyon air temperature
REAL, DIMENSION(:), INTENT(IN)    :: PQ_CANYON      ! canyon air specific humidity.
REAL, DIMENSION(:), INTENT(IN)    :: PU_CANYON      ! hor. wind in canyon
REAL, DIMENSION(:), INTENT(IN)    :: PU_LOWCAN     ! wind near the road
REAL, DIMENSION(:), INTENT(IN)    :: PT_LOWCAN     ! temp. near the road
REAL, DIMENSION(:), INTENT(IN)    :: PQ_LOWCAN     ! hum. near the road
REAL, DIMENSION(:), INTENT(IN)    :: PZ_LOWCAN     ! height of atm. var. near the road
REAL, DIMENSION(:), INTENT(IN)    :: PTS_ROOF       ! surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_ROAD       ! surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_WALL       ! surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GARDEN     ! surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_SNOW_ROOF! fraction of snow on roof
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_SNOW_ROAD! fraction of snow on road
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS          ! surface exner function
REAL, DIMENSION(:), INTENT(IN)    :: PTA            ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA            ! specific humidity
                                                    ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD          ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PPS            ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PEXNA          ! exner function
                                                    ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA          ! air density
REAL, DIMENSION(:), INTENT(IN)    :: PZREF          ! reference height of the first
                                                    ! atmospheric level (temperature)
REAL, DIMENSION(:), INTENT(IN)    :: PUREF          ! reference height of the first
                                                    ! atmospheric level (wind)
REAL, DIMENSION(:), INTENT(IN)    :: PWS_ROOF_MAX   ! maximum deepness of roof
REAL, DIMENSION(:), INTENT(IN)    :: PWS_ROAD_MAX   ! and water reservoirs (kg/m2)
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_A_COEF    ! implicit coefficients (m2s/kg)
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_B_COEF    ! for wind coupling     (m/s)
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_A_COEF_LOWCAN ! implicit coefficients for wind coupling (m2s/kg)
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_B_COEF_LOWCAN ! between low canyon wind and road (m/s)
!
REAL, DIMENSION(:), INTENT(OUT)   :: PQSAT_ROOF     ! qsat(Ts)
REAL, DIMENSION(:), INTENT(OUT)   :: PQSAT_ROAD     ! qsat(Ts)
REAL, DIMENSION(:), INTENT(OUT)   :: PDELT_ROOF     ! water fraction on
REAL, DIMENSION(:), INTENT(OUT)   :: PDELT_ROAD     ! snow-free surfaces
REAL, DIMENSION(:), INTENT(OUT)   :: PCD            ! drag coefficient
REAL, DIMENSION(:), INTENT(OUT)   :: PCDN           ! neutral drag coefficient
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROOF       ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROOF_WAT   ! aerodynamical conductance (for water)
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_WALL       ! aerodynamical conductance
!                                                   ! between canyon air and
!                                                   ! walls 
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROAD       ! aerodynamical conductance
!                                                   ! between canyon air and
!                                                   ! roads
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROAD_WAT   ! aerodynamical conductance
!                                                   ! between canyon air and
!                                                   ! road (for water)
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_TOP        ! aerodynamical conductance
!                                                   ! between canyon top and atm.
REAL, DIMENSION(:), INTENT(IN)    :: PAC_GARDEN     ! aerodynamical conductance
!                                                   ! between canyon air and GARDEN areas
REAL, DIMENSION(:), INTENT(OUT)   :: PRI            ! Town Richardson number
!
REAL, DIMENSION(:), INTENT(OUT)   :: PUW_ROAD       ! Momentum flux for roads
REAL, DIMENSION(:), INTENT(OUT)   :: PUW_ROOF       ! Momentum flux for roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PDUWDU_ROAD    ! 
REAL, DIMENSION(:), INTENT(OUT)   :: PDUWDU_ROOF    ! 
REAL, DIMENSION(:), INTENT(OUT)   :: PUSTAR_TOWN    ! Fraction velocity for town
!
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_WIN        ! aerodynamical conductance for window
!
!*      0.2    declarations of local variables
!
!
REAL, DIMENSION(SIZE(PTA)) :: ZTS_TOWN     ! town averaged temp.
REAL, DIMENSION(SIZE(PTA)) :: ZQ_TOWN      ! town averaged hum.
REAL, DIMENSION(SIZE(PTA)) :: ZAVDELT_ROOF ! averaged water frac.
REAL, DIMENSION(SIZE(PTA)) :: ZQ_ROOF      ! roof spec. hum.
REAL, DIMENSION(SIZE(PTA)) :: ZZ0_ROOF     ! roof roughness length
REAL, DIMENSION(SIZE(PTA)) :: ZZ0_ROAD     ! road roughness length
REAL, DIMENSION(SIZE(PTA)) :: ZW_CAN       ! ver. wind in canyon
REAL, DIMENSION(SIZE(PTA)) :: ZRI          ! Richardson number
REAL, DIMENSION(SIZE(PTA)) :: ZLE_MAX      ! maximum latent heat flux available
REAL, DIMENSION(SIZE(PTA)) :: ZLE          ! actual latent heat flux
REAL, DIMENSION(SIZE(PTA)) :: ZRA_ROOF     ! aerodynamical resistance
REAL, DIMENSION(SIZE(PTA)) :: ZCH_ROOF     ! drag coefficient for heat
REAL, DIMENSION(SIZE(PTA)) :: ZRA_TOP      ! aerodynamical resistance
REAL, DIMENSION(SIZE(PTA)) :: ZCH_TOP      ! drag coefficient for heat
REAL, DIMENSION(SIZE(PTA)) :: ZRA_ROAD     ! aerodynamical resistance
REAL, DIMENSION(SIZE(PTA)) :: ZCH_ROAD     ! drag coeifficient for heat
REAL, DIMENSION(SIZE(PTA)) :: ZCD_ROAD     ! road  surf. exchange coefficient
REAL, DIMENSION(SIZE(PTA)) :: ZAC          ! town aerodynamical conductance (not used)
REAL, DIMENSION(SIZE(PTA)) :: ZRA          ! town aerodynamical resistance  (not used)
REAL, DIMENSION(SIZE(PTA)) :: ZCH          ! town drag coefficient for heat (not used)
REAL, DIMENSION(SIZE(PTA)) :: ZCD          ! any surf. exchange coefficient (not used)
REAL, DIMENSION(SIZE(PTA)) :: ZCDN         ! any surf. neutral exch. coef.  (not used)
!
REAL, DIMENSION(SIZE(PTA)) :: ZU_STAR, ZW_STAR !! 
REAL, DIMENSION(SIZE(PTA)) :: ZQ0              !! 
!
REAL, DIMENSION(SIZE(PTA)) :: ZUSTAR2      ! square of friction velocity (m2/s2
REAL, DIMENSION(SIZE(PTA)) :: ZVMOD        ! module of the horizontal wind at t+1
!
! for calculation of momentum fluxes
REAL, DIMENSION(SIZE(PTA)) :: ZLMO         ! Monin-Obukhov length
REAL, DIMENSION(SIZE(PTA)) :: ZUSTAR_ROAD  ! friction velocity for roads
REAL, DIMENSION(SIZE(PTA)) :: ZUSTAR_ROOF  ! friction velocity for roofs
REAL, DIMENSION(SIZE(PTA)) :: ZUSTAR_TOWN  !
!
REAL, DIMENSION(SIZE(PTA)) :: ZZ0_TOP      ! roughness length for zac_top calculation
REAL, DIMENSION(SIZE(PTA)) :: ZCHTCN_WIN   ! natural convective heat transfer coef. for window [W/(m2.K)]
REAL, DIMENSION(SIZE(PTA)) :: ZCHTCN_ROOF  ! natural convective heat transfer coef. for roof [W/(m2.K)]
REAL, DIMENSION(SIZE(PTA)) :: ZCHTCS_ROOF  ! forced convective heat transfer coef. for smooth roof [W/(m2.K)]
REAL, DIMENSION(SIZE(PTA)) :: ZCHTCN_WALL  ! natural convective heat transfer coef. for wall [W/(m2.K)]
REAL, DIMENSION(SIZE(PTA)) :: ZCHTCS_WALL  ! forced natural convective heat transfer coef. for smooth wall [W/(m2.K)]
!
INTEGER                   ::  JLOOP, JJ            !! 
!
REAL :: ZZ0_O_Z0H = 200.  ! z0/z0h ratio used in Mascart (1995) formulation.
!                         ! It is set to the maximum value acceptable by
!                         ! formulation. Observed values are often larger in cities.
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('URBAN_DRAG',0,ZHOOK_HANDLE)
!
ZZ0_ROOF(:)    = 0.15                      ! z0 for roofs
ZZ0_ROAD(:)    = MIN(0.05,0.1*PZ_LOWCAN(:))! z0 for roads
!
ZZ0_TOP(:) = T%XZ0_TOWN(:)
!
PCD    (:) = XUNDEF
PCDN   (:) = XUNDEF
PAC_TOP(:) = XUNDEF
PRI    (:) = XUNDEF
!
PUW_ROAD   (:) = XUNDEF
PUW_ROOF   (:) = XUNDEF
PDUWDU_ROAD(:) = XUNDEF
PDUWDU_ROOF(:) = XUNDEF
PUSTAR_TOWN(:) = XUNDEF
!
!-------------------------------------------------------------------------------
!
!*      1.     roof and road saturation specific humidity
!              ------------------------------------------
!
PQSAT_ROOF(:) = QSAT (PTS_ROOF(:), PPS(:))
!
PQSAT_ROAD(:) = QSAT (PTS_ROAD(:), PPS(:))
!
!-------------------------------------------------------------------------------
!
!*      2.     fraction of water on roofs
!              --------------------------
!
PDELT_ROOF=1.
!
!*      2.1    general case
!              ------------
!
WHERE (PQSAT_ROOF(:) >= PQA(:) )
  PDELT_ROOF(:) = (T%XWS_ROOF(:)/PWS_ROOF_MAX)**(2./3.)
END WHERE
!
!*      2.2    dew deposition on roofs (PDELT_ROOF=1)
!              -----------------------
!
!-------------------------------------------------------------------------------
!
!*      3.     fraction of water on roads
!              --------------------------
!
PDELT_ROAD=1.
!
!*      3.1    general case
!              ------------
!
WHERE (PQSAT_ROAD(:) >= PQ_CANYON(:) )
  PDELT_ROAD(:) = (T%XWS_ROAD(:)/PWS_ROAD_MAX)**(2./3.)
END WHERE
!
!*      3.2    dew deposition on roads (PDELT_ROAD=1)
!              -----------------------
!
!-------------------------------------------------------------------------------
!
!*      4.     Drag coefficient for momentum between roof level and atmosphere
!              ---------------------------------------------------------------
!
!
DO JJ=1,SIZE(PTA)
!
!*      4.1    Averaged temperature at roof level
!              ----------------------------------
!
  ZTS_TOWN(JJ) = T%XBLD(JJ) * PTS_ROOF(JJ) + (1.-T%XBLD(JJ)) * PT_CANYON(JJ)
!
!*      4.2    Averaged water fraction on roofs
!              -------------------------------
!
  ZAVDELT_ROOF(JJ) = PDELT_ROOF(JJ) * PDELT_SNOW_ROOF(JJ)
!
!*      4.3    Roof specific humidity
!              ----------------------
!
  ZQ_ROOF(JJ) = PQSAT_ROOF(JJ) * ZAVDELT_ROOF(JJ)
!
!*      4.4    Averaged Saturation specific humidity
!              -------------------------------------
!
  ZQ_TOWN(JJ) =  T%XBLD(JJ) * ZQ_ROOF(JJ) + (1.-T%XBLD(JJ)) * PQ_CANYON(JJ)
!
ENDDO
!
!-------------------------------------------------------------------------------
!
!*      5.     Momentum drag coefficient
!              -------------------------
!
IF (.NOT. TOP%LCANOPY) THEN
  CALL URBAN_EXCH_COEF(TOP%CZ0H, ZZ0_O_Z0H, ZTS_TOWN, ZQ_TOWN, PEXNS, PEXNA, PTA, PQA,     &
                       PZREF+ T%XBLD_HEIGHT/3., PUREF+T%XBLD_HEIGHT/3., PVMOD, T%XZ0_TOWN, &
                       PRI, PCD, PCDN, ZAC, ZRA, ZCH                              )
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      6.     Drag coefficient for heat fluxes between roofs and atmosphere
!              -------------------------------------------------------------
!
IF (TOP%CCH_BEM == "DOE-2") THEN
   ZCHTCN_ROOF = CHTC_UP_DOE(PTS_ROOF, PTA)
   ZCHTCS_ROOF = CHTC_SMOOTH_WIND_DOE(ZCHTCN_ROOF, PVMOD)
   PAC_ROOF = CHTC_ROUGH_DOE(ZCHTCN_ROOF, ZCHTCS_ROOF, T%XROUGH_ROOF) / PRHOA / XCPD
ELSE
   CALL URBAN_EXCH_COEF(TOP%CZ0H, ZZ0_O_Z0H, PTS_ROOF, ZQ_ROOF, PEXNS, PEXNA, PTA, PQA, &
                        PZREF, PUREF, PVMOD, ZZ0_ROOF, ZRI, ZCD, ZCDN, PAC_ROOF,        &
                        ZRA_ROOF, ZCH_ROOF               )
ENDIF
!
!
DO JJ=1,SIZE(PTA)
  ZLE_MAX(JJ)     = T%XWS_ROOF(JJ) / PTSTEP * XLVTT
  ZLE    (JJ)     =(PQSAT_ROOF(JJ) - PQA(JJ)) &
                 * PAC_ROOF(JJ) * PDELT_ROOF(JJ) * XLVTT * PRHOA(JJ)
!
  PAC_ROOF_WAT(JJ) = PAC_ROOF(JJ)
!
  IF (PDELT_ROOF(JJ)==0.) PAC_ROOF_WAT(JJ)=0.
!
  IF (ZLE(JJ)>0.) PAC_ROOF_WAT(JJ) = PAC_ROOF(JJ) * MIN ( 1. , ZLE_MAX(JJ)/ZLE(JJ) )
!
ENDDO
!-------------------------------------------------------------------------------
!
!*      7.     Drag coefficient for heat fluxes between canyon and atmosphere
!              --------------------------------------------------------------
!
!* Because air/air exchanges are considered, roughness length for heat is set
!  equal to roughness length for momentum.
!
IF (.NOT. TOP%LCANOPY) THEN
  CALL URBAN_EXCH_COEF('MASC95', 1., PT_CANYON, PQ_CANYON, PEXNS, PEXNA, PTA, PQA,    &
                        PZREF+T%XBLD_HEIGHT-PZ_LOWCAN, PUREF+T%XBLD_HEIGHT-PZ_LOWCAN, &
                        PVMOD, ZZ0_TOP,  ZRI, ZCD, ZCDN, PAC_TOP, ZRA_TOP, ZCH_TOP  )
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      8.     Drag coefficient for heat fluxes between walls, road and canyon
!              ---------------------------------------------------------------
!
!*      8.1    aerodynamical conductance for walls
!              -----------------------------------
!
IF (TOP%CCH_BEM == "DOE-2") THEN
  DO JJ=1,SIZE(PTA)
    ZCHTCN_WALL(JJ) = CHTC_VERT_DOE(PTS_WALL(JJ), PT_CANYON(JJ))
    ZCHTCS_WALL(JJ) = 0.5 * (CHTC_SMOOTH_LEE_DOE (ZCHTCN_WALL(JJ), PU_CANYON(JJ)) + &
                             CHTC_SMOOTH_WIND_DOE(ZCHTCN_WALL(JJ), PU_CANYON(JJ)) )
                      
    PAC_WALL(JJ) = CHTC_ROUGH_DOE(ZCHTCN_WALL(JJ), ZCHTCS_WALL(JJ), T%XROUGH_WALL(JJ)) / XCPD / PRHOA(JJ)
  END DO
ELSE
  PAC_WALL(:) = ( 11.8 + 4.2 * PU_CANYON(:) ) / XCPD / PRHOA(:)
END IF
!
!*      8.2    aerodynamical conductance for roads
!              -----------------------------------
!
ZW_STAR(:) = 0.
ZQ0(:)     = 0.
!
!
DO JLOOP=1,3
 !
  ZW_CAN(:)   = ZW_STAR(:)
  !
  !
  CALL URBAN_EXCH_COEF(TOP%CZ0H, ZZ0_O_Z0H, PTS_ROAD, PQ_LOWCAN, PEXNS, PEXNA,  &
                       PT_LOWCAN, PQ_LOWCAN, PZ_LOWCAN, PZ_LOWCAN,              &
                       PU_LOWCAN+ZW_CAN, ZZ0_ROAD, ZRI, ZCD_ROAD, ZCDN,         &
                       PAC_ROAD, ZRA_ROAD, ZCH_ROAD   )
  !
  DO JJ=1,SIZE(PTA)

    ZQ0(JJ)     = (PTS_WALL  (JJ) - PT_CANYON(JJ)) * PAC_WALL  (JJ) * T%XWALL_O_GRND(JJ)

    IF (T%XROAD(JJ) .GT. 0.) THEN
      ZQ0(JJ)   = ZQ0(JJ) &
            + (PTS_ROAD  (JJ) - PT_LOWCAN(JJ)) * PAC_ROAD  (JJ) * T%XROAD  (JJ)/(T%XROAD(JJ)+T%XGARDEN(JJ)) 
    ENDIF
    IF (T%XGARDEN(JJ) .GT. 0.) THEN
      ZQ0(JJ)   = ZQ0(JJ) &
            + (PTS_GARDEN(JJ) - PT_LOWCAN(JJ)) * PAC_GARDEN(JJ) * T%XGARDEN(JJ)/(T%XROAD(JJ)+T%XGARDEN(JJ))
    ENDIF    
    !
    IF (ZQ0(JJ) >= 0.) THEN
      ZW_STAR(JJ) = ( (XG * PEXNA(JJ) / PTA(JJ)) * ZQ0(JJ) * T%XBLD_HEIGHT(JJ)) ** (1/3.)
    ELSE
      ZW_STAR(JJ) = 0.
    ENDIF
!
  ENDDO
!
END DO
!
!
!*      8.4    aerodynamical conductance for water limited by available water
!              --------------------------------------------------------------
!
DO JJ=1,SIZE(PTA)
  !
  ZLE_MAX(JJ)     = T%XWS_ROAD(JJ) / PTSTEP * XLVTT
  ZLE    (JJ)     = ( PQSAT_ROAD(JJ) - PQ_LOWCAN(JJ) )                   &
                   *   PAC_ROAD(JJ) * PDELT_ROAD(JJ) * XLVTT * PRHOA(JJ)
  !
  PAC_ROAD_WAT(JJ) = PAC_ROAD(JJ)
  !
  IF (PDELT_ROAD(JJ)==0.) PAC_ROAD_WAT(JJ) = 0.
  !
  IF (ZLE(JJ)>0.) PAC_ROAD_WAT(JJ) = PAC_ROAD(JJ) * MIN ( 1. , ZLE_MAX(JJ)/ZLE(JJ) )
  !
  !
  !*      8.5    aerodynamical conductance for window
  !              ------------------------------------
  !
  ZCHTCN_WIN(JJ) = CHTC_VERT_DOE(B%XT_WIN1(JJ), PT_CANYON(JJ))
  !
  PAC_WIN(JJ) = 0.5 * (CHTC_SMOOTH_LEE_DOE(ZCHTCN_WIN(JJ), PU_CANYON(JJ)) + &
                      CHTC_SMOOTH_WIND_DOE(ZCHTCN_WIN(JJ), PU_CANYON(JJ)) ) &
                   / PRHOA(JJ) / XCPD
  !
  !-------------------------------------------------------------------------------
  !
  !*      9.     Momentum fluxes
  !              ---------------
  !
  !*      9.1    For roads
  !              ---------
  !
  !* road friction
  !
  IF (TOP%LCANOPY) THEN
    !
    ZUSTAR2(JJ)=XUNDEF
    !
    IF(HIMPLICIT_WIND=='OLD')THEN
      !   old implicitation
      ZUSTAR2(JJ) = (ZCD_ROAD(JJ)*PU_LOWCAN(JJ)*PPEW_B_COEF_LOWCAN(JJ))/              &
                    (1.0-PRHOA(JJ)*ZCD_ROAD(JJ)*PU_LOWCAN(JJ)*PPEW_A_COEF_LOWCAN(JJ))
    ELSE
      !   new implicitation
      ZUSTAR2(JJ) = (ZCD_ROAD(JJ)*PU_LOWCAN(JJ)*(2.*PPEW_B_COEF_LOWCAN(JJ)-PU_LOWCAN(JJ)))/  &
                    (1.0-2.0*PRHOA(JJ)*ZCD_ROAD(JJ)*PU_LOWCAN(JJ)*PPEW_A_COEF_LOWCAN(JJ))
      !                   
      ZVMOD(JJ) = PRHOA(JJ)*PPEW_A_COEF_LOWCAN(JJ)*ZUSTAR2(JJ) + PPEW_B_COEF_LOWCAN(JJ)
      ZVMOD(JJ) = MAX(ZVMOD(JJ),0.)
      !
      IF(PPEW_A_COEF_LOWCAN(JJ)/= 0.)THEN
        ZUSTAR2(JJ) = MAX( ( ZVMOD(JJ) - PPEW_B_COEF_LOWCAN(JJ) ) / (PRHOA(JJ)*PPEW_A_COEF_LOWCAN(JJ)), 0.)
      ENDIF
      !              
    ENDIF
    !
    PUW_ROAD(JJ) = - ZUSTAR2(JJ)
    !
    PDUWDU_ROAD(JJ) = 0. ! implicitation already taken into account in PUW_ROAD
    !
    !*      9.2    For roofs
    !              ---------
    !
    !* roof friction
    !* neutral case, as guess
    !
    !
    ZUSTAR_ROOF(JJ) = PVMOD(JJ) * XKARMAN / LOG(PZREF(JJ)/ZZ0_ROOF(JJ))
    !
    PUW_ROOF(JJ)    = - ZUSTAR_ROOF(JJ)**2
    PDUWDU_ROOF(JJ) = 0.
    IF (PVMOD(JJ)/=0.) PDUWDU_ROOF(JJ) = 2. * PUW_ROOF(JJ) / PVMOD(JJ)
    !
  ELSE
    !
    !*      9.3    For town
    !              --------
    !
    ZUSTAR2(JJ)=XUNDEF
    !  
    IF(HIMPLICIT_WIND=='OLD')THEN
      !   old implicitation
      ZUSTAR2(JJ) = (PCD(JJ)*PVMOD(JJ)*PPEW_B_COEF(JJ))/            &
                    (1.0-PRHOA(JJ)*PCD(JJ)*PVMOD(JJ)*PPEW_A_COEF(JJ))
    ELSE
      !   new implicitation
      ZUSTAR2(JJ) = (PCD(JJ)*PVMOD(JJ)*(2.*PPEW_B_COEF(JJ)-PVMOD(JJ)))/ &
                    (1.0-2.0*PRHOA(JJ)*PCD(JJ)*PVMOD(JJ)*PPEW_A_COEF(JJ)) 
      !                   
      ZVMOD(JJ) = PRHOA(JJ)*PPEW_A_COEF(JJ)*ZUSTAR2(JJ) + PPEW_B_COEF(JJ)
      ZVMOD(JJ) = MAX(ZVMOD(JJ),0.)
      !
      IF(PPEW_A_COEF(JJ)/= 0.)THEN
        ZUSTAR2(JJ) = MAX( ( ZVMOD(JJ) - PPEW_B_COEF(JJ) ) / (PRHOA(JJ)*PPEW_A_COEF(JJ)), 0.)
      ENDIF
      !                        
    ENDIF
    !
    PUSTAR_TOWN(JJ) = SQRT(ZUSTAR2(JJ))
    !
  ENDIF
  !
ENDDO
!
IF (LHOOK) CALL DR_HOOK('URBAN_DRAG',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE URBAN_DRAG
