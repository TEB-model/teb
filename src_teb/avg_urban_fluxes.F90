!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE AVG_URBAN_FLUXES(                                                &
                       PTS_TOWN, PEMIS_TOWN,                                    &
                       PT_CANYON, PQ_CANYON,                                    &
                       PT_LOWCAN, PQ_LOWCAN,                                    &
                       PTS_ROOF,PTS_ROAD,PTS_WALL_A,PTS_WALL_B,PTS_GARDEN,      &
                       PTA, PQA, PRHOA, PPS,                                    &
                       PH_TRAFFIC, PLE_TRAFFIC, PH_INDUSTRY, PLE_INDUSTRY,      &
                       PBLD, PROAD, PGARDEN, PWALL_O_HOR, PWALL_O_GRND,         &
                       PFRAC_GR,                                                &
                       PEMIS_ROOF, PESNOW_ROOF, PEMIS_GREENROOF,                &
                       PLW_RAD,                                                 &
                       PABS_LW_ROOF, PABS_LW_WALL_A, PABS_LW_WALL_B,            &
                       PABS_LW_ROAD, PABS_LW_GARDEN, PABS_LW_GREENROOF,         &
                       PABS_LW_SNOW_ROOF, PABS_LW_SNOW_ROAD,                    &
                       PAC_ROOF, PAC_ROOF_WAT,                                  &
                       PAC_WALL, PAC_ROAD, PAC_ROAD_WAT, PAC_TOP,               &
                       PAC_GARDEN,                                              &
                       PQSAT_GARDEN, PAC_AGG_GARDEN, PHU_AGG_GARDEN,            &
                       PQSAT_ROOF, PQSAT_ROAD,                                  &
                       PDELT_ROOF, PDELT_ROAD,                                  &
                       PROOF_FRAC, PWALL_FRAC, PROAD_FRAC, PGARDEN_FRAC,        &
                       PTOTS_O_HORS,                                            &
                       PDF_ROOF, PDN_ROOF, PDF_ROAD, PDN_ROAD,                  &
                       PRN_ROOF, PH_ROOF, PLE_ROOF, PGFLUX_ROOF,                &
                       PRN_ROAD, PH_ROAD, PLE_ROAD, PGFLUX_ROAD,                &
                       PRN_GARDEN, PH_GARDEN, PLE_GARDEN, PGFLUX_GARDEN,        &
                       PRN_WALL_A, PH_WALL_A, PLE_WALL_A, PGFLUX_WALL_A,        &
                       PRN_WALL_B, PH_WALL_B, PLE_WALL_B, PGFLUX_WALL_B,        &
                       PLEW_ROOF, PLESNOW_ROOF,                                 &
                       PLEW_ROAD, PLESNOW_ROAD, PHSNOW_ROAD,                    &
                       PEVAP_GARDEN, PEVAP_GREENROOF,                           &
                       PRN_GRND, PH_GRND, PLE_GRND, PGFLUX_GRND,                &
                       PRN_TOWN, PH_TOWN, PLE_TOWN, PGFLUX_TOWN, PEVAP_TOWN,    &
                       PRUNOFF_GARDEN,PRUNOFF_ROAD,PRUNOFF_STRLROOF,            &
                       PRUNOFF_GREENROOF, PDRAIN_GREENROOF, PRUNOFF_TOWN,       &
                       PABS_LW_PANEL, PEMIS_PANEL, PFRAC_PANEL, PRN_PANEL,      &
                       PH_PANEL,                                                &
                       PH_WASTE, PLE_WASTE, PF_WASTE_CAN,                       &
                       PABS_LW_WIN, PT_WIN1, PGR, PEMIT_LW_ROAD,                &
                       PEMIT_LW_GARDEN, PEMIT_LW_GRND, HBEM,                    &
                       PSVF_ROAD, PSVF_GARDEN, PSVF_WALL, PGARDEN_O_GRND,       &
                       PROAD_O_GRND,                                            &
                       PEMIS_ROAD, PESNOW_ROAD, PEMIS_WALL, PEMIS_GARDEN,       &
                       OCANOPY                                                  )  
!   ##########################################################################
!
!!****  *AVG_URBAN_FLUXES* computes fluxes on urbanized surfaces  
!!
!!    PURPOSE
!!    -------
!         
!     
!!**  METHOD
!     ------
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
!!                     12/02 (A. Lemonsu) modifications of emissivity and Tstown
!!                     07/07 (P. LeMoigne) expression of latent heat fluxes as 
!!                           functions of w'theta' instead of w'T'
!!                     17/10 (G. Pigeon)  computation of anthropogenic heat due
!!                           to domestic heating
!!                     10/11 (G. Pigeon) simplification for road, garden, roof,
!!                           wall fractions
!!                     08/13 (V. Masson) adds solar panels
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XCPD, XLVTT, XLSTT, XSTEFAN
!
USE MODE_THERMOS
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
REAL, DIMENSION(:), INTENT(OUT)   :: PTS_TOWN          ! town surface temperature
REAL, DIMENSION(:), INTENT(OUT)   :: PEMIS_TOWN        ! town equivalent emissivity
REAL, DIMENSION(:), INTENT(INOUT) :: PT_CANYON         ! canyon air temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PQ_CANYON         ! canyon air specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PT_LOWCAN         ! low canyon air temperature
REAL, DIMENSION(:), INTENT(IN)    :: PQ_LOWCAN         ! low canyon air specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PTS_ROOF          ! roof surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_ROAD          ! road surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_WALL_A        ! wall surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_WALL_B        ! wall surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GARDEN        ! green area surface temperature

REAL, DIMENSION(:), INTENT(IN)    :: PTA               ! temperature at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PQA               ! specific humidity
                                                       ! at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA             ! air density
                                                       ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PPS               ! surface pressure
REAL, DIMENSION(:), INTENT(IN)    :: PH_TRAFFIC        ! anthropogenic sensible
!                                                      ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PLE_TRAFFIC       ! anthropogenic latent
!                                                      ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PH_INDUSTRY       ! anthropogenic sensible
!                                                      ! heat fluxes due to factories
REAL, DIMENSION(:), INTENT(IN)    :: PLE_INDUSTRY      ! anthropogenic latent
!                                                      ! heat fluxes due to factories
REAL, DIMENSION(:), INTENT(IN)    :: PBLD              ! fraction of buildings
REAL, DIMENSION(:), INTENT(IN)    :: PROAD             ! fraction of roads
REAL, DIMENSION(:), INTENT(IN)    :: PGARDEN           ! fraction of green areas
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_O_HOR       ! wall Surf. / (bld+road+green) Surf.
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_O_GRND      ! wall Surf. / ground (road+green) Surf.
REAL, DIMENSION(:), INTENT(IN)    :: PFRAC_GR          ! fraction of green roofs
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_ROOF        ! roof emissivity
! 
REAL, DIMENSION(:), INTENT(IN)    :: PESNOW_ROOF       ! snow roof emissivity
!
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_GREENROOF   ! green roof emissivity
!
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD           ! incoming longwave rad.
!
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_ROOF      ! absorbed LW rad. by roof
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_WALL_A    ! absorbed LW rad. by wall
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_WALL_B    ! absorbed LW rad. by wall
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_ROAD      ! absorbed LW rad. by road
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_GARDEN    ! absorbed LW rad. by green areas
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_GREENROOF ! absorbed LW rad. by green roofs
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_SNOW_ROOF ! abs. LW rad. by snow
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_SNOW_ROAD ! abs. LW rad. by snow
!
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROOF          ! surface conductance
!                                                      ! for heat transfers
!                                                      ! above roofs
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROOF_WAT      ! surface conductance
!                                                      ! for heat transfers
!                                                      ! above roofs (for water)
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WALL          ! surface conductance
!                                                      ! for heat transfer
!                                                      ! between wall and canyon
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROAD          ! surface conductance
!                                                      ! for heat transfers
!                                                      ! between road and canyon
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROAD_WAT      ! surface conductance
!                                                      ! for heat transfers
!                                                      ! inside canyon (for water)
REAL, DIMENSION(:), INTENT(IN)    :: PAC_TOP           ! aerodynamical conductance
!                                                      ! between atmosphere and
!                                                      ! canyon top
REAL, DIMENSION(:), INTENT(IN)    :: PAC_GARDEN        ! aerodynamical conductance
!                                                      ! between atmosphere and
!                                                      ! green areas
REAL, DIMENSION(:), INTENT(IN)    :: PQSAT_GARDEN      ! q_sat(Ts)
REAL, DIMENSION(:), INTENT(IN)    :: PAC_AGG_GARDEN    ! aggregated aerodyn resistance for green areas
REAL, DIMENSION(:), INTENT(IN)    :: PHU_AGG_GARDEN    ! aggregated relative humidity for green areas
!
REAL, DIMENSION(:), INTENT(IN)    :: PQSAT_ROOF        ! q_sat(Ts)
REAL, DIMENSION(:), INTENT(IN)    :: PQSAT_ROAD        ! q_sat(Ts)
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_ROOF        ! water fraction on snow-free
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_ROAD        ! roof and roads
REAL, DIMENSION(:), INTENT(IN)    :: PROOF_FRAC        ! roof, wall,
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_FRAC        ! road, and green area
REAL, DIMENSION(:), INTENT(IN)    :: PROAD_FRAC        ! fractions
REAL, DIMENSION(:), INTENT(IN)    :: PGARDEN_FRAC      ! of exchange surf.
REAL, DIMENSION(:), INTENT(IN)    :: PTOTS_O_HORS      ! total canyon+roof surface
!                                                      ! over horizontal surface
REAL, DIMENSION(:), INTENT(IN)    :: PDF_ROOF          ! snow-free    roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROOF          ! snow-covered roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDF_ROAD          ! snow-free    road fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROAD          ! snow-covered road fraction
!
REAL, DIMENSION(:), INTENT(IN)    :: PRN_ROOF          ! net radiation over roof
REAL, DIMENSION(:), INTENT(IN)    :: PH_ROOF           ! sensible heat flux over roof
REAL, DIMENSION(:), INTENT(IN)    :: PLE_ROOF          ! latent heat flux over roof
REAL, DIMENSION(:), INTENT(IN)    :: PGFLUX_ROOF       ! flux through the roof
REAL, DIMENSION(:), INTENT(IN)    :: PRN_ROAD          ! net radiation over road
REAL, DIMENSION(:), INTENT(IN)    :: PH_ROAD           ! sensible heat flux over road
REAL, DIMENSION(:), INTENT(IN)    :: PLE_ROAD          ! latent heat flux over road
REAL, DIMENSION(:), INTENT(IN)    :: PGFLUX_ROAD       ! flux through the road
REAL, DIMENSION(:), INTENT(IN)    :: PRN_GARDEN        ! net radiation over green areas
REAL, DIMENSION(:), INTENT(IN)    :: PH_GARDEN         ! sensible heat flux over green areas
REAL, DIMENSION(:), INTENT(IN)    :: PLE_GARDEN        ! latent heat flux over green areas
REAL, DIMENSION(:), INTENT(IN)    :: PGFLUX_GARDEN     ! flux through the green areas
REAL, DIMENSION(:), INTENT(IN)    :: PRN_WALL_A        ! net radiation over wall
REAL, DIMENSION(:), INTENT(IN)    :: PH_WALL_A         ! sensible heat flux over wall
REAL, DIMENSION(:), INTENT(IN)    :: PLE_WALL_A        ! latent heat flux over wall
REAL, DIMENSION(:), INTENT(IN)    :: PGFLUX_WALL_A     ! flux through the wall
REAL, DIMENSION(:), INTENT(IN)    :: PRN_WALL_B        ! net radiation over wall
REAL, DIMENSION(:), INTENT(IN)    :: PH_WALL_B         ! sensible heat flux over wall
REAL, DIMENSION(:), INTENT(IN)    :: PLE_WALL_B        ! latent heat flux over wall
REAL, DIMENSION(:), INTENT(IN)    :: PGFLUX_WALL_B     ! flux through the wall
!
REAL, DIMENSION(:), INTENT(IN)    :: PLEW_ROOF         ! latent heat flux of snowfree roof
REAL, DIMENSION(:), INTENT(IN)    :: PLESNOW_ROOF      ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PLEW_ROAD         ! latent heat flux of snowfree road
REAL, DIMENSION(:), INTENT(IN)    :: PLESNOW_ROAD      ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PHSNOW_ROAD       ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PEVAP_GARDEN      ! evaporation over gardens
REAL, DIMENSION(:), INTENT(IN)    :: PEVAP_GREENROOF   ! evaporation over green roofs
!
REAL, DIMENSION(:), INTENT(IN)    :: PRUNOFF_GARDEN    ! surface runoff over green areas      (kg/m2/s)
REAL, DIMENSION(:), INTENT(IN)    :: PRUNOFF_ROAD      ! surface runoff over roads            (kg/m2/s)
REAL, DIMENSION(:), INTENT(IN)    :: PRUNOFF_STRLROOF  ! surface runoff over structural roofs (kg/m2/s)
REAL, DIMENSION(:), INTENT(IN)    :: PRUNOFF_GREENROOF ! surface runoff over green roofs      (kg/m2/s)
REAL, DIMENSION(:), INTENT(IN)    :: PDRAIN_GREENROOF  ! outlet drainage at green roof base   (kg/m2/s)
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_GRND          ! net radiation over ground
REAL, DIMENSION(:), INTENT(OUT)   :: PH_GRND           ! sensible heat flux over ground
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_GRND          ! latent heat flux over ground
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_GRND       ! flux through the ground
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_TOWN          ! net radiation over town
REAL, DIMENSION(:), INTENT(OUT)   :: PH_TOWN           ! sensible heat flux over town
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_TOWN          ! latent heat flux over town
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_TOWN       ! flux through the ground for town
REAL, DIMENSION(:), INTENT(OUT)   :: PEVAP_TOWN        ! evaporation (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT)   :: PRUNOFF_TOWN      ! aggregated runoff for town (kg/m2/s)

REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_PANEL     ! absorbed LW radiation by solar panels
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_PANEL       ! emissivity of solar panels
REAL, DIMENSION(:), INTENT(IN)    :: PFRAC_PANEL       ! fraction of solar panels on roofs
REAL, DIMENSION(:), INTENT(IN)    :: PRN_PANEL         ! net radiation of solar panels
REAL, DIMENSION(:), INTENT(IN)    :: PH_PANEL          ! sensible heat flux of solar panels
!
REAL, DIMENSION(:), INTENT(IN)    :: PH_WASTE     ! sensible waste heat released by HVAC systems
REAL, DIMENSION(:), INTENT(IN)    :: PLE_WASTE    ! latent waste heat released by HVAC systems
REAL, DIMENSION(:), INTENT(IN)    :: PF_WASTE_CAN ! fraction of waste heat released into the canyon
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_WIN  ! absorbed LW radiation by windows
REAL, DIMENSION(:), INTENT(IN)    :: PT_WIN1      ! window outdoor temperature
REAL, DIMENSION(:), INTENT(IN)    :: PGR          ! glazing ratio
!
REAL, DIMENSION(:), INTENT(IN)    :: PEMIT_LW_ROAD ! LW emitted by the road (W/m2 road)
REAL, DIMENSION(:), INTENT(IN)    :: PEMIT_LW_GARDEN ! LW emitted by the garden (W/m2 garden)
REAL, DIMENSION(:), INTENT(OUT)   :: PEMIT_LW_GRND ! LW emitted by the ground (road+garden) (W/m2 ground)
 CHARACTER(LEN=3), INTENT(IN)      :: HBEM ! Building Energy model 'DEF' or 'BEM'
REAL, DIMENSION(:), INTENT(IN)    :: PGARDEN_O_GRND! garden surf. / (road+garden surf.) 
REAL, DIMENSION(:), INTENT(IN)    :: PROAD_O_GRND  ! road surf.   / (road+garden surf.) 

REAL, DIMENSION(:), INTENT(IN)    :: PSVF_ROAD     ! road sky view factor
REAL, DIMENSION(:), INTENT(IN)    :: PSVF_GARDEN   ! garden sky view factor
REAL, DIMENSION(:), INTENT(IN)    :: PSVF_WALL     ! wall sky view factor
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_ROAD    ! snow free road emissivity
REAL, DIMENSION(:), INTENT(IN)    :: PESNOW_ROAD   ! snow (road) emissivity
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_WALL    ! wall emissivity
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_GARDEN  ! garden emissivity
LOGICAL,            INTENT(IN)    :: OCANOPY ! is canopy active ?
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PLW_RAD))    :: ZLW_UP            ! upwards radiations
REAL, DIMENSION(SIZE(PROAD)) :: ZQSAT_CANYON
REAL, DIMENSION(SIZE(PROAD)) :: ZROAD, ZGARDEN
!
INTEGER :: JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('AVG_URBAN_FLUXES',0,ZHOOK_HANDLE)
!
ZROAD(:)=0.
ZGARDEN(:)=0.
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
!*      1.     Averaged fluxes for ground (green areas + road)
!              -----------------------------------------------
!
  PRN_GRND(JJ)    = ZROAD  (JJ) * PRN_ROAD  (JJ)  +  ZGARDEN(JJ) * PRN_GARDEN(JJ)
!
  PH_GRND (JJ)    = ZROAD  (JJ) * PH_ROAD   (JJ)  +  ZGARDEN(JJ) * PH_GARDEN (JJ)       &
                 +  PH_TRAFFIC(JJ) /  (1.-PBLD (JJ))
!
  PLE_GRND(JJ)    = ZROAD  (JJ) * PLE_ROAD  (JJ)  +  ZGARDEN(JJ) * PLE_GARDEN(JJ)       &
                 +  PLE_TRAFFIC(JJ) / (1.-PBLD (JJ))
!
  PGFLUX_GRND(JJ) = ZROAD  (JJ) * PGFLUX_ROAD  (JJ)  +  ZGARDEN(JJ) * PGFLUX_GARDEN(JJ)
!
!
  PEMIT_LW_GRND(JJ) = PROAD_O_GRND(JJ) * PEMIT_LW_ROAD(JJ) + PGARDEN_O_GRND(JJ) * PEMIT_LW_GARDEN(JJ)
!-------------------------------------------------------------------------------
!
!*      2.     Averaged fluxes JJ built + green areas
!              -------------------------------------
!
  PRN_TOWN(JJ)    = PTOTS_O_HORS(JJ) * ( &
                   PROOF_FRAC  (JJ) * PRN_PANEL   (JJ) * PFRAC_PANEL(JJ)&
                 + PROOF_FRAC  (JJ) * PRN_ROOF    (JJ)                  &
                 + PROAD_FRAC  (JJ) * PRN_ROAD    (JJ)                  &
                 + PGARDEN_FRAC(JJ) * PRN_GARDEN  (JJ)                  &
                 + PWALL_FRAC  (JJ) * PRN_WALL_A  (JJ) * 0.5            &
                 + PWALL_FRAC  (JJ) * PRN_WALL_B  (JJ) * 0.5 )  
!
  PH_TOWN (JJ)    = PTOTS_O_HORS(JJ) * ( &
                   PROOF_FRAC  (JJ) * PH_PANEL   (JJ) * PFRAC_PANEL(JJ) &
                 + PROOF_FRAC  (JJ) * PH_ROOF    (JJ)                   &
                 + PROAD_FRAC  (JJ) * PH_ROAD    (JJ)                   &
                 + PGARDEN_FRAC(JJ) * PH_GARDEN  (JJ)                   &
                 + PWALL_FRAC  (JJ) * PH_WALL_A  (JJ) * 0.5             &
                 + PWALL_FRAC  (JJ) * PH_WALL_B  (JJ) * 0.5 )           &
                 + PH_TRAFFIC(JJ)                                       &
                 + PH_INDUSTRY(JJ)  
!
  PLE_TOWN(JJ)    = PTOTS_O_HORS(JJ) * ( &
                   PROOF_FRAC  (JJ) * PLE_ROOF  (JJ)                     &
                 + PROAD_FRAC  (JJ) * PLE_ROAD  (JJ)                     &
                 + PGARDEN_FRAC(JJ) * PLE_GARDEN(JJ)                     &
                 + PWALL_FRAC  (JJ) * PLE_WALL_A(JJ) * 0.5               &
                 + PWALL_FRAC  (JJ) * PLE_WALL_B(JJ) * 0.5 )             &
                 + PLE_TRAFFIC (JJ)                                      &
                 + PLE_INDUSTRY(JJ)  
!
  PGFLUX_TOWN(JJ)= PTOTS_O_HORS(JJ) * ( &
                  PROOF_FRAC  (JJ) * PGFLUX_ROOF  (JJ)                   &
                + PROAD_FRAC  (JJ) * PGFLUX_ROAD  (JJ)                   &
                + PGARDEN_FRAC(JJ) * PGFLUX_GARDEN(JJ)                   &
                + PWALL_FRAC  (JJ) * PGFLUX_WALL_A(JJ) * 0.5             &
                + PWALL_FRAC  (JJ) * PGFLUX_WALL_B(JJ) * 0.5 )   
!
!-------------------------------------------------------------------------------
!
!*      3.     Infra-red Radiative properties
!              ------------------------------
!
!*      3.1    Upward IR radiation for town
!              ----------------------------
!
  ZLW_UP(JJ) = PLW_RAD      (JJ)                                                        &
            - ( PROOF_FRAC  (JJ)*(1.-PFRAC_GR(JJ))*PDF_ROOF (JJ)*PABS_LW_ROOF      (JJ) &
               +PROOF_FRAC  (JJ)*(1.-PFRAC_GR(JJ))*PDN_ROOF (JJ)*PABS_LW_SNOW_ROOF (JJ) &
               +PROOF_FRAC  (JJ)*    PFRAC_GR(JJ)               *PABS_LW_GREENROOF (JJ) &
               +PROOF_FRAC  (JJ)*    PFRAC_PANEL(JJ)            *PABS_LW_PANEL     (JJ) &
               +PROAD_FRAC  (JJ)                  *PDF_ROAD (JJ)*PABS_LW_ROAD      (JJ) &
               +PROAD_FRAC  (JJ)                  *PDN_ROAD (JJ)*PABS_LW_SNOW_ROAD (JJ) &
               +PGARDEN_FRAC(JJ)                                *PABS_LW_GARDEN    (JJ) &
               +PWALL_FRAC  (JJ)* 0.5                           *PABS_LW_WALL_A    (JJ) &
               +PWALL_FRAC  (JJ)* 0.5                           *PABS_LW_WALL_B    (JJ) &
              )*PTOTS_O_HORS(JJ)  
!
!*      3.2    Town emissivity
!              ---------------
!
    PEMIS_TOWN(JJ) =  PBLD       (JJ)*(1.-PFRAC_GR   (JJ))*PDF_ROOF(JJ)*PEMIS_ROOF     (JJ)* (1.-PFRAC_PANEL(JJ)) &
                    + PBLD       (JJ)*(1.-PFRAC_GR   (JJ))*PDN_ROOF(JJ)*PESNOW_ROOF    (JJ)* (1.-PFRAC_PANEL(JJ)) &
                    + PBLD       (JJ)*    PFRAC_GR   (JJ)              *PEMIS_GREENROOF(JJ)* (1.-PFRAC_PANEL(JJ)) &
                    + PBLD       (JJ)                                  *PEMIS_PANEL    (JJ)*     PFRAC_PANEL(JJ)  &
                    + PROAD      (JJ)*(   PSVF_ROAD  (JJ) *PDF_ROAD(JJ)*PEMIS_ROAD     (JJ)                       &
                                      +   PSVF_ROAD  (JJ) *PDN_ROAD(JJ)*PESNOW_ROAD    (JJ))                      &
                    + PGARDEN    (JJ)*    PSVF_GARDEN(JJ)              *PEMIS_GARDEN   (JJ)                       &
                    + PWALL_O_HOR(JJ)*    PSVF_WALL  (JJ)              *PEMIS_WALL     (JJ) 
!*      3.3    Town radiative surface temperature
!              ----------------------------------
!
  PTS_TOWN(JJ)   = ((ZLW_UP(JJ) - PLW_RAD(JJ)*(1.-PEMIS_TOWN(JJ))) /PEMIS_TOWN(JJ)/XSTEFAN)**0.25
!
!-------------------------------------------------------------------------------
!
!*      4.     Averaged evaporative flux (kg/m2/s)
!              -----------------------------------
!
  PEVAP_TOWN(JJ) = PTOTS_O_HORS(JJ)*(                                                           &
                   PROOF_FRAC  (JJ)*PDF_ROOF(JJ)*(1.-PFRAC_GR(JJ))*PLEW_ROOF       (JJ) /XLVTT  &
                 + PROOF_FRAC  (JJ)*PDN_ROOF(JJ)*(1.-PFRAC_GR(JJ))*PLESNOW_ROOF    (JJ) /XLSTT  &
                 + PROOF_FRAC  (JJ)             *    PFRAC_GR(JJ) *PEVAP_GREENROOF (JJ)         &
                 + PROAD_FRAC  (JJ)*PDF_ROAD(JJ)                  *PLEW_ROAD       (JJ) /XLVTT  &
                 + PROAD_FRAC  (JJ)*PDN_ROAD(JJ)                  *PLESNOW_ROAD    (JJ) /XLSTT  &
                 + PGARDEN_FRAC(JJ)                               *PEVAP_GARDEN    (JJ)         &
                 + PWALL_FRAC  (JJ)*0.5               * (PLE_WALL_A(JJ)+PLE_WALL_B (JJ))/XLVTT) &
                 + PROOF_FRAC  (JJ)         *(1.-PF_WASTE_CAN(JJ))*PLE_WASTE       (JJ) /XLVTT  &
                 +                                                     PLE_TRAFFIC (JJ) /XLVTT  &
                 +                                                     PLE_INDUSTRY(JJ) /XLVTT
!
!-------------------------------------------------------------------------------
!
!*      5.     Averaged runoff flux (kg/m2/s)
!              -----------------------------------
  PRUNOFF_TOWN(JJ) =  ((1.-PFRAC_GR(JJ))* PRUNOFF_STRLROOF (JJ)               &
                      +    PFRAC_GR(JJ) *(PRUNOFF_GREENROOF(JJ)+PDRAIN_GREENROOF(JJ))) * PBLD(JJ)   &
                    +     PROAD    (JJ) * PRUNOFF_ROAD     (JJ)               &
                    +     PGARDEN  (JJ) * PRUNOFF_GARDEN   (JJ)                 
!-------------------------------------------------------------------------------
!
!*      6.    Air canyon temperature at time t+dt
!             -----------------------------------
!
  IF (.NOT. OCANOPY) THEN
    PT_CANYON(JJ) =                                                                      &
           (  PTS_ROAD    (JJ) * PAC_ROAD (JJ) * PDF_ROAD (JJ) * ZROAD (JJ)              &
            + PTS_GARDEN  (JJ) * PAC_GARDEN(JJ)                * ZGARDEN(JJ)             &
            + PTS_WALL_A  (JJ) * PAC_WALL (JJ) * (1.-PGR(JJ))   * PWALL_O_GRND(JJ) * 0.5 &
            + PTS_WALL_B  (JJ) * PAC_WALL (JJ) * (1.-PGR(JJ))   * PWALL_O_GRND(JJ) * 0.5 &
            + PT_WIN1     (JJ) * PAC_WALL (JJ) *     PGR(JJ)    * PWALL_O_GRND(JJ)       &
            + PTA         (JJ) * PAC_TOP  (JJ)                                           &
            + PH_TRAFFIC  (JJ) / (1.-PBLD (JJ)) / PRHOA(JJ) / XCPD                       &
            + PHSNOW_ROAD (JJ) * PDN_ROAD (JJ)  / PRHOA(JJ) / XCPD                       &
            + PH_WASTE(JJ) * PF_WASTE_CAN(JJ)/(1-PBLD(JJ))/ PRHOA(JJ) / XCPD )           &
         / (                    PAC_ROAD (JJ) * PDF_ROAD (JJ) * ZROAD (JJ)               &
            +                   PAC_GARDEN(JJ)                * ZGARDEN(JJ)              &
            +                   PAC_WALL (JJ) * PWALL_O_GRND(JJ)                         &
            +                   PAC_TOP  (JJ)                     )  
!
!-------------------------------------------------------------------------------
!
!*      7.     Air canyon specific humidity
!              ----------------------------
!
  PQ_CANYON(JJ) = &
      (  PQSAT_ROAD   (JJ) * PAC_ROAD_WAT (JJ) * PDF_ROAD(JJ) * ZROAD (JJ) * PDELT_ROAD(JJ) &
       + PQSAT_GARDEN  (JJ) * PAC_AGG_GARDEN(JJ) * PHU_AGG_GARDEN(JJ) * ZGARDEN(JJ)         &
       + PQA          (JJ) * PAC_TOP(JJ)                                                    &
       + PLE_TRAFFIC  (JJ) / (1.-PBLD(JJ)) / PRHOA(JJ) / XLVTT                              &
       + PLESNOW_ROAD (JJ) * PDN_ROAD(JJ)  / PRHOA(JJ) / XLVTT     * ZROAD (JJ)             &
       + PLE_WASTE(JJ) * PF_WASTE_CAN(JJ) / (1-PBLD(JJ)) / PRHOA(JJ) / XLVTT       )        &
    / (  PAC_ROAD_WAT (JJ) * PDF_ROAD(JJ) * PDELT_ROAD(JJ)         * ZROAD (JJ)             &
       + PAC_AGG_GARDEN(JJ) * PHU_AGG_GARDEN(JJ)                    * ZGARDEN(JJ)           &
       + PAC_TOP(JJ)                                                               )
  ENDIF
!
ENDDO
!
!*check that canyon specif. humidity lower than saturation
!
ZQSAT_CANYON(:) = QSAT (PT_CANYON(:), PPS(:))
WHERE (PQ_CANYON(:) > ZQSAT_CANYON(:))
   PQ_CANYON(:) = ZQSAT_CANYON(:)
END WHERE
IF (LHOOK) CALL DR_HOOK('AVG_URBAN_FLUXES',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVG_URBAN_FLUXES
