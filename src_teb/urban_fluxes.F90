!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     #########
    SUBROUTINE URBAN_FLUXES(HIMPLICIT_WIND, OCANOPY, PT_CANYON,                 &
                       PPEW_A_COEF, PPEW_B_COEF,                                &
                       PEXNS, PRHOA, PVMOD,                                     &
                       PH_TRAFFIC, PLE_TRAFFIC, PH_INDUSTRY, PLE_INDUSTRY,      &
                       PBLD, PROAD, PWALL_O_HOR, PFRAC_GR,                      &
                       PABS_SW_ROOF, PABS_LW_ROOF,                              &
                       PABS_SW_WALL_A, PABS_LW_WALL_A,                          &
                       PABS_SW_WALL_B, PABS_LW_WALL_B,                          &
                       PABS_SW_ROAD, PABS_LW_ROAD,                              &
                       PAC_WALL, PCD,                                           &
                       PDF_ROOF, PDN_ROOF, PDF_ROAD, PDN_ROAD,                  &
                       PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,    &
                       PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,    & 
                       PMELT_ROOF, PDQS_ROOF,                                   &
                       PMELT_ROAD, PDQS_ROAD,                                   &
                       PDQS_WALL_A, PDQS_WALL_B,                                &
                       PFLX_BLD_ROOF, PFLX_BLD_WALL_A, PFLX_BLD_WALL_B,         &
                       PFLX_BLD_FLOOR,                                          &
                       PFLX_BLD_MASS, PE_SHADING, HBEM, PGR, PT_WIN1, PH_WASTE, &
                       PLE_WASTE, PF_WASTE_CAN, PABS_SW_WIN, PABS_LW_WIN, PQIN, &
                       PN_FLOOR, PLEW_ROOF, PLEW_ROAD,                          &
                       PRN_GREENROOF, PH_GREENROOF, PLE_GREENROOF,              &
                       PGFLUX_GREENROOF,                                        &
                       PH_ROOF, PH_ROAD, PH_WALL_A, PH_WALL_B,                  &
                       PRN_ROOF, PLE_ROOF,                                      &
                       PGFLUX_ROOF, PRN_ROAD, PLE_ROAD, PGFLUX_ROAD,            &
                       PRN_WALL_A, PLE_WALL_A, PGFLUX_WALL_A,                   &
                       PRN_WALL_B, PLE_WALL_B, PGFLUX_WALL_B,                   &
                       PRN_STRLROOF, PH_STRLROOF, PLE_STRLROOF,                 & 
                       PGFLUX_STRLROOF,                                         &
                       PRN_BLT,  PH_BLT,  PLE_BLT,                              &
                       PGFLUX_BLT, PMELT_BLT, PQF_BLD, PFLX_BLD, PDQS_BLT,      &
                       PQF_TOWN, PUSTAR_TOWN, PHVAC_COOL, PHVAC_HEAT            )
!   ##########################################################################
!
!!****  *URBAN_FLUXES* computes fluxes on urbanized surfaces  
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
!!	V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    23/01/98 
!!                     12/02 (A. Lemonsu) modifications of emissivity and Tstown
!!                     07/07 (P. LeMoigne) expression of latent heat fluxes as 
!!                           functions of w'theta' instead of w'T'
!!                     17/10 (G. Pigeon)  computation of anthropogenic heat due
!!                            to domestic heating
!!      Modified    09/2012 : B. Decharme New wind implicitation
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_CSTS,     ONLY : XCPD
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
!                                                     ! 'OLD' = direct
!                                                     ! 'NEW' = Taylor serie, order 1
LOGICAL,            INTENT(IN)    :: OCANOPY      ! is canopy active ?
REAL, DIMENSION(:), INTENT(IN)    :: PT_CANYON    ! canyon air temperature

REAL, DIMENSION(:), INTENT(IN)    :: PPEW_A_COEF  ! implicit coefficients
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_B_COEF  ! for wind coupling
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS        ! surface exner function

REAL, DIMENSION(:), INTENT(IN)    :: PRHOA        ! air density
                                                  ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD        ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PH_TRAFFIC   ! anthropogenic sensible
!                                                 ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PLE_TRAFFIC  ! anthropogenic latent
!                                                 ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PH_INDUSTRY  ! anthropogenic sensible
!                                                 ! heat fluxes due to factories
REAL, DIMENSION(:), INTENT(IN)    :: PLE_INDUSTRY ! anthropogenic latent
!                                                 ! heat fluxes due to factories
REAL, DIMENSION(:), INTENT(IN)    :: PBLD         ! fraction of buildings
REAL, DIMENSION(:), INTENT(IN)    :: PROAD        ! fraction of roads
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_O_HOR  ! wall Surf. / (bld+road+green) Surf.
REAL, DIMENSION(:), INTENT(IN)    :: PFRAC_GR     ! Fraction of green roofs
! 
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROOF ! absorbed SW rad. by roof
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_ROOF ! absorbed LW rad. by roof
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WALL_A! absorbed SW rad. by wall A
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_WALL_A! absorbed LW rad. by wall A
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WALL_B! absorbed SW rad. by wall B
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_WALL_B! absorbed LW rad. by wall B
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROAD ! absorbed SW rad. by road
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_ROAD ! absorbed LW rad. by road
!
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WALL     ! surface conductance
!                                                 ! for heat transfer
!                                                 ! between wall and canyon
REAL, DIMENSION(:), INTENT(IN)    :: PCD          ! drag coefficient
REAL, DIMENSION(:), INTENT(IN)    :: PDF_ROOF     ! snow-free    roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROOF     ! snow-covered roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDF_ROAD     ! snow-free    road fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROAD     ! snow-covered road fraction
!
REAL, DIMENSION(:), INTENT(IN)    :: PRNSNOW_ROOF ! net radiation over snow
REAL, DIMENSION(:), INTENT(IN)    :: PHSNOW_ROOF  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PLESNOW_ROOF ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PGSNOW_ROOF  ! flux under the snow
REAL, DIMENSION(:), INTENT(IN)    :: PRNSNOW_ROAD ! net radiation over snow
REAL, DIMENSION(:), INTENT(IN)    :: PHSNOW_ROAD  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PLESNOW_ROAD ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PGSNOW_ROAD  ! flux under the snow
!
REAL, DIMENSION(:), INTENT(IN)    :: PMELT_ROOF   ! snow melting on roof
REAL, DIMENSION(:), INTENT(IN)    :: PDQS_ROOF    ! storage inside roofs
REAL, DIMENSION(:), INTENT(IN)    :: PMELT_ROAD   ! snow melting on road
REAL, DIMENSION(:), INTENT(IN)    :: PDQS_ROAD    ! storage inside roads
REAL, DIMENSION(:), INTENT(IN)    :: PDQS_WALL_A  ! storage inside walls
REAL, DIMENSION(:), INTENT(IN)    :: PDQS_WALL_B  ! storage inside walls
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_ROOF! heat flx from bld to roof
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_WALL_A! heat flx from bld to wall A
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_WALL_B! heat flx from bld to wall B
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_FLOOR! heat flx from bld to floor
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_MASS! heat flx from bld to mass
REAL, DIMENSION(:), INTENT(IN)    :: PE_SHADING   ! energy not ref., nor absorbed, nor
                                                  !trans. by glazing [W
                                                  !m-2(win)]
 CHARACTER(LEN=3), INTENT(IN)      :: HBEM         ! Building Energy model 'DEF' or 'BEM'
!
REAL, DIMENSION(:), INTENT(IN)    :: PGR          ! glazing ratio
REAL, DIMENSION(:), INTENT(IN)    :: PT_WIN1      ! Outdoor window temperature [K]
REAL, DIMENSION(:), INTENT(IN)    :: PH_WASTE     ! sensible waste heat released by HVAC systems [W m-2(urb)]
REAL, DIMENSION(:), INTENT(IN)    :: PLE_WASTE    ! latent waste heat released by HVAC systems [W m-2(urb)]
REAL, DIMENSION(:), INTENT(IN)    :: PF_WASTE_CAN ! fraction of waste heat released into the canyon
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WIN  ! absorbed SW radiation by windows [W m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_WIN  ! absorbed LW radiation by windows [W m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PQIN         ! Internal heat gains [W m-2(floor)]
REAL, DIMENSION(:), INTENT(IN)    :: PN_FLOOR     ! Number of Building floors [-]
!
REAL, DIMENSION(:), INTENT(IN)   :: PLEW_ROOF     ! latent heat flux over snow-free roof
REAL, DIMENSION(:), INTENT(IN)   :: PLEW_ROAD     ! latent heat flux of snow-free road
!
REAL, DIMENSION(:), INTENT(IN)   :: PRN_GREENROOF ! net radiation of green roofs
REAL, DIMENSION(:), INTENT(IN)   :: PH_GREENROOF  ! sensible heat flux of green roofs
REAL, DIMENSION(:), INTENT(IN)   :: PLE_GREENROOF ! latent heat flux of green roofs
REAL, DIMENSION(:), INTENT(IN)   :: PGFLUX_GREENROOF ! heat storage flux of green roofs
!
REAL, DIMENSION(:), INTENT(INOUT) :: PH_ROOF      ! sensible heat flux over roof
REAL, DIMENSION(:), INTENT(INOUT) :: PH_ROAD      ! sensible heat flux from snow free road
                                                  ! to air (in) and from all road to air (out)
REAL, DIMENSION(:), INTENT(INOUT) :: PH_WALL_A    ! sensible heat flux wall to
                                                  ! air (in) and facade to air (out)     
REAL, DIMENSION(:), INTENT(INOUT) :: PH_WALL_B    ! sensible heat flux wall to
                                                  ! air (in) and facade to air (out)     
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_ROOF     ! net radiation over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_ROOF     ! latent heat flux over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_ROOF  ! flux through the roof
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_ROAD     ! net radiation over road
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_ROAD     ! latent heat flux over road

REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_ROAD  ! flux through the road
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_WALL_A   ! net radiation over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WALL_A   ! latent heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_WALL_A! flux through the wall
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_WALL_B   ! net radiation over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WALL_B   ! latent heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_WALL_B! flux through the wall
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_STRLROOF  ! net radiation of structural roof
REAL, DIMENSION(:), INTENT(OUT)   :: PH_STRLROOF   ! sensible heat flux of structural roof
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_STRLROOF  ! latent heat flux of structural roof
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_STRLROOF ! heat storage flux of structural roof
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_BLT      ! net radiation over built covers
REAL, DIMENSION(:), INTENT(OUT)   :: PH_BLT       ! sensible heat flux over built covers
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_BLT      ! latent heat flux over built covers
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_BLT   ! flux through the built covers
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_BLT    ! snow melting for town
!
REAL, DIMENSION(:), INTENT(OUT)   :: PQF_BLD      ! anthropogenic heat due to domestic heating [W m-2(bld)]
REAL, DIMENSION(:), INTENT(OUT)   :: PFLX_BLD     ! heat flx from bld to its structure         [W m-2(bld)]
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_BLT     ! storage inside urban materials
REAL, DIMENSION(:), INTENT(OUT)   :: PQF_TOWN     ! total anthropogenic heat                   [W m-2(urb)]
REAL, DIMENSION(:), INTENT(OUT)   :: PUSTAR_TOWN  ! friction velocity over town
!
REAL, DIMENSION(:), INTENT(IN)   :: PHVAC_COOL   ! Energy consumption of the cooling system
                                                  ! [W m-2(bld)]
REAL, DIMENSION(:), INTENT(IN)   :: PHVAC_HEAT   ! Energy consumption of the heating system
                                                  ! [W m-2(bld)]
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PRHOA))  :: ZH_ROOF_SNOWFREE
REAL, DIMENSION(SIZE(PRHOA))  :: ZRN_ROOF_SNOWFREE
REAL, DIMENSION(SIZE(PRHOA))  :: ZUSTAR2 ! square of friction velocity (m2/s2)
REAL, DIMENSION(SIZE(PRHOA))  :: ZVMOD   ! Wind
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('URBAN_FLUXES',0,ZHOOK_HANDLE)
!
!*      1.     Fluxes at snow-free roofs
!              -------------------------
!
!                                            net radiation
!
!
ZRN_ROOF_SNOWFREE(:) = PABS_SW_ROOF(:) + PABS_LW_ROOF(:)
!
!                                            sensible heat flux
!
ZH_ROOF_SNOWFREE(:) = PH_ROOF(:)
!
!-------------------------------------------------------------------------------
!
!*      2.     Fluxes at snow-free roads
!              -------------------------
!
!                                            net radiation
!
PRN_ROAD(:) = PABS_SW_ROAD(:) + PABS_LW_ROAD(:)
!
!-------------------------------------------------------------------------------
!
!*      3.     Fluxes at walls
!              ---------------
!
!                                            net radiation
!
PRN_WALL_A(:) = (PABS_SW_WALL_A(:) + PABS_LW_WALL_A(:)) * (1.-PGR(:))  &
              + (PABS_SW_WIN (:) + PABS_LW_WIN (:)) * PGR(:)
PRN_WALL_B(:) = (PABS_SW_WALL_B(:) + PABS_LW_WALL_B(:)) * (1.-PGR(:))  &
              + (PABS_SW_WIN (:) + PABS_LW_WIN (:)) * PGR(:)
!
!                                            heat flux into the ground
!
PGFLUX_WALL_A(:) = PRN_WALL_A(:) - PH_WALL_A(:)            
PGFLUX_WALL_B(:) = PRN_WALL_B(:) - PH_WALL_B(:)            
!
!                                            sensible heat flux
!
!before -> PH_WALL in [W.m-2(wall)]
PH_WALL_A(:) = (1. - PGR(:)) * PH_WALL_A(:) + PGR(:) * &
  ( PAC_WALL(:)*XCPD*PRHOA(:)/PEXNS(:) * (PT_WIN1(:)-PT_CANYON(:)) + PE_SHADING(:) ) &
  + PF_WASTE_CAN(:) * PH_WASTE(:) / PWALL_O_HOR(:) 
PH_WALL_B(:) = (1. - PGR(:)) * PH_WALL_B(:) + PGR(:) * &
  ( PAC_WALL(:)*XCPD*PRHOA(:)/PEXNS(:) * (PT_WIN1(:)-PT_CANYON(:)) + PE_SHADING(:) ) &
  + PF_WASTE_CAN(:) * PH_WASTE(:) / PWALL_O_HOR(:) 
!after PH_WALL in [W.m-2(facade=wall + win)]
PLE_WALL_A(:) = PF_WASTE_CAN(:) * PLE_WASTE(:) / PWALL_O_HOR(:)
PLE_WALL_B(:) = PF_WASTE_CAN(:) * PLE_WASTE(:) / PWALL_O_HOR(:)
!
!-------------------------------------------------------------------------------
!
!*      4.     Snow-free and snow-covered surfaces averaging
!              ---------------------------------------------
!
!*      4.1    Roads
!              -----
!
!                                            heat flux into the ground
!
PGFLUX_ROAD (:) =  PDF_ROAD(:) * (PRN_ROAD(:) - PH_ROAD (:) - PLEW_ROAD(:) ) &
                   + PDN_ROAD(:) * PGSNOW_ROAD(:)  
!
!                                            net radiation
!
PRN_ROAD(:) = PRN_ROAD(:) * PDF_ROAD(:) + PRNSNOW_ROAD(:) * PDN_ROAD(:)
!
!                                            sensible heat flux
!                                            total latent heat of evaporation from
!                                            the road (snow free + snow)
!
! sensible heat flux
!
PH_ROAD  (:) = PH_ROAD (:) * PDF_ROAD(:) + PHSNOW_ROAD(:) * PDN_ROAD(:)
! total latent heat of evaporation from  the road (snow free + snow)
!
PLE_ROAD (:) = PLEW_ROAD(:) * PDF_ROAD(:) + PLESNOW_ROAD(:) * PDN_ROAD(:)
!
!*      4.2    Roofs
!              -----
!
!                                            heat flux into the ground
!
PGFLUX_STRLROOF(:) =  PDF_ROOF(:) * (ZRN_ROOF_SNOWFREE(:) - ZH_ROOF_SNOWFREE(:) - PLEW_ROOF(:)) &
                    + PDN_ROOF(:) *  PGSNOW_ROOF(:)  
PGFLUX_ROOF    (:) = (1.-PFRAC_GR(:)) * PGFLUX_STRLROOF(:) + PFRAC_GR(:) * PGFLUX_GREENROOF(:)
!
!
!                                            net radiation
! 
PRN_STRLROOF   (:) = ZRN_ROOF_SNOWFREE(:) * PDF_ROOF(:) + PRNSNOW_ROOF(:) * PDN_ROOF(:)
PRN_ROOF       (:) = (1.-PFRAC_GR(:)) * PRN_STRLROOF(:) + PFRAC_GR(:) * PRN_GREENROOF(:)
!
!                                            sensible heat flux
!                                            total latent heat of evaporation from
!                                            the roof (snow free + snow)
!
! sensible heat flux
!
PH_STRLROOF    (:) = PH_ROOF(:) * PDF_ROOF(:) + PHSNOW_ROOF(:)  * PDN_ROOF(:)
PH_ROOF        (:) = (1.-PFRAC_GR(:)) * PH_STRLROOF(:) + PFRAC_GR(:) * PH_GREENROOF(:) &
                   + (1 - PF_WASTE_CAN(:)) * PH_WASTE(:) / PBLD(:)
!
! total latent heat of evaporation from the roof (snow free + snow)
!
PLE_STRLROOF   (:) = PLEW_ROOF(:) * PDF_ROOF(:)  + PLESNOW_ROOF(:) * PDN_ROOF(:)
PLE_ROOF       (:) = (1.-PFRAC_GR(:)) * PLE_STRLROOF(:) + PFRAC_GR(:) * PLE_GREENROOF(:) &
                   + (1 - PF_WASTE_CAN(:)) * PLE_WASTE(:)/PBLD(:)
!
!-------------------------------------------------------------------------------
!
!*      5.     Momentum fluxes
!              ---------------
!
PUSTAR_TOWN(:) = XUNDEF
!
IF (.NOT. OCANOPY) THEN
  !
  ZUSTAR2(:) = 0.0
  ZVMOD  (:) = PVMOD(:)
  !
  IF(HIMPLICIT_WIND=='OLD')THEN 
    ! old implicitation
    ZUSTAR2(:) = (PCD(:)*PVMOD(:)*PPEW_B_COEF(:))/    &
                 (1.0-PRHOA(:)*PCD(:)*PVMOD(:)*PPEW_A_COEF(:))
  ELSE
    ! new implicitation
    ZUSTAR2(:) = (PCD(:)*PVMOD(:)*(2.*PPEW_B_COEF(:)-PVMOD(:)))/ &
                 (1.0-2.0*PRHOA(:)*PCD(:)*PVMOD(:)*PPEW_A_COEF(:))
    !                   
    ZVMOD(:) = PRHOA(:)*PPEW_A_COEF(:)*ZUSTAR2(:) + PPEW_B_COEF(:)
    ZVMOD(:) = MAX(ZVMOD(:),0.)
    !
    WHERE(PPEW_A_COEF(:)/= 0.)
      ZUSTAR2(:) = MAX( ( ZVMOD(:) - PPEW_B_COEF(:) ) / (PRHOA(:)*PPEW_A_COEF(:)), 0.)
    ENDWHERE
    !               
  ENDIF
  !
  PUSTAR_TOWN(:) = SQRT(ZUSTAR2(:))     
  !  
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      6.     Averaged fluxes
!              ---------------
!
!*      6.1    Built fraction
!              --------------
!
PRN_BLT (:)    = (      PBLD(:)        * PRN_ROOF(:)        &
                  +     PROAD(:)       * PRN_ROAD(:)        &
                  + 0.5*PWALL_O_HOR(:) * PRN_WALL_A(:)      &
                  + 0.5*PWALL_O_HOR(:) * PRN_WALL_B(:))     &
                  / (PROAD(:) + PBLD(:))
!
PH_BLT  (:)    = (      PBLD(:)        * PH_ROOF(:)         &
                  +     PROAD(:)       * PH_ROAD(:)         &
                  + 0.5*PWALL_O_HOR(:) * PH_WALL_A(:)       &   
                  + 0.5*PWALL_O_HOR(:) * PH_WALL_B(:))      &   
                  / (PROAD(:) + PBLD(:))
!
PLE_BLT (:)    = (      PBLD(:)        * PLE_ROOF (:)       &
                  +     PROAD(:)       * PLE_ROAD (:)       &
                  + 0.5*PWALL_O_HOR(:) * PLE_WALL_A (:)     & 
                  + 0.5*PWALL_O_HOR(:) * PLE_WALL_B (:))    & 
                  / (PROAD(:) + PBLD(:))
!
PGFLUX_BLT (:) = (      PBLD(:)        * PGFLUX_ROOF (:)    &
                  +     PROAD(:)       * PGFLUX_ROAD (:)    &
                  + 0.5*PWALL_O_HOR(:) * PGFLUX_WALL_A (:)  &
                  + 0.5*PWALL_O_HOR(:) * PGFLUX_WALL_B (:)) &
                  / (PROAD(:) + PBLD(:))
!
PMELT_BLT  (:) = (     PBLD(:)         * PMELT_ROOF(:) * PDN_ROOF(:)  &
                  +    PROAD(:)        * PMELT_ROAD(:) * PDN_ROAD(:)) &
                  / (PROAD(:) + PBLD(:))
!
PDQS_BLT  (:) = (      PBLD(:)         * PDQS_ROOF (:)      &
                  +      PROAD(:)      * PDQS_ROAD (:)      &
                  + 0.5*PWALL_O_HOR(:) * PDQS_WALL_A (:)    &
                  + 0.5*PWALL_O_HOR(:) * PDQS_WALL_B (:) )  &
                / (PROAD(:) + PBLD(:))
!
IF (HBEM == "DEF") THEN
  PQF_BLD(:) = ( ZH_ROOF_SNOWFREE(:) + PLEW_ROOF(:) +                 & 
                 PDQS_ROOF(:) - ZRN_ROOF_SNOWFREE(:)  ) * PDF_ROOF(:) &
             + ( PDQS_ROOF(:) - PGSNOW_ROOF      (:)  ) * PDN_ROOF(:) &
             + 0.5*PWALL_O_HOR(:)/PBLD(:) * ( PH_WALL_A(:) + PLE_WALL_A(:) + PDQS_WALL_A(:) - PRN_WALL_A(:) ) &
             + 0.5*PWALL_O_HOR(:)/PBLD(:) * ( PH_WALL_B(:) + PLE_WALL_B(:) + PDQS_WALL_B(:) - PRN_WALL_B(:) )
  PFLX_BLD(:)= XUNDEF
ELSEIF (HBEM == "BEM") THEN
  PQF_BLD(:) = PQIN(:)*PN_FLOOR(:) + PHVAC_COOL(:) + PHVAC_HEAT(:)
  PFLX_BLD(:)=  PFLX_BLD_ROOF(:) + 0.5*PWALL_O_HOR(:)/PBLD(:)*PFLX_BLD_WALL_A(:) &
            +                    0.5*PWALL_O_HOR(:)/PBLD(:)*PFLX_BLD_WALL_B(:) &
            + PFLX_BLD_FLOOR(:) + PFLX_BLD_MASS(:)  
ENDIF
!
PQF_TOWN(:)= PBLD(:)*PQF_BLD(:) + PH_TRAFFIC(:) + PH_INDUSTRY(:) + PLE_TRAFFIC(:) + PLE_INDUSTRY(:)
!
!Flux from the building to its structure -> we need to add the component to the
!floor, the mass and the window
!
!
IF (LHOOK) CALL DR_HOOK('URBAN_FLUXES',1,ZHOOK_HANDLE)
!
END SUBROUTINE URBAN_FLUXES
