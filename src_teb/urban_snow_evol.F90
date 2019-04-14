!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     #########
    SUBROUTINE URBAN_SNOW_EVOL(                                                 &
                       PT_LOWCAN, PQ_LOWCAN, PU_LOWCAN,                         &
                       PTS_ROOF,PTS_ROAD,PTS_WALL_A, PTS_WALL_B,                &
                       PT_ROOF, PD_ROOF, PTC_ROOF, PHC_ROOF,                    &
                       HSNOW_ROOF,                                              &
                       PWSNOW_ROOF, PTSNOW_ROOF, PRSNOW_ROOF, PASNOW_ROOF,      &
                       PTSSNOW_ROOF, PESNOW_ROOF,                               &
                       HSNOW_ROAD,                                              &
                       PWSNOW_ROAD, PTSNOW_ROAD, PRSNOW_ROAD, PASNOW_ROAD,      &
                       PTSSNOW_ROAD, PESNOW_ROAD,                               &
                       PPS, PTA, PQA, PRHOA,                                    &
                       PLW_RAD,                                                 &
                       PSR, PZREF, PUREF, PVMOD,                                &
                       PTSTEP,                                                  &
                       PZ_LOWCAN,                                               &
                       PDN_ROOF, PABS_SW_SNOW_ROOF, PABS_LW_SNOW_ROOF,          &
                       PDN_ROAD, PABS_SW_SNOW_ROAD, PABS_LW_SNOW_ROAD,          &
                       PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,    &
                       PMELT_ROOF,                                              &
                       PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,    &
                       PMELT_ROAD,                                              &
                       PLW_WA_TO_NR , PLW_WB_TO_NR, PLW_S_TO_NR, PLW_WIN_TO_NR, &
                       PDQS_SNOW_ROOF, PDQS_SNOW_ROAD, PT_WIN1                  )  
!   ##########################################################################
!
!!****  *URBAN_SNOW_EVOL*  
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
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SNOW_PAR, ONLY : XZ0SN, XZ0HSN,                                    &
                            XANSMIN_ROOF, XANSMAX_ROOF, XANS_TODRY_ROOF,      &
                            XANS_T_ROOF, XRHOSMIN_ROOF, XRHOSMAX_ROOF,        &
                            XWCRN_ROOF,                                       &
                            XANSMIN_ROAD, XANSMAX_ROAD, XANS_TODRY_ROAD,      &
                            XANS_T_ROAD, XRHOSMIN_ROAD, XRHOSMAX_ROAD,        &
                            XWCRN_ROAD  
USE MODD_CSTS,     ONLY : XSTEFAN
!
USE MODE_SURF_SNOW_FRAC
!
USE MODI_ROOF_IMPL_COEF
USE MODI_SNOW_COVER_1LAYER
!
USE MODD_SURF_PAR, ONLY : XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
REAL, DIMENSION(:),   INTENT(IN)    :: PT_LOWCAN  ! LOWCAN air temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PQ_LOWCAN  ! LOWCAN air specific humidity
REAL, DIMENSION(:),   INTENT(IN)    :: PU_LOWCAN  ! LOWCAN hor. wind
REAL, DIMENSION(:),   INTENT(IN)    :: PTS_ROOF   ! roof surface temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PTS_ROAD   ! road surface temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PTS_WALL_A ! wall surface temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PTS_WALL_B ! wall surface temperature
REAL, DIMENSION(:,:), INTENT(IN)    :: PT_ROOF    ! roof temperature profile
REAL, DIMENSION(:,:), INTENT(IN)    :: PD_ROOF    ! roof layer thickness
REAL, DIMENSION(:,:), INTENT(IN)    :: PTC_ROOF   ! roof layer thermal conductivity
REAL, DIMENSION(:,:), INTENT(IN)    :: PHC_ROOF   ! roof layer heat capacity
 CHARACTER(LEN=*),     INTENT(IN)    :: HSNOW_ROOF ! snow roof scheme
!                                                 ! 'NONE'
!                                                 ! 'D95 '
!                                                 ! '1-L '
 CHARACTER(LEN=*),     INTENT(IN)    :: HSNOW_ROAD ! snow road scheme
!                                                 ! 'NONE'
!                                                 ! 'D95 '
!                                                 ! '1-L '
REAL, DIMENSION(:,:), INTENT(INOUT) :: PWSNOW_ROOF ! snow layers reservoir
REAL, DIMENSION(:,:), INTENT(INOUT) :: PTSNOW_ROOF ! snow layers temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRSNOW_ROOF ! snow layers density
REAL, DIMENSION(:),   INTENT(INOUT) :: PASNOW_ROOF ! snow albedo
REAL, DIMENSION(:),   INTENT(INOUT) :: PESNOW_ROOF ! snow emissivity
REAL, DIMENSION(:),   INTENT(INOUT) :: PTSSNOW_ROOF! snow surface temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PWSNOW_ROAD ! snow layers reservoir
REAL, DIMENSION(:,:), INTENT(INOUT) :: PTSNOW_ROAD ! snow layers temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRSNOW_ROAD ! snow layers density
REAL, DIMENSION(:),   INTENT(INOUT) :: PASNOW_ROAD ! snow albedo
REAL, DIMENSION(:),   INTENT(INOUT) :: PESNOW_ROAD ! snow emissivity
REAL, DIMENSION(:),   INTENT(INOUT) :: PTSSNOW_ROAD! snow surface temperature

REAL, DIMENSION(:), INTENT(IN)    :: PPS      ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PTA      ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA      ! specific humidity
                                              ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD    ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA    ! air density at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD  ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PSR      ! snow rate
REAL, DIMENSION(:), INTENT(IN)    :: PZREF    ! reference height of the first
                                              ! atmospheric level (temperature)
REAL, DIMENSION(:), INTENT(IN)    :: PUREF    ! reference height of the first
                                              ! atmospheric level (wind)
                                              ! at first atmospheric level
REAL,               INTENT(IN)    :: PTSTEP   ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PZ_LOWCAN  ! height of forcing
!
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROOF          ! snow-covered roof frac.
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_SNOW_ROOF ! SW absorbed by roof snow
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_SNOW_ROOF ! absorbed IR rad by snow on roof
REAL, DIMENSION(:), INTENT(INOUT) :: PDN_ROAD          ! snow-covered road frac.
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_SNOW_ROAD ! SW absorbed by road snow
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_SNOW_ROAD ! absorbed IR rad by snow on road
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSNOW_ROOF ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSNOW_ROOF  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESNOW_ROOF ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSNOW_ROOF  ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_ROOF   ! snow melt
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSNOW_ROAD ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSNOW_ROAD  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESNOW_ROAD ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSNOW_ROAD  ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_ROAD   ! snow melt
!
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_NR        ! LW contrib. wall       -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WB_TO_NR        ! LW contrib. wall       -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_NR         ! LW contrib. sky        -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WIN_TO_NR       ! LW contrib. win       -> road(snow)
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_SNOW_ROOF ! Heat storage in snowpack on roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_SNOW_ROAD ! Heat storage in snowpack on roads
REAL, DIMENSION(:), INTENT(IN)    :: PT_WIN1        ! Window surface temperature
!
!*      0.2    declarations of local variables
!
!
REAL, DIMENSION(SIZE(PTA)) :: ZLW1_ROAD   ! independant from
REAL, DIMENSION(SIZE(PTA)) :: ZLW1_ROOF   ! surface temperature
!
REAL, DIMENSION(SIZE(PTA)) :: ZLW2_ROAD   ! to be multiplied by
REAL, DIMENSION(SIZE(PTA)) :: ZLW2_ROOF   ! 4th power of
!                                         ! surface temperature

REAL, DIMENSION(SIZE(PTA)) :: ZSR_ROOF    ! snow fall on roof snow (kg/s/m2 of snow)
REAL, DIMENSION(SIZE(PTA)) :: ZSR_ROAD    ! snow fall on road snow (kg/s/m2 of snow)
!
REAL, DIMENSION(SIZE(PTA)) :: ZT_SKY      ! sky temperature
REAL, DIMENSION(SIZE(PTA)) :: ZTS_COEFA   ! Coefficient A for implicit coupling
!                                         ! of snow with the underlying surface
REAL, DIMENSION(SIZE(PTA)) :: ZTS_COEFB   ! Coefficient B for implicit coupling
!                                         ! of snow with the underlying surface
!
! flags to call to snow routines
!
LOGICAL :: GSNOW_ROOF, GSNOW_ROAD
!
! loop counters
!
INTEGER :: JL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('URBAN_SNOW_EVOL',0,ZHOOK_HANDLE)
PRNSNOW_ROOF(:)=0.
PHSNOW_ROOF (:)=0.
PLESNOW_ROOF(:)=0.
PGSNOW_ROOF (:)=0.
PMELT_ROOF  (:)=0.
PRNSNOW_ROAD(:)=0.
PHSNOW_ROAD (:)=0.
PLESNOW_ROAD(:)=0.
PGSNOW_ROAD (:)=0.
PMELT_ROAD  (:)=0.
PABS_LW_SNOW_ROOF(:)=0.
PABS_LW_SNOW_ROAD(:)=0.
!
!-------------------------------------------------------------------------------
!
GSNOW_ROOF = ANY( PSR(:)>0. .OR. PWSNOW_ROOF(:,1)>0. )
GSNOW_ROAD = ANY( PSR(:)>0. .OR. PWSNOW_ROAD(:,1)>0. )
!
!-------------------------------------------------------------------------------
!
!*      5.     Snow mantel model
!              -----------------
!
!*      5.1    roofs
!              -----
!
IF ( GSNOW_ROOF ) THEN
!
!* initializes LW radiative coefficients
!
  ZLW1_ROOF(:) =   PESNOW_ROOF(:) * PLW_RAD(:)
  ZLW2_ROOF(:) = - PESNOW_ROOF(:) * XSTEFAN
!
!* The global amount of snow on roofs is supposed located on a
!  fraction of the roof surface. All computations are then
!  done only for each m2 of snow, and not for each m2 of roof.
!
  DO JL=1,SIZE(PWSNOW_ROOF,2)
    WHERE (PDN_ROOF(:)>0.) PWSNOW_ROOF(:,JL) = PWSNOW_ROOF(:,JL) / PDN_ROOF(:)
  END DO
  ZSR_ROOF=0.
  WHERE (PDN_ROOF(:)>0.) ZSR_ROOF   (:) = PSR   (:) / PDN_ROOF(:)
!
!* Estimates implicit coupling between snow and roof
! (strictly equal to an implicit formulation for 100% snow coverage)
!
  CALL ROOF_IMPL_COEF(PTSTEP, SIZE(PT_ROOF,2), PD_ROOF, PTC_ROOF, PHC_ROOF, PT_ROOF, ZTS_COEFA,ZTS_COEFB)
!
!* call to snow mantel scheme
!
  IF (HSNOW_ROOF=='1-L')  &
   CALL SNOW_COVER_1LAYER(PTSTEP, XANSMIN_ROOF, XANSMAX_ROOF, XANS_TODRY_ROOF, &
                           XRHOSMIN_ROOF, XRHOSMAX_ROOF, XANS_T_ROOF, .TRUE.,   &
                           0., XWCRN_ROOF,                                      &
                           XZ0SN,XZ0HSN,                                        &
                           PTSNOW_ROOF(:,1), PASNOW_ROOF,                       &
                           PRSNOW_ROOF(:,1), PWSNOW_ROOF(:,1), PTSSNOW_ROOF,    &
                           PESNOW_ROOF,                                         &
                           PTS_ROOF,  ZTS_COEFA, ZTS_COEFB, PABS_SW_SNOW_ROOF,  &
                           ZLW1_ROOF, ZLW2_ROOF,                                &
                           PTA, PQA, PVMOD, PPS, PRHOA, ZSR_ROOF, PZREF, PUREF, &
                           PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,&
                           PMELT_ROOF, PDQS_SNOW_ROOF, PABS_LW_SNOW_ROOF        )  
!

!
!* The global amount of snow on roofs is reported to total roof surface.
!
  DO JL=1,SIZE(PWSNOW_ROOF,2)
    PWSNOW_ROOF(:,JL) = PWSNOW_ROOF(:,JL) * PDN_ROOF(:)
  END DO
!           
END IF
!
!*      5.2    roads
!              -----
!
IF ( GSNOW_ROAD ) THEN
  !
  ZT_SKY(:) = (PLW_RAD(:)/XSTEFAN)**0.25
!
  ZLW1_ROAD(:) = PLW_S_TO_NR  (:) * (ZT_SKY    (:) - PTSSNOW_ROAD(:)) &
               + PLW_WA_TO_NR (:) * (PTS_WALL_A(:) - PTSSNOW_ROAD(:)) &
               + PLW_WB_TO_NR (:) * (PTS_WALL_B(:) - PTSSNOW_ROAD(:)) &
               + PLW_WIN_TO_NR(:) * (PT_WIN1   (:) - PTSSNOW_ROAD(:))
  ZLW2_ROAD(:) =  0.0
  !
  !* The global amount of snow on roads is supposed located on a
  !  fraction of the road surface. All computations are then
  !  done only for each m2 of snow, and not for each m2 of road.
  !
  DO JL=1,SIZE(PWSNOW_ROAD,2)
    WHERE (PDN_ROAD(:)>0.) PWSNOW_ROAD(:,JL) = PWSNOW_ROAD(:,JL) / PDN_ROAD(:)
  END DO
  ZSR_ROAD=0.
  WHERE (PDN_ROAD(:)>0.) ZSR_ROAD   (:) = PSR   (:) / PDN_ROAD(:)
  !
  !* no implicit coupling necessary with road
  !
  ZTS_COEFA = 0.
  ZTS_COEFB = PTS_ROAD
  !
  !* call to snow mantel scheme
  !
  IF (HSNOW_ROAD=='1-L')                                                        &
    CALL SNOW_COVER_1LAYER(PTSTEP, XANSMIN_ROAD, XANSMAX_ROAD, XANS_TODRY_ROAD, &
                           XRHOSMIN_ROAD, XRHOSMAX_ROAD, XANS_T_ROAD, .FALSE.,  &
                           0., XWCRN_ROAD,                                      &
                           XZ0SN,XZ0HSN,                                        &
                           PTSNOW_ROAD(:,1), PASNOW_ROAD,                       &
                           PRSNOW_ROAD(:,1), PWSNOW_ROAD(:,1), PTSSNOW_ROAD,    &
                           PESNOW_ROAD,                                         &
                           PTS_ROAD, ZTS_COEFA, ZTS_COEFB,                      &
                           PABS_SW_SNOW_ROAD, ZLW1_ROAD, ZLW2_ROAD,             &
                           PT_LOWCAN, PQ_LOWCAN, PU_LOWCAN, PPS, PRHOA,         &
                           ZSR_ROAD, PZ_LOWCAN, PZ_LOWCAN,                      &
                           PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,&
                           PMELT_ROAD, PDQS_SNOW_ROAD ,PABS_LW_SNOW_ROAD        )  
!
!* The global amount of snow on roads is reported to total road surface.
!
  DO JL=1,SIZE(PWSNOW_ROAD,2)
    PWSNOW_ROAD(:,JL) = PWSNOW_ROAD(:,JL) * PDN_ROAD(:)
  END DO
!
  WHERE (PTSNOW_ROAD(:,1) .EQ. XUNDEF) PDN_ROAD(:) = 0.0
!
END IF
IF (LHOOK) CALL DR_HOOK('URBAN_SNOW_EVOL',1,ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE URBAN_SNOW_EVOL
