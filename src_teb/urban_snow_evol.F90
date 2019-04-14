!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE URBAN_SNOW_EVOL(T, B, PT_LWCN, PQ_LWCN, PU_LWCN, PTS_RF, PTS_RD, PTS_WL_A,   &
                               PTS_WL_B, PPS, PTA, PQA, PRHOA, PLW_RAD, PSR, PZREF, PUREF,  &
                               PVMOD, PTSTEP, PZ_LWCN, PDN_RF, PABS_SW_SN_RF, PABS_LW_SN_RF,& 
                               PDN_RD, PABS_SW_SN_RD, PABS_LW_SN_RD, PRNSN_RF, PHSN_RF,     &
                               PLESN_RF, PGSN_RF, PMELT_RF, PRNSN_RD, PHSN_RD, PLESN_RD,    &
                               PGSN_RD, PMELT_RD, PLW_WA_TO_NR , PLW_WB_TO_NR, PLW_S_TO_NR, &
                               PLW_WIN_TO_NR, PDQS_SN_RF, PDQS_SN_RD   )  
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
!!      V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    23/01/98 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_BEM_n, ONLY : BEM_t
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
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(BEM_t), INTENT(INOUT) :: B
!
REAL, DIMENSION(:),   INTENT(IN)    :: PT_LWCN  ! LWCN air temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PQ_LWCN  ! LWCN air specific humidity
REAL, DIMENSION(:),   INTENT(IN)    :: PU_LWCN  ! LWCN hor. wind
REAL, DIMENSION(:),   INTENT(IN)    :: PTS_RF   ! roof surface temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PTS_RD   ! road surface temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PTS_WL_A ! wall surface temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PTS_WL_B ! wall surface temperature
!
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
REAL, DIMENSION(:), INTENT(IN)    :: PZ_LWCN  ! height of forcing
!
REAL, DIMENSION(:), INTENT(IN)    :: PDN_RF          ! snow-covered roof frac.
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_SN_RF ! SW absorbed by roof snow
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_SN_RF ! absorbed IR rad by snow on roof
REAL, DIMENSION(:), INTENT(INOUT) :: PDN_RD          ! snow-covered road frac.
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_SN_RD ! SW absorbed by road snow
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_SN_RD ! absorbed IR rad by snow on road
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSN_RF ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSN_RF  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESN_RF ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSN_RF  ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_RF   ! snow melt
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSN_RD ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSN_RD  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESN_RD ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSN_RD  ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_RD   ! snow melt
!
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_NR        ! LW contrib. wall       -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WB_TO_NR        ! LW contrib. wall       -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_NR         ! LW contrib. sky        -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WIN_TO_NR       ! LW contrib. win       -> road(snow)
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_SN_RF ! Heat storage in snowpack on roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_SN_RD ! Heat storage in snowpack on roads
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PTA)) :: ZLW1_RD   ! independant from
REAL, DIMENSION(SIZE(PTA)) :: ZLW1_RF   ! surface temperature
!
REAL, DIMENSION(SIZE(PTA)) :: ZLW2_RD   ! to be multiplied by
REAL, DIMENSION(SIZE(PTA)) :: ZLW2_RF   ! 4th power of
!                                         ! surface temperature

REAL, DIMENSION(SIZE(PTA)) :: ZSR_RF    ! snow fall on roof snow (kg/s/m2 of snow)
REAL, DIMENSION(SIZE(PTA)) :: ZSR_RD    ! snow fall on road snow (kg/s/m2 of snow)
!
REAL, DIMENSION(SIZE(PTA)) :: ZT_SKY      ! sky temperature
REAL, DIMENSION(SIZE(PTA)) :: ZTS_COEFA   ! Coefficient A for implicit coupling
!                                         ! of snow with the underlying surface
REAL, DIMENSION(SIZE(PTA)) :: ZTS_COEFB   ! Coefficient B for implicit coupling
!                                         ! of snow with the underlying surface
!
! flags to call to snow routines
!
LOGICAL :: GSN_RF, GSN_RD
!
! loop counters
!
INTEGER :: JL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('URBAN_SNOW_EVOL',0,ZHOOK_HANDLE)
PRNSN_RF(:)=0.
PHSN_RF (:)=0.
PLESN_RF(:)=0.
PGSN_RF (:)=0.
PMELT_RF(:)=0.
PRNSN_RD(:)=0.
PHSN_RD (:)=0.
PLESN_RD(:)=0.
PGSN_RD (:)=0.
PMELT_RD(:)=0.
PABS_LW_SN_RF(:)=0.
PABS_LW_SN_RD(:)=0.
!
!-------------------------------------------------------------------------------
!
GSN_RF = ANY( PSR(:)>0. .OR. T%TSNOW_ROOF%WSNOW(:,1)>0. )
GSN_RD = ANY( PSR(:)>0. .OR. T%TSNOW_ROAD%WSNOW(:,1)>0. )
!
!-------------------------------------------------------------------------------
!
!*      5.     Snow mantel model
!              -----------------
!
!*      5.1    roofs
!              -----
!
IF ( GSN_RF ) THEN
!
!* initializes LW radiative coefficients
!
  ZLW1_RF(:) =   T%TSNOW_ROOF%EMIS(:) * PLW_RAD(:)
  ZLW2_RF(:) = - T%TSNOW_ROOF%EMIS(:) * XSTEFAN
!
!* The global amount of snow on roofs is supposed located on a
!  fraction of the roof surface. All computations are then
!  done only for each m2 of snow, and not for each m2 of roof.
!
  DO JL=1,SIZE(T%TSNOW_ROOF%WSNOW,2)
    WHERE (PDN_RF(:)>0.) T%TSNOW_ROOF%WSNOW(:,JL) = T%TSNOW_ROOF%WSNOW(:,JL) / PDN_RF(:)
  END DO
  ZSR_RF=0.
  WHERE (PDN_RF(:)>0.) ZSR_RF   (:) = PSR   (:) / PDN_RF(:)
!
!* Estimates implicit coupling between snow and roof
! (strictly equal to an implicit formulation for 100% snow coverage)
!
  CALL ROOF_IMPL_COEF(T, PTSTEP, ZTS_COEFA, ZTS_COEFB)
!
!* call to snow mantel scheme
!
  IF (T%TSNOW_ROOF%SCHEME=='1-L')  &
   CALL SNOW_COVER_1LAYER(PTSTEP, XANSMIN_ROOF, XANSMAX_ROOF, XANS_TODRY_ROOF,   &
                          XRHOSMIN_ROOF, XRHOSMAX_ROOF, XANS_T_ROOF, .TRUE., 0., &
                          XWCRN_ROOF, XZ0SN, XZ0HSN, T%TSNOW_ROOF, PTS_RF,       &
                          ZTS_COEFA, ZTS_COEFB, PABS_SW_SN_RF, ZLW1_RF,  ZLW2_RF,&
                          PTA, PQA, PVMOD, PPS, PRHOA, ZSR_RF, PZREF, PUREF,     &
                          PRNSN_RF, PHSN_RF, PLESN_RF, PGSN_RF, PMELT_RF,        &
                          PDQS_SN_RF, PABS_LW_SN_RF)  
!

!
!* The global amount of snow on roofs is reported to total roof surface.
!
  DO JL=1,SIZE(T%TSNOW_ROOF%WSNOW,2)
    T%TSNOW_ROOF%WSNOW(:,JL) = T%TSNOW_ROOF%WSNOW(:,JL) * PDN_RF(:)
  END DO
!           
END IF
!
!*      5.2    roads
!              -----
!
IF ( GSN_RD ) THEN
  !
  ZT_SKY(:) = (PLW_RAD(:)/XSTEFAN)**0.25
!
  ZLW1_RD(:) = PLW_S_TO_NR  (:)   * (ZT_SKY   (:) - T%TSNOW_ROAD%TS(:)) &
               + PLW_WA_TO_NR (:) * (PTS_WL_A (:) - T%TSNOW_ROAD%TS(:)) &
               + PLW_WB_TO_NR (:) * (PTS_WL_B (:) - T%TSNOW_ROAD%TS(:)) &
               + PLW_WIN_TO_NR(:) * (B%XT_WIN1(:) - T%TSNOW_ROAD%TS(:))
  ZLW2_RD(:) =  0.0
  !
  !* The global amount of snow on roads is supposed located on a
  !  fraction of the road surface. All computations are then
  !  done only for each m2 of snow, and not for each m2 of road.
  !
  DO JL=1,SIZE(T%TSNOW_ROAD%WSNOW,2)
    WHERE (PDN_RD(:)>0.) T%TSNOW_ROAD%WSNOW(:,JL) = T%TSNOW_ROAD%WSNOW(:,JL) / PDN_RD(:)
  END DO
  ZSR_RD=0.
  WHERE (PDN_RD(:)>0.) ZSR_RD   (:) = PSR   (:) / PDN_RD(:)
  !
  !* no implicit coupling necessary with road
  !
  ZTS_COEFA = 0.
  ZTS_COEFB = PTS_RD
  !
  !* call to snow mantel scheme
  !
  IF (T%TSNOW_ROAD%SCHEME=='1-L')                &
    CALL SNOW_COVER_1LAYER(PTSTEP, XANSMIN_ROAD, XANSMAX_ROAD, XANS_TODRY_ROAD,   &
                           XRHOSMIN_ROAD, XRHOSMAX_ROAD, XANS_T_ROAD, .FALSE.,    &
                           0., XWCRN_ROAD, XZ0SN, XZ0HSN, T%TSNOW_ROAD, PTS_RD,   & 
                           ZTS_COEFA, ZTS_COEFB,  PABS_SW_SN_RD, ZLW1_RD, ZLW2_RD,&
                           PT_LWCN, PQ_LWCN, PU_LWCN, PPS, PRHOA, ZSR_RD, PZ_LWCN,&
                           PZ_LWCN, PRNSN_RD, PHSN_RD, PLESN_RD, PGSN_RD,         &
                           PMELT_RD, PDQS_SN_RD, PABS_LW_SN_RD        )  
!
!* The global amount of snow on roads is reported to total road surface.
!
  DO JL=1,SIZE(T%TSNOW_ROAD%WSNOW,2)
    T%TSNOW_ROAD%WSNOW(:,JL) = T%TSNOW_ROAD%WSNOW(:,JL) * PDN_RD(:)
  END DO
!
  WHERE (T%TSNOW_ROAD%T(:,1) .EQ. XUNDEF) PDN_RD(:) = 0.0
!
END IF
IF (LHOOK) CALL DR_HOOK('URBAN_SNOW_EVOL',1,ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE URBAN_SNOW_EVOL
