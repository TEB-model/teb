!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!   ##########################################################################
    SUBROUTINE ROAD_LAYER_E_BUDGET(T, B, PTSTEP, PDN_ROAD, PRHOA, PAC_ROAD, PAC_ROAD_WAT, &
                                   PLW_RAD, PPS, PQSAT_ROAD, PDELT_ROAD, PEXNS,           &
                                   PABS_SW_ROAD, PGSNOW_ROAD, PQ_LOWCAN, PT_LOWCAN,       &
                                   PTS_WALL_A, PTS_WALL_B, PTSNOW_ROAD, PTS_GARDEN,       &
                                   PLW_WA_TO_R, PLW_WB_TO_R, PLW_S_TO_R, PLW_WIN_TO_R,    &
                                   PEMIT_LW_ROAD, PDQS_ROAD, PABS_LW_ROAD, PHFREE_ROAD,   &
                                   PLEFREE_ROAD, PIMB_ROAD, PRR )
!   ##########################################################################
!
!!****  *ROAD_LAYER_E_BUDGET*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the evoultion of roads surface temperatures
!         
!     
!!**  METHOD
!     ------
!
!    6 : equations for evolution of Ts_road 
!        **********************************
!
!
!     dTr_1(t) / dt = 1/(dr_1*Cr_1) * (  Rn_r - H_r - LE_r 
!                                      - 2*Kr_1*(Tr_1-Tr_2)/(dr_1 +dr_2)       )
!
!     dTr_k(t) / dt = 1/(dr_k*Cr_k) * (- 2*Kr_k-1*(Tr_k-Tr_k-1)/(dr_k-1 +dr_k) 
!                                      - 2*Kr_k  *(Tr_k-Tr_k+1)/(dr_k+1 +dr_k) )
!
!       with
!
!   K*_k  = (d*_k+ d*_k+1)/(d*_k/k*_k+ d*_k+1/k*_k+1)
!
!   Rn_r = abs_Rg_r
!  - sigma * emis_r                                                   * Ts_r**4 (t+dt)
!  +         emis_r                       *    SVF_r                  * LWR
!  + sigma * emis_r * emis_w              * (1-SVF_r)                 * Ts_w**4 (t+dt)
!  +         emis_r            (1-emis_w) * (1-SVF_r)   *      SVF_w  * LWR
!  + sigma * emis_r * emis_w * (1-emis_w) * (1-SVF_r)   * (1-2*SVF_w) * Ts_w**4 (t+dt)
!  + sigma * emis_r * emis_r * (1-emis_w) * (1-SVF_r)   *      SVF_w  * Ts_r**4 (t+dt)
!
!  H_r  = rho Cp CH V ( Ts_r (t+dt) - Ta_canyon )
!
!  LE_r = rho Lv CH V ( qs_r (t+dt) - qa_canyon )
!
!
! The system is implicited (or semi-implicited).
!
! ZIMPL=1    ---> implicit system
! ZIMPL=0.5  ---> semi-implicit system
! ZIMPL=0    ---> explicit system
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
!!                  21/11/01 (V. Masson and A. Lemonsu) bug of latent flux
!!                           for very strong evaporation (all reservoir emptied
!!                           in one time-step)
!!                     02/11 (V. Masson) split of the routine for roads and walls separately
!!      G. Pigeon      09/2012: add heating/cooling of rain from air temperature
!!                             to surface road temp. for the road energy budget 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_BEM_n, ONLY : BEM_t
!
USE MODD_CSTS,ONLY : XCPD, XLVTT, XSTEFAN, XCL
!
USE MODE_THERMOS
!
USE MODI_LAYER_E_BUDGET
USE MODI_LAYER_E_BUDGET_GET_COEF
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
REAL,               INTENT(IN)    :: PTSTEP       ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROAD     ! road snow fraction
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA        ! rho
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROAD     ! aerodynamical conductance
!                                                 ! between road and canyon
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROAD_WAT ! aerodynamical conductance
!                                                 ! between road and canyon
!                                                 ! (for water)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD      ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PPS          ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PQSAT_ROAD   ! q_sat(Ts)
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_ROAD   ! fraction of water
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS        ! surface Exner function
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROAD ! absorbed solar radiation
REAL, DIMENSION(:), INTENT(IN)    :: PGSNOW_ROAD  ! road snow conduction
!                                                 ! heat fluxes at mantel
!                                                 ! base
REAL, DIMENSION(:), INTENT(IN)    :: PQ_LOWCAN    ! and specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PT_LOWCAN    ! low canyon air temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_WALL_A   ! wall surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_WALL_B   ! wall surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTSNOW_ROAD  ! road snow temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GARDEN   ! green area surface temperature
!
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_R   ! LW interactions wall  -> road 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WB_TO_R   ! LW interactions wall  -> road 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_R    ! LW interactions sky   -> road 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WIN_TO_R ! LW interactions window -> road 
!
REAL, DIMENSION(:), INTENT(OUT)   :: PEMIT_LW_ROAD! LW flux emitted by the road (W/m2 of road)
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_ROAD    !heat storage inside the road
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_ROAD ! absorbed infrared rad.
REAL, DIMENSION(:), INTENT(OUT)   :: PHFREE_ROAD  ! sensible heat flux on the
                                                  ! snow free part of the road [W m-2]
REAL, DIMENSION(:), INTENT(OUT)   :: PLEFREE_ROAD ! latent heat flux on the
                                                  ! snow free part of the road [W m-2]
REAL, DIMENSION(:), INTENT(OUT)   :: PIMB_ROAD    ! road residual energy imbalance 
                                                  ! for verification [W m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PRR          ! rain rate [kg m-2 s-1]

!
!*      0.2    declarations of local variables
!
REAL :: ZIMPL=1.0      ! implicit coefficient
REAL :: ZEXPL=0.0      ! explicit coefficient
!
REAL, DIMENSION(SIZE(T%XT_ROAD,1),SIZE(T%XT_ROAD,2)) :: ZA,& ! lower diag.
                                                    ZB,& ! main  diag.
                                                    ZC,& ! upper diag.
                                                    ZY   ! r.h.s.
!
REAL, DIMENSION(SIZE(PPS)) :: ZDN_ROAD    ! snow-covered surface fraction on road
REAL, DIMENSION(SIZE(PPS)) :: ZDF_ROAD    ! snow-free surface fraction on road
!
REAL, DIMENSION(SIZE(PPS)) :: ZDQSAT_ROAD ! dq_sat/dTs
REAL, DIMENSION(SIZE(PPS)) :: ZRHO_ACF_R  ! rho * conductance
!                                         !     * snow-free f.
REAL, DIMENSION(SIZE(PPS)) :: ZRHO_ACF_R_WAT ! rho * conductance for water
!                                         !     * snow-free f.
! thermal capacity times layer depth
REAL, DIMENSION(SIZE(PPS)) :: ZTS_ROAD    ! road surface temperature
REAL, DIMENSION(SIZE(PPS)) :: ZHEAT_RR    ! heat used too cool/heat the rain from the roof
REAL, DIMENSION(SIZE(PPS)) :: ZT_SKY      ! road surface temperature
!
INTEGER :: IROAD_LAYER           ! number of road layers
INTEGER :: JJ            ! loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ROAD_LAYER_E_BUDGET',0,ZHOOK_HANDLE)
!
 CALL LAYER_E_BUDGET_GET_COEF( T%XT_ROAD, PTSTEP, ZIMPL, T%XHC_ROAD, T%XTC_ROAD, T%XD_ROAD, &
                              ZA, ZB, ZC, ZY )
!
!*      1.     Layer thermal properties
!              ------------------------
!
IROAD_LAYER = SIZE(T%XT_ROAD,2)
!
DO JJ=1, SIZE(PDN_ROAD) 
  !
  ZDN_ROAD(JJ) = PDN_ROAD (JJ)
  ZDF_ROAD(JJ) = 1. - ZDN_ROAD (JJ)
  !
  !*      2.3    Surface temperatures
  !              --------------------
  !
  ZTS_ROAD(JJ) = T%XT_ROAD(JJ,1)
  !
  !*      2.2    flux properties
  !              ---------------
  !
  ZRHO_ACF_R    (JJ) = PRHOA(JJ) * PAC_ROAD(JJ)    * ZDF_ROAD(JJ)
  ZRHO_ACF_R_WAT(JJ) = PRHOA(JJ) * PAC_ROAD_WAT(JJ) * ZDF_ROAD(JJ)
  !
  !*     2.4   Sky temperature
  !            ---------------
  !
  ZT_SKY(JJ) = (PLW_RAD(JJ)/XSTEFAN)**0.25
  !
ENDDO
!
!*      2.4    qsat, dqsat/dTs, and humidity for roads
!              ---------------------------------------
! 
ZDQSAT_ROAD(:) = DQSAT(ZTS_ROAD(:),PPS(:),PQSAT_ROAD(:))
!
!-------------------------------------------------------------------------------
!
!*      3.     First road layers coefficients (in contact with outdoor env.)
!              -------------------------------------------------------------
!
DO JJ=1,SIZE(T%XT_ROAD,1)
  !
  ZB(JJ,1) = ZB(JJ,1) + ZIMPL * XCPD/PEXNS(JJ) * ZRHO_ACF_R(JJ) &
                      + ZIMPL * XLVTT * ZRHO_ACF_R_WAT(JJ) * PDELT_ROAD(JJ) * ZDQSAT_ROAD(JJ)
  !
  ZY(JJ,1) = ZY(JJ,1)  &
             + XCPD/PEXNS(JJ) * ZRHO_ACF_R(JJ) * ( PT_LOWCAN(JJ) - ZEXPL * ZTS_ROAD(JJ) ) &
             + ZDF_ROAD(JJ)*PABS_SW_ROAD(JJ) + ZDN_ROAD(JJ)*PGSNOW_ROAD(JJ)               &
             + XLVTT * ZRHO_ACF_R_WAT(JJ) * PDELT_ROAD(JJ)                                &
               * ( PQ_LOWCAN(JJ) - PQSAT_ROAD(JJ) + ZIMPL * ZDQSAT_ROAD(JJ) * ZTS_ROAD(JJ) )     
  !
  ZB(JJ,1) = ZB(JJ,1) &
             + ZIMPL * ZDF_ROAD(JJ) * ( PLW_S_TO_R(JJ) + PLW_WA_TO_R(JJ) + &
                                        PLW_WB_TO_R(JJ) + PLW_WIN_TO_R(JJ) + &
                                        PRR(JJ) * XCL ) ! heat/cool rain
  !
  ZY(JJ,1) = ZY(JJ,1) &
             + ZDF_ROAD(JJ) * (                                             &
               PLW_S_TO_R  (JJ) * (ZT_SKY    (JJ) - ZEXPL * ZTS_ROAD(JJ))   &
             + PLW_WIN_TO_R(JJ) * (B%XT_WIN1 (JJ) - ZEXPL * ZTS_ROAD(JJ))   &
             + PLW_WA_TO_R (JJ) * (PTS_WALL_A(JJ) - ZEXPL * ZTS_ROAD(JJ))   & 
             + PLW_WB_TO_R (JJ) * (PTS_WALL_B(JJ) - ZEXPL * ZTS_ROAD(JJ))   &
             + PRR(JJ) * XCL *    (PT_LOWCAN (JJ) - ZEXPL * ZTS_ROAD(JJ) ))   !heat/cool rain     
  !     
ENDDO
!
!
 CALL LAYER_E_BUDGET( T%XT_ROAD, PTSTEP, ZIMPL, T%XHC_ROAD, T%XTC_ROAD, T%XD_ROAD, &
                     ZA, ZB, ZC, ZY, PDQS_ROAD )
!
!-------------------------------------------------------------------------------
!
!*     12.    Road and wall absorbed infra-red radiation on snow-free surfaces
!             ----------------------------------------------------------------
!
!* absorbed LW
DO JJ=1,SIZE(T%XT_ROAD,1)
  !
  ! surface temperature used in energy balance
  ZTS_ROAD(JJ) = ZEXPL *  ZTS_ROAD(JJ) + ZIMPL * T%XT_ROAD(JJ,1)
  PABS_LW_ROAD(JJ) = PLW_S_TO_R  (JJ) * (ZT_SKY(JJ)     - ZTS_ROAD(JJ)) + &
                     PLW_WA_TO_R (JJ) * (PTS_WALL_A(JJ) - ZTS_ROAD(JJ)) + &
                     PLW_WB_TO_R (JJ) * (PTS_WALL_B(JJ) - ZTS_ROAD(JJ)) + &
                     PLW_WIN_TO_R(JJ) * (B%XT_WIN1(JJ)  - ZTS_ROAD(JJ))
  !
  !*     9.    Road emitted LW radiation on snow-free surfaces
  !            -----------------------------------------------
  PEMIT_LW_ROAD(JJ) = XSTEFAN * T%XT_ROAD(JJ,1)**4 + &
                      (1 - T%XEMIS_ROAD(JJ))/T%XEMIS_ROAD(JJ) * PABS_LW_ROAD(JJ)
  !
  !*      10.     road and wall sensible heat flux
  !              --------------------------------
  !
  PHFREE_ROAD(JJ) = ZRHO_ACF_R(JJ) * XCPD/PEXNS(JJ) * &
                   ( ZIMPL*T%XT_ROAD(JJ,1) + ZEXPL*ZTS_ROAD(JJ) - PT_LOWCAN(JJ) )
  !
  !*      11     road latent heat flux
  !              ---------------------
  !
  PLEFREE_ROAD(JJ) = ZRHO_ACF_R_WAT(JJ) * XLVTT * PDELT_ROAD(JJ) * &
                    ( PQSAT_ROAD(JJ) - PQ_LOWCAN(JJ) +             &
                     ZIMPL * ZDQSAT_ROAD(JJ) * (T%XT_ROAD(JJ,1) - ZTS_ROAD(JJ)) )
  ZHEAT_RR(JJ) = PRR(JJ) * XCL * (ZTS_ROAD(JJ) - PT_LOWCAN(JJ))
  !
  !*      12     heat storage inside roads
  !              -------------------------
  !
  !*      13     road energy residual imbalance for verification
  !              -----------------------------------------------
  !
  PIMB_ROAD(JJ) = PABS_SW_ROAD(JJ) + PABS_LW_ROAD(JJ) - PDQS_ROAD(JJ) &
               - ZDF_ROAD(JJ) * ( PHFREE_ROAD(JJ) + PLEFREE_ROAD(JJ)) &
               - ZDN_ROAD(JJ) *   PGSNOW_ROAD(JJ)
  !
ENDDO
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ROAD_LAYER_E_BUDGET',1,ZHOOK_HANDLE)
!
END SUBROUTINE ROAD_LAYER_E_BUDGET

