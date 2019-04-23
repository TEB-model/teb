!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!   ##########################################################################
    SUBROUTINE WALL_LAYER_E_BUDGET(TOP, T, B, PT_WL, PTS_WL_B, PTI_WL_B, PTSTEP, PDN_RD,     &
                                   PRHOA, PAC_WL, PAC_BLD, PLW_RAD, PPS, PEXNS, PABS_SW_WL,  &
                                   PT_CANYON, PTS_RD, PTSNOW_RD, PTS_GD, PTS_FL,             &
                                   PLW_WA_TO_WB, PLW_R_TO_W, PLW_G_TO_W, PLW_NR_TO_W,        &
                                   PLW_WIN_TO_W, PLW_S_TO_W, PFLX_BLD_WL, PDQS_WL,           &
                                   PABS_LW_WL, PEMIT_LW_WL, PH_WL, PIMB_WL, PRADHT_IN,       &
                                   PRAD_RF_WL, PRAD_WL_WIN, PRAD_WL_FL, PRAD_WL_MA,          &
                                   PCONV_WL_BLD, PLOAD_IN_WL ) 
!   ##########################################################################
!
!!****  *ROAD_WALL_LAYER_E_BUDGET*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the evoultion of roads and walls surface temperatures
!         
!     
!!**  METHOD
!     ------
!
!    6 : equations for evolution of Ts_road and Ts_wall simultaneously
!        *************************************************************
!
!     dTw_k(t) / dt = 1/(dw_k*Cw_k) * (- 2*Kw_k-1*(Tw_k-Tw_k-1)/(dw_k-1 +dw_k) 
!                                      - 2*Kw_k  *(Tw_k-Tw_k+1)/(dw_k+1 +dw_k) )
!
!     dTw_1(t) / dt = 1/(dw_1*Cw_1) * (  Rn_w - H_w - LE_w 
!                                      - 2*Kw_1*(Tw_1-Tw_2)/(dw_1 +dw_2)       )
!
!
!       with
!
!   K*_k  = (d*_k+ d*_k+1)/(d*_k/k*_k+ d*_k+1/k*_k+1)
!
!   Rn_w = abs_Rg_w 
!  - sigma * emis_w                                                   * Ts_w**4 (t+dt)
!  +         emis_w                       *      SVF_w                * LWR
!  + sigma * emis_w * emis_r              *      SVF_w                * Ts_r**4 (t+dt)
!  + sigma * emis_w * emis_w              * (1-2*SVF_w)               * Ts_w**4 (t+dt)
!  +         emis_w            (1-emis_r) *      SVF_r  *      SVF_w  * LWR
!  +         emis_w            (1-emis_w) *      SVF_w  * (1-2*SVF_w) * LWR
!  + sigma * emis_w * emis_w * (1-emis_r) *      SVF_w  * (1-  SVF_r) * Ts_w**4 (t+dt)
!  + sigma * emis_w * emis_w * (1-emis_w) * (1-2*SVF_w) * (1-2*SVF_w) * Ts_w**4 (t+dt)
!  + sigma * emis_w * emis_r * (1-emis_w) *      SVF_w  * (1-2*SVF_w) * Ts_r**4 (t+dt)
!
!  H_w  = rho Cp CH V ( Ts_w (t+dt) - Ta_canyon )
!
!  LE_w = rho Lv CH V ( qs_w (t+dt) - qa_canyon )
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
!!                     02/11 (V. Masson) splits the routine for road and walls separately
!!                     01/12 (V. Masson) separates the 2 walls
!!                     09/12 (G. Pigeon) modif internal convective coef convection
!!                     10/12 (G. Pigeon) add solar heat gain of indoor wall
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
USE MODD_CSTS,ONLY : XCPD, XSTEFAN
!
USE MODI_LAYER_E_BUDGET_GET_COEF
USE MODI_LAYER_E_BUDGET
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
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_WL    ! wall layers temperatures
REAL, DIMENSION(:),   INTENT(IN)  :: PTS_WL_B ! opposite wall surface temperature
REAL, DIMENSION(:),   INTENT(IN)  :: PTI_WL_B ! opposite wall internal temperature
REAL,               INTENT(IN)    :: PTSTEP       ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PDN_RD     ! snow-covered fraction on roads
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA        ! rho
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WL     ! aerodynamical conductance [m/s]
!                                                 ! between wall and canyon
REAL, DIMENSION(:), INTENT(IN)    :: PAC_BLD      ! aerodynamical conductance
                                                  ! inside the building itself
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD      ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PPS          ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS        ! surface Exner function
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WL ! absorbed solar radiation
REAL, DIMENSION(:), INTENT(IN)    :: PT_CANYON    ! air canyon temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_RD     ! road surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTSNOW_RD  ! road snow temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GD   ! green area surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_FL    ! floor layers temperatures [K]
!
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_WB ! LW interactions wall  -> opposite wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_W   ! LW interactions road -> wall ; DEF formulation 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_W   ! LW interactions garden -> wall ; DEF formulation
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_W   ! LW interactions sky   -> wall 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_W  ! LW interactions road(snow) -> wall 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WIN_TO_W ! Radiative heat trasfer coeff wall-window 
                                                  ! [W K-1 m-2]
!
REAL, DIMENSION(:), INTENT(OUT)   :: PFLX_BLD_WL! flux from bld to wall
REAL, DIMENSION(:), INTENT(INOUT) :: PDQS_WL    ! heat storage inside the wall 
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_WL ! absorbed infrared rad. [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PEMIT_LW_WL  ! LW flux emitted by the wall [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PH_WL      ! Sensible heat flux from wall to air [W/m2(wall)]
                                                  ! wall = facade - glazing
REAL, DIMENSION(:), INTENT(OUT)   :: PIMB_WL    ! wall residual energy imbalance 
                                                  ! for verification [W m-2]
REAL, DIMENSION(:),   INTENT(IN)  :: PRADHT_IN     ! Indoor radiant heat transfer coefficient
                                                    ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_RF_WL ! rad. fluxes from roof to wall [W m-2(roof)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WL_WIN  ! rad. fluxes from wall to win  [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WL_FL! rad. fluxes from wall to floor [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WL_MA ! rad. fluxes from wall to mass [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_WL_BLD ! conv. fluxes from wall to bld [W m-2(wall)]

REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_IN_WL  ! LOAD from solar heat gain + rad int. gains  W/m2 [Wall]


!
!*      0.2    declarations of local variables
!
!
REAL :: ZIMPL=1.0      ! implicit coefficient
REAL :: ZEXPL=0.0      ! explicit coefficient
!
REAL, DIMENSION(SIZE(PPS),SIZE(PT_WL,2)) ::  ZA,& ! lower diag.
                                               ZB,& ! main  diag.
                                               ZC,& ! upper diag.
                                               ZY   ! r.h.s.                       
!
REAL, DIMENSION(SIZE(PPS)) :: ZMTC_O_D_WL_IN
REAL, DIMENSION(SIZE(PPS)) :: ZDF_RD    ! Road snow free fraction
REAL, DIMENSION(SIZE(PPS)) :: ZRHO_ACF_W  ! rho * conductance
!                                         !     * snow-free f.
!
! thermal capacity times layer depth
REAL, DIMENSION(SIZE(PPS)) :: ZTS_WL       ! wall surface temperature
REAL, DIMENSION(SIZE(PPS)) :: ZTI_WL       ! wall indoor surface temperature
REAL, DIMENSION(SIZE(PPS)) :: ZTI_WL_CONV  ! wall indoor surface temperature for conv. flux
REAL, DIMENSION(SIZE(PPS)) :: ZT_SKY         ! sky temperature [K]
!
REAL, DIMENSION(SIZE(PPS)) :: ZTI_RF       ! Indoor roof temperature [K]
REAL, DIMENSION(SIZE(PPS)) :: ZDIF_RAD_WL_RF !diff between the rad flux that should receive the wall from the roof
                                                 ! and what it really receives [W m-2(bld)]
REAL, DIMENSION(SIZE(PPS)) :: ZRAD_WL_RF     ! rad flux between the wall and the roof computed for the wall balance
REAL, DIMENSION(SIZE(PPS)) :: ZF_WL_WL       ! View factor wall-wall inside the building
REAL, DIMENSION(SIZE(PPS)) :: ZCHTC_IN_WL  ! indoor convective heat transfer coeff wall [W m-2 K-1]
INTEGER :: IWL_LAYER           ! number of wall layers
INTEGER :: JJ                    ! loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WALL_LAYER_E_BUDGET',0,ZHOOK_HANDLE)
!
PFLX_BLD_WL  = XUNDEF
PIMB_WL      = XUNDEF
PRAD_WL_WIN  = XUNDEF
PRAD_WL_FL= XUNDEF
PRAD_WL_MA = XUNDEF
PCONV_WL_BLD = XUNDEF

!
 CALL LAYER_E_BUDGET_GET_COEF( PT_WL, PTSTEP, ZIMPL, T%XHC_WALL, T%XTC_WALL, T%XD_WALL, &
                              ZA, ZB, ZC, ZY )
!
IWL_LAYER = SIZE(PT_WL,2)
!
DO JJ=1,SIZE(PDN_RD)
  !
  ZDF_RD(JJ) = 1. - PDN_RD(JJ)
  !
  ZTS_WL(JJ) = PT_WL(JJ,1)
  ZTI_WL(JJ) = PT_WL(JJ, IWL_LAYER)
  !
  !*      2.1    outdoor convective flux properties 
  !              ----------------------------------
  !
  ZRHO_ACF_W (JJ) = PRHOA(JJ) * PAC_WL(JJ)
  !
  !*      2.2    Sky temperature
  !              ---------------
  !
  ZT_SKY(JJ) = (PLW_RAD(JJ)/XSTEFAN)**0.25
  !  
  !*      2.3    indoor average thermal conductivity
  !              -----------------------------------
  !
  IF (TOP%CBEM .EQ. "DEF") THEN
    ZMTC_O_D_WL_IN(JJ) = 2. * T%XTC_WALL(JJ,IWL_LAYER) / T%XD_WALL (JJ,IWL_LAYER)
    ZMTC_O_D_WL_IN(JJ) = 1./(  1./ZMTC_O_D_WL_IN(JJ)  + 1./(XCPD*PRHOA(JJ)*PAC_BLD(JJ)))
  ENDIF
ENDDO  
  !  
  !*      2.4    indoor convective coefficient
  !              -----------------------------
  !
  ZCHTC_IN_WL(:) = CHTC_VERT_DOE(PT_WL(:,IWL_LAYER), B%XTI_BLD(:))
  DO JJ=1,SIZE(ZCHTC_IN_WL)
     ZCHTC_IN_WL(JJ) = MAX(1., ZCHTC_IN_WL(JJ))
  ENDDO

!
!-------------------------------------------------------------------------------
!
!*      3.    Outer wall layer coefficients
!             ------------------------------
!
DO JJ=1,SIZE(PT_WL,1)
  !
  ZB(JJ,1) = ZB(JJ,1) + ZIMPL * XCPD/PEXNS(JJ) * ZRHO_ACF_W(JJ)
  !
  ZY(JJ,1) = ZY(JJ,1) + PABS_SW_WL(JJ)  &
                      + XCPD/PEXNS(JJ) * ZRHO_ACF_W(JJ) * ( PT_CANYON(JJ) - ZEXPL * ZTS_WL(JJ) )
  !
  !
  ZB(JJ,1) = ZB(JJ,1) &
             + ZIMPL * ( PLW_S_TO_W(JJ) + PLW_WA_TO_WB(JJ)                &
                       + ZDF_RD(JJ)*PLW_R_TO_W(JJ) +  PLW_G_TO_W(JJ)    &
                       + PDN_RD(JJ) *   PLW_NR_TO_W(JJ)                 &
                       + PLW_WIN_TO_W(JJ)  )    
  !
  ZY(JJ,1) = ZY(JJ,1) + &
                      PLW_S_TO_W(JJ) * (ZT_SKY(JJ)     - ZEXPL * ZTS_WL(JJ))  &
                    + PLW_WA_TO_WB(JJ) * (PTS_WL_B(JJ) - ZEXPL * ZTS_WL(JJ))  &
       + ZDF_RD(JJ) *  PLW_R_TO_W(JJ)  * (PTS_RD(JJ)   - ZEXPL * ZTS_WL(JJ))  &
       + PDN_RD(JJ) *  PLW_NR_TO_W(JJ) * (PTSNOW_RD(JJ)- ZEXPL * ZTS_WL(JJ))  &
                + PLW_WIN_TO_W(JJ) * (B%XT_WIN1(JJ)    - ZEXPL * ZTS_WL(JJ))
  !
  IF (SIZE(PTS_GD)>0) THEN
    ZY(JJ,1) = ZY(JJ,1) + PLW_G_TO_W(JJ)  * (PTS_GD(JJ)   - ZEXPL * ZTS_WL(JJ))
  ENDIF 
  !
ENDDO
!
!-------------------------------------------------------------------------------
!
!*      4.    Inside wall layer coefficients
!             -----------------------------
!
DO JJ=1,SIZE(PT_WL,1)
  !                
  IF (TOP%CBEM=="DEF") THEN
    !
    ZB(JJ,IWL_LAYER) = ZB(JJ,IWL_LAYER) + ZIMPL * ZMTC_O_D_WL_IN(JJ)
    !
    ZY(JJ,IWL_LAYER) = ZY(JJ,IWL_LAYER) &
                        + ZMTC_O_D_WL_IN(JJ) * B%XTI_BLD(JJ) &
                        - ZEXPL * ZMTC_O_D_WL_IN(JJ) * PT_WL(JJ,IWL_LAYER)
    !
  ELSEIF (TOP%CBEM=="BEM") THEN
    !
    ZF_WL_WL (JJ) = 1. - B%XF_WALL_MASS(JJ) - B%XF_WALL_WIN(JJ) - 2.*B%XF_WALL_FLOOR(JJ) 
    !
     ZB(JJ,IWL_LAYER) = ZB(JJ,IWL_LAYER) + ZIMPL *               &
                        (ZCHTC_IN_WL(JJ) * 4./3. + PRADHT_IN(JJ) *           &
                      (  B%XF_WALL_MASS(JJ) +     B%XF_WALL_WIN  (JJ) &
                          + ZF_WL_WL(JJ) + 2 * B%XF_WALL_FLOOR(JJ)))
    !
    ZTI_RF(JJ) = PRAD_RF_WL(JJ) / PRADHT_IN(JJ) + PT_WL(JJ,IWL_LAYER)
    !
    ZY(JJ,IWL_LAYER) = ZY(JJ,IWL_LAYER) +  &
        ZCHTC_IN_WL(JJ) * (B%XTI_BLD(JJ) - 1./3. * PT_WL(JJ, IWL_LAYER) * (4 * ZEXPL -1)) + &
        PRADHT_IN(JJ) * ( &
           B%XF_WALL_MASS (JJ) * (B%XT_MASS(JJ,1) - ZEXPL * PT_WL(JJ,IWL_LAYER)) + &
           B%XF_WALL_WIN  (JJ) * (B%XT_WIN2 (JJ)  - ZEXPL * PT_WL(JJ,IWL_LAYER)) + &
           B%XF_WALL_FLOOR(JJ) * (PTS_FL(JJ)      - ZEXPL * PT_WL(JJ,IWL_LAYER)) + &
           ZF_WL_WL (JJ) * (PTI_WL_B(JJ)          - ZEXPL * PT_WL(JJ,IWL_LAYER)) + &
           B%XF_WALL_FLOOR(JJ) * (ZTI_RF(JJ)      - ZEXPL * PT_WL(JJ,IWL_LAYER)) )+ &
           PLOAD_IN_WL(JJ)
    !
  ENDIF
  !
END DO
!
!-------------------------------------------------------------------------------
!
!*      5.    heat conduction calculation
!             ---------------------------
!
 CALL LAYER_E_BUDGET( PT_WL, PTSTEP, ZIMPL, T%XHC_WALL, T%XTC_WALL, T%XD_WALL, &
                     ZA, ZB, ZC, ZY, PDQS_WL )
!
!-------------------------------------------------------------------------------
!
!*   6.   diagnostics of flux echanged with the wall
!         ------------------------------------------
!
!
!* radiative surface temperature used during the energy balance
ZTS_WL(:) = ZIMPL * PT_WL(:,1) + ZEXPL * ZTS_WL(:)
!
PABS_LW_WL(:) = PLW_S_TO_W  (:) * (ZT_SKY   (:) - ZTS_WL(:)) + &
     ZDF_RD(:) *PLW_R_TO_W  (:) * (PTS_RD   (:) - ZTS_WL(:)) + &
                PLW_WA_TO_WB(:) * (PTS_WL_B (:) - ZTS_WL(:)) + &
                PLW_WIN_TO_W(:) * (B%XT_WIN1(:) - ZTS_WL(:)) + &
     PDN_RD(:) *PLW_NR_TO_W (:) * (PTSNOW_RD(:) - ZTS_WL(:))
!
IF (SIZE(PTS_GD)>0) THEN
  PABS_LW_WL(:) = PABS_LW_WL(:) + PLW_G_TO_W  (:) * (PTS_GD   (:) - ZTS_WL(:))
ENDIF
!
!* emitted lw flux
PEMIT_LW_WL(:) = XSTEFAN * PT_WL(:,1)**4 + &
                   (1 - T%XEMIS_WALL(:))/T%XEMIS_WALL(:) * PABS_LW_WL(:)
!
!* sensible heat flux to outdoor
PH_WL(:) = ZRHO_ACF_W(:) * XCPD/PEXNS(:) *  &
             ( ZIMPL*PT_WL(:,1) + ZEXPL*ZTS_WL(:) - PT_CANYON(:) )
!
IF (TOP%CBEM=='BEM') THEN
    !
    !compute ZTI_WALL used in flux calculation
    ZTI_WL_CONV(:) = 4./3. * ZIMPL * PT_WL(:,IWL_LAYER) + 1./3. * ZTI_WL(:) * (4 * ZEXPL -1.)
    ZTI_WL(:) = ZEXPL * ZTI_WL(:) + ZIMPL * PT_WL(:,IWL_LAYER) 
    !
    !compute IR exchanged fluxes with the roof in the wall balance
    ZRAD_WL_RF(:)  = PRADHT_IN(:)     * (ZTI_WL(:) - ZTI_RF(:))
    !compute the difference with the IR flux echanged between roof and wall for
    !the roof balance, this difference will be added to the floor 
    ZDIF_RAD_WL_RF(:)  = ZRAD_WL_RF(:) + PRAD_RF_WL(:)

    !compute exchanged fluxes with other surfaces for which the balance is done after
    PRAD_WL_FL(:) = PRADHT_IN(:)     * (ZTI_WL(:) - PTS_FL(:))
    ! contribution of the difference between the flux from the wall/flux from
    ! the roof
    PRAD_WL_FL(:)   = PRAD_WL_FL(:) + ZDIF_RAD_WL_RF(:) 
    PRAD_WL_MA(:)   = PRADHT_IN(:)     * (ZTI_WL(:) - B%XT_MASS(:,1))
    PRAD_WL_WIN(:)  = PRADHT_IN(:)     * (ZTI_WL(:) - B%XT_WIN2(:))
    PCONV_WL_BLD(:) = ZCHTC_IN_WL  (:) * (ZTI_WL_CONV(:) - B%XTI_BLD (:))
    !
    !
    PFLX_BLD_WL  (:) = -PCONV_WL_BLD(:) &
       + B%XF_WALL_FLOOR(:) * PRAD_RF_WL(:)                             &
       + PRADHT_IN(:) *(B%XF_WALL_MASS (:) * (B%XT_MASS(:,1) - ZTI_WL(:)) &
                      + B%XF_WALL_WIN  (:) * (B%XT_WIN2 (:)  - ZTI_WL(:)) &
                      + B%XF_WALL_FLOOR(:) * (PTS_FL(:)      - ZTI_WL(:)) )
!
!* wall energy residual imbalance for verification
!  -----------------------------------------------
  PIMB_WL(:) = PABS_SW_WL(:) + PABS_LW_WL(:)  &
               - PDQS_WL(:) - PH_WL(:)         &
               + PFLX_BLD_WL(:)
!
!
END IF
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WALL_LAYER_E_BUDGET',1,ZHOOK_HANDLE)
!
END SUBROUTINE WALL_LAYER_E_BUDGET
