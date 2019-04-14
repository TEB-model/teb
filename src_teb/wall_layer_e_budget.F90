!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ##########################################################################
    SUBROUTINE WALL_LAYER_E_BUDGET(HBEM, PT_WALL, PTS_WALL_B, PTI_WALL_B, PTSTEP,          &
                                   PHC_WALL, PTC_WALL, PD_WALL, PDN_ROAD, PRHOA, PAC_WALL, &
                                   PAC_BLD, PTI_BLD, PLW_RAD, PPS, PEXNS,     &
                                   PABS_SW_WALL, PT_CANYON, PTS_ROAD, PTSNOW_ROAD,         &
                                   PTS_GARDEN, PTS_MASS, PTS_FLOOR, PEMIS_WALL,            &
                                   PLW_WA_TO_WB, PLW_R_TO_W, PLW_G_TO_W, PLW_NR_TO_W,      &
                                   PLW_WIN_TO_W, PLW_S_TO_W,                               &
                                   PT_WIN1,                                                &
                                   PFLX_BLD_WALL, PDQS_WALL, PT_WIN2, PABS_LW_WALL,        &
                                   PEMIT_LW_WALL, PH_WALL, PIMB_WALL, PF_WALL_MASS,        &
                                   PF_WALL_FLOOR, PF_WALL_WIN, PRADHT_IN, PRAD_ROOF_WALL,  &
                                   PRAD_WALL_WIN, PRAD_WALL_FLOOR,                         &
                                   PRAD_WALL_MASS, PCONV_WALL_BLD, PLOAD_IN_WALL ) 
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
!!	V. Masson           * Meteo-France *
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
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_WALL    ! wall layers temperatures
REAL, DIMENSION(:),   INTENT(IN)  :: PTS_WALL_B ! opposite wall surface temperature
REAL, DIMENSION(:),   INTENT(IN)  :: PTI_WALL_B ! opposite wall internal temperature
REAL,               INTENT(IN)    :: PTSTEP       ! time step
 CHARACTER(LEN=3), INTENT(IN)      :: HBEM         ! Building Energy model 'DEF' or 'BEM'
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_WALL     ! heat capacity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_WALL     ! thermal conductivity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_WALL      ! depth of wall layers
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROAD     ! snow-covered fraction on roads
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA        ! rho
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WALL     ! aerodynamical conductance [m/s]
!                                                 ! between wall and canyon
REAL, DIMENSION(:), INTENT(IN)    :: PAC_BLD      ! aerodynamical conductance
                                                  ! inside the building itself
REAL, DIMENSION(:), INTENT(IN)    :: PTI_BLD      ! inside building temperature
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD      ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PPS          ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS        ! surface Exner function
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WALL ! absorbed solar radiation
REAL, DIMENSION(:), INTENT(IN)    :: PT_CANYON    ! air canyon temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_ROAD     ! road surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTSNOW_ROAD  ! road snow temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GARDEN   ! green area surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_MASS     ! surface mass temperature  [K]
REAL, DIMENSION(:), INTENT(IN)    :: PTS_FLOOR    ! floor layers temperatures [K]
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_WALL   ! wall emissivity
!
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_WB ! LW interactions wall  -> opposite wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_W   ! LW interactions road -> wall ; DEF formulation 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_W   ! LW interactions garden -> wall ; DEF formulation
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_W   ! LW interactions sky   -> wall 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_W  ! LW interactions road(snow) -> wall 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WIN_TO_W ! Radiative heat trasfer coeff wall-window 
                                                  ! [W K-1 m-2]
!
REAL, DIMENSION(:), INTENT(IN)    :: PT_WIN1      ! outdoor window temperature [K]
REAL, DIMENSION(:), INTENT(OUT)   :: PFLX_BLD_WALL! flux from bld to wall
REAL, DIMENSION(:), INTENT(INOUT) :: PDQS_WALL    ! heat storage inside the wall 
REAL, DIMENSION(:), INTENT(IN)    :: PT_WIN2      ! indoor window temperature [K]
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_WALL ! absorbed infrared rad. [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PEMIT_LW_WALL  ! LW flux emitted by the wall [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PH_WALL      ! Sensible heat flux from wall to air [W/m²(wall)]
                                                  ! wall = facade - glazing
REAL, DIMENSION(:), INTENT(OUT)   :: PIMB_WALL    ! wall residual energy imbalance 
                                                  ! for verification [W m-2]
REAL, DIMENSION(:),   INTENT(IN)  :: PF_WALL_FLOOR ! View factor wall-floor
REAL, DIMENSION(:),   INTENT(IN)  :: PF_WALL_MASS  ! View factor wall-mass
REAL, DIMENSION(:),   INTENT(IN)  :: PF_WALL_WIN   ! View factor wall-win
REAL, DIMENSION(:),   INTENT(IN)  :: PRADHT_IN     ! Indoor radiant heat transfer coefficient
                                                    ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_ROOF_WALL ! rad. fluxes from roof to wall [W m-2(roof)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WALL_WIN  ! rad. fluxes from wall to win  [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WALL_FLOOR! rad. fluxes from wall to floor [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WALL_MASS ! rad. fluxes from wall to mass [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_WALL_BLD ! conv. fluxes from wall to bld [W m-2(wall)]

REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_IN_WALL  ! LOAD from solar heat gain + rad int. gains  W/m² [Wall]


!
!*      0.2    declarations of local variables
!
!
REAL :: ZIMPL=1.0      ! implicit coefficient
REAL :: ZEXPL=0.0      ! explicit coefficient
!
REAL, DIMENSION(SIZE(PPS),SIZE(PT_WALL,2)) ::  ZA,& ! lower diag.
                                               ZB,& ! main  diag.
                                               ZC,& ! upper diag.
                                               ZY   ! r.h.s.                       
!
REAL, DIMENSION(SIZE(PPS)) :: ZMTC_O_D_WALL_IN
REAL, DIMENSION(SIZE(PPS)) :: ZDF_ROAD    ! Road snow free fraction
REAL, DIMENSION(SIZE(PPS)) :: ZRHO_ACF_W  ! rho * conductance
!                                         !     * snow-free f.
!
! thermal capacity times layer depth
REAL, DIMENSION(SIZE(PPS)) :: ZTS_WALL       ! wall surface temperature
REAL, DIMENSION(SIZE(PPS)) :: ZTI_WALL       ! wall indoor surface temperature
REAL, DIMENSION(SIZE(PPS)) :: ZTI_WALL_CONV  ! wall indoor surface temperature for conv. flux
REAL, DIMENSION(SIZE(PPS)) :: ZT_SKY         ! sky temperature [K]
!
REAL, DIMENSION(SIZE(PPS)) :: ZTI_ROOF       ! Indoor roof temperature [K]
REAL, DIMENSION(SIZE(PPS)) :: ZDIF_RAD_WALL_ROOF !diff between the rad flux that should receive the wall from the roof
                                                 ! and what it really receives [W m-2(bld)]
REAL, DIMENSION(SIZE(PPS)) :: ZRAD_WALL_ROOF     ! rad flux between the wall and the roof computed for the wall balance
REAL, DIMENSION(SIZE(PPS)) :: ZF_WALL_WALL       ! View factor wall-wall inside the building
REAL, DIMENSION(SIZE(PPS)) :: ZCHTC_IN_WALL  ! indoor convective heat transfer coeff wall [W m-2 K-1]
INTEGER :: IWALL_LAYER           ! number of wall layers
INTEGER :: JJ                    ! loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WALL_LAYER_E_BUDGET',0,ZHOOK_HANDLE)
!
PFLX_BLD_WALL  = XUNDEF
PIMB_WALL      = XUNDEF
PRAD_WALL_WIN  = XUNDEF
PRAD_WALL_FLOOR= XUNDEF
PRAD_WALL_MASS = XUNDEF
PCONV_WALL_BLD = XUNDEF

!
 CALL LAYER_E_BUDGET_GET_COEF( PT_WALL, PTSTEP, ZIMPL, PHC_WALL, PTC_WALL, PD_WALL, &
                              ZA, ZB, ZC, ZY )
!
IWALL_LAYER = SIZE(PT_WALL,2)
!
DO JJ=1,SIZE(PDN_ROAD)
  !
  ZDF_ROAD(JJ) = 1. - PDN_ROAD(JJ)
  !
  ZTS_WALL(JJ) = PT_WALL(JJ,1)
  ZTI_WALL(JJ) = PT_WALL(JJ, IWALL_LAYER)
  !
  !*      2.1    outdoor convective flux properties 
  !              ----------------------------------
  !
  ZRHO_ACF_W (JJ) = PRHOA(JJ) * PAC_WALL(JJ)
  !
  !*      2.2    Sky temperature
  !              ---------------
  !
  ZT_SKY(JJ) = (PLW_RAD(JJ)/XSTEFAN)**0.25
  !  
  !*      2.3    indoor average thermal conductivity
  !              -----------------------------------
  !
  IF (HBEM .EQ. "DEF") THEN
    ZMTC_O_D_WALL_IN(JJ) = 2. * PTC_WALL(JJ,IWALL_LAYER) / PD_WALL (JJ,IWALL_LAYER)
    ZMTC_O_D_WALL_IN(JJ) = 1./(  1./ZMTC_O_D_WALL_IN(JJ)  + 1./(XCPD*PRHOA(JJ)*PAC_BLD(JJ)))
  ENDIF
ENDDO  
  !  
  !*      2.4    indoor convective coefficient
  !              -----------------------------
  !
  ZCHTC_IN_WALL(:) = CHTC_VERT_DOE(PT_WALL(:,IWALL_LAYER), PTI_BLD(:))
  DO JJ=1,SIZE(ZCHTC_IN_WALL)
     ZCHTC_IN_WALL(JJ) = MAX(1., ZCHTC_IN_WALL(JJ))
  ENDDO

!
!-------------------------------------------------------------------------------
!
!*      3.    Outer wall layer coefficients
!             ------------------------------
!
DO JJ=1,SIZE(PT_WALL,1)
  !
  ZB(JJ,1) = ZB(JJ,1) + ZIMPL * XCPD/PEXNS(JJ) * ZRHO_ACF_W(JJ)
  !
  ZY(JJ,1) = ZY(JJ,1) + PABS_SW_WALL(JJ)  &
                      + XCPD/PEXNS(JJ) * ZRHO_ACF_W(JJ) * ( PT_CANYON(JJ) - ZEXPL * ZTS_WALL(JJ) )
  !
  !
  ZB(JJ,1) = ZB(JJ,1) &
             + ZIMPL * ( PLW_S_TO_W(JJ) + PLW_WA_TO_WB(JJ)                &
                       + ZDF_ROAD(JJ)*PLW_R_TO_W(JJ) +  PLW_G_TO_W(JJ)    &
                       + PDN_ROAD(JJ) *   PLW_NR_TO_W(JJ)                 &
                       + PLW_WIN_TO_W(JJ)  )    
  !
  ZY(JJ,1) = ZY(JJ,1) + &
                         PLW_S_TO_W(JJ) * (ZT_SKY(JJ)     - ZEXPL * ZTS_WALL(JJ))  &
                       + PLW_WA_TO_WB(JJ) * (PTS_WALL_B(JJ) - ZEXPL * ZTS_WALL(JJ))  &
       + ZDF_ROAD(JJ) *  PLW_R_TO_W(JJ)  * (PTS_ROAD(JJ)   - ZEXPL * ZTS_WALL(JJ))  &
                       + PLW_G_TO_W(JJ)  * (PTS_GARDEN(JJ) - ZEXPL * ZTS_WALL(JJ)) &
       + PDN_ROAD(JJ) *  PLW_NR_TO_W(JJ) * (PTSNOW_ROAD(JJ)- ZEXPL * ZTS_WALL(JJ))  &
                       + PLW_WIN_TO_W(JJ) * (PT_WIN1(JJ)    - ZEXPL * ZTS_WALL(JJ))  

  !
ENDDO
!
!-------------------------------------------------------------------------------
!
!*      4.    Inside wall layer coefficients
!             -----------------------------
!
DO JJ=1,SIZE(PT_WALL,1)
  !                
  IF (HBEM=="DEF") THEN
    !
    ZB(JJ,IWALL_LAYER) = ZB(JJ,IWALL_LAYER) + ZIMPL * ZMTC_O_D_WALL_IN(JJ)
    !
    ZY(JJ,IWALL_LAYER) = ZY(JJ,IWALL_LAYER) &
                        + ZMTC_O_D_WALL_IN(JJ) * PTI_BLD(JJ) &
                        - ZEXPL * ZMTC_O_D_WALL_IN(JJ) * PT_WALL(JJ,IWALL_LAYER)
    !
  ELSEIF (HBEM=="BEM") THEN
    !
    ZF_WALL_WALL (JJ) = 1. - PF_WALL_MASS(JJ) - PF_WALL_WIN(JJ) - 2.*PF_WALL_FLOOR(JJ) 
    !
     ZB(JJ,IWALL_LAYER) = ZB(JJ,IWALL_LAYER) + ZIMPL *               &
                        (ZCHTC_IN_WALL(JJ) * 4./3. + PRADHT_IN(JJ) *           &
                         (  PF_WALL_MASS(JJ) +     PF_WALL_WIN  (JJ) &
                          + ZF_WALL_WALL(JJ) + 2 * PF_WALL_FLOOR(JJ)))
    !
    ZTI_ROOF(JJ) = PRAD_ROOF_WALL(JJ) / PRADHT_IN(JJ) + PT_WALL(JJ,IWALL_LAYER)
    ZY(JJ,IWALL_LAYER) = ZY(JJ,IWALL_LAYER) +  &
        ZCHTC_IN_WALL(JJ) * (PTI_BLD(JJ) - 1./3. * PT_WALL(JJ, IWALL_LAYER) * (4 * ZEXPL -1)) + &
        PRADHT_IN(JJ) * ( &
           PF_WALL_MASS (JJ) * (PTS_MASS(JJ)  - ZEXPL * PT_WALL(JJ,IWALL_LAYER)) + &
           PF_WALL_WIN  (JJ) * (PT_WIN2 (JJ)  - ZEXPL * PT_WALL(JJ,IWALL_LAYER)) + &
           PF_WALL_FLOOR(JJ) * (PTS_FLOOR(JJ) - ZEXPL * PT_WALL(JJ,IWALL_LAYER)) + &
           ZF_WALL_WALL (JJ) * (PTI_WALL_B(JJ)- ZEXPL * PT_WALL(JJ,IWALL_LAYER)) + &
           PF_WALL_FLOOR(JJ) * (ZTI_ROOF(JJ) - ZEXPL * PT_WALL(JJ,IWALL_LAYER)) )+ &
           PLOAD_IN_WALL(JJ)
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
 CALL LAYER_E_BUDGET( PT_WALL, PTSTEP, ZIMPL, PHC_WALL, PTC_WALL, PD_WALL, &
                     ZA, ZB, ZC, ZY, PDQS_WALL )
!
!-------------------------------------------------------------------------------
!
!*   6.   diagnostics of flux echanged with the wall
!         ------------------------------------------
!
!
!* radiative surface temperature used during the energy balance
ZTS_WALL(:) = ZIMPL * PT_WALL(:,1) + ZEXPL * ZTS_WALL(:)
!
PABS_LW_WALL(:) = PLW_S_TO_W  (:) * (ZT_SKY     (:) - ZTS_WALL(:)) + &
     ZDF_ROAD(:) *PLW_R_TO_W  (:) * (PTS_ROAD   (:) - ZTS_WALL(:)) + &
                  PLW_G_TO_W  (:) * (PTS_GARDEN (:) - ZTS_WALL(:)) + &
                  PLW_WA_TO_WB(:) * (PTS_WALL_B (:) - ZTS_WALL(:)) + &
                  PLW_WIN_TO_W(:) * (PT_WIN1    (:) - ZTS_WALL(:)) + &
     PDN_ROAD(:) *PLW_NR_TO_W (:) * (PTSNOW_ROAD(:) - ZTS_WALL(:))
!
!* emitted lw flux
PEMIT_LW_WALL(:) = XSTEFAN * PT_WALL(:,1)**4 + &
                   (1 - PEMIS_WALL(:))/PEMIS_WALL(:) * PABS_LW_WALL(:)
!
!* sensible heat flux to outdoor
PH_WALL(:) = ZRHO_ACF_W(:) * XCPD/PEXNS(:) *  &
             ( ZIMPL*PT_WALL(:,1) + ZEXPL*ZTS_WALL(:) - PT_CANYON(:) )
!
IF (HBEM=='BEM') THEN
    !
    !compute ZTI_WALL used in flux calculation
    ZTI_WALL_CONV(:) = 4./3. * ZIMPL * PT_WALL(:,IWALL_LAYER) + 1./3. * ZTI_WALL(:) * (4 * ZEXPL -1.)
    ZTI_WALL(:) = ZEXPL * ZTI_WALL(:) + ZIMPL * PT_WALL(:,IWALL_LAYER) 
    !
    !compute IR exchanged fluxes with the roof in the wall balance
    ZRAD_WALL_ROOF(:)  = PRADHT_IN(:)     * (ZTI_WALL(:) - ZTI_ROOF(:))
    !compute the difference with the IR flux echanged between roof and wall for
    !the roof balance, this difference will be added to the floor 
    ZDIF_RAD_WALL_ROOF(:)  = ZRAD_WALL_ROOF(:) + PRAD_ROOF_WALL(:)

    !compute exchanged fluxes with other surfaces for which the balance is done after
    PRAD_WALL_FLOOR(:) = PRADHT_IN(:)     * (ZTI_WALL(:) - PTS_FLOOR(:))
    ! contribution of the difference between the flux from the wall/flux from
    ! the roof
    PRAD_WALL_FLOOR(:) = PRAD_WALL_FLOOR(:) + ZDIF_RAD_WALL_ROOF(:) 
    PRAD_WALL_MASS(:)  = PRADHT_IN(:)     * (ZTI_WALL(:) - PTS_MASS(:))
    PRAD_WALL_WIN(:)   = PRADHT_IN(:)     * (ZTI_WALL(:) - PT_WIN2(:))
    PCONV_WALL_BLD(:)  = ZCHTC_IN_WALL  (:) * (ZTI_WALL_CONV(:) - PTI_BLD (:))
    !
    !
    PFLX_BLD_WALL  (:) = -PCONV_WALL_BLD(:) &
       + PF_WALL_FLOOR(:) * PRAD_ROOF_WALL(:)                             &
       + PRADHT_IN(:) *(PF_WALL_MASS (:) * (PTS_MASS (:) - ZTI_WALL(:)) &
                       + PF_WALL_WIN  (:) * (PT_WIN2 (:)  - ZTI_WALL(:)) &
                       + PF_WALL_FLOOR(:) * (PTS_FLOOR(:) - ZTI_WALL(:)) )
!
!* wall energy residual imbalance for verification
!  -----------------------------------------------
  PIMB_WALL(:) = PABS_SW_WALL(:) + PABS_LW_WALL(:)  &
               - PDQS_WALL(:) - PH_WALL(:)         &
               + PFLX_BLD_WALL(:)
!
!
END IF
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WALL_LAYER_E_BUDGET',1,ZHOOK_HANDLE)
!
END SUBROUTINE WALL_LAYER_E_BUDGET
