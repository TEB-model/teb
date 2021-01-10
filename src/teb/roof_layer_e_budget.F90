!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE ROOF_LAYER_E_BUDGET(TOP, T, B, PQSAT_ROOF, PAC_BLD, PTSTEP, PDN_ROOF,   &
                                   PRHOA, PAC_ROOF, PAC_ROOF_WAT, PLW_RAD, PPS,        &
                                   PDELT_ROOF, PTA, PQA, PEXNA, PEXNS, PABS_SW_ROOF,   &
                                   PGSNOW_ROOF, PFLX_BLD_ROOF, PDQS_ROOF, PABS_LW_ROOF,&
                                   PHFREE_ROOF, PLEFREE_ROOF, PIMB_ROOF,               &
                                   PG_GREENROOF_ROOF, PRADHT_IN, PTS_FLOOR, PTI_WALL,  &
                                   PRAD_ROOF_WALL, PRAD_ROOF_WIN, PRAD_ROOF_FLOOR,     &
                                   PRAD_ROOF_MASS, PCONV_ROOF_BLD, PRR, PLOAD_IN_ROOF )
!   ##########################################################################
!
!!****  *ROOF_LAYER_E_BUDGET*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the evoultion of surface temperature of roofs
!         
!     
!!**  METHOD
!     ------
!
!
!
!
!    5 : equation for evolution of Ts_roof
!        *********************************
!
!     dTt_1(t) / dt = 1/(dt_1*Ct_1) * (  Rn - H - LE
!                                      - 2*Kt_1*(Tt_1-Tt_2)/(dt_1 +dt_2)       )
!
!     dTt_k(t) / dt = 1/(dt_k*Ct_k) * (- 2*Kt_k-1*(Tt_k-Tt_k-1)/(dt_k-1 +dt_k) 
!                                      - 2*Kt_k  *(Tt_k-Tt_k+1)/(dt_k+1 +dt_k) )
!
!       with
!
!       K*_k  = (d*_k+ d*_k+1)/(d*_k/k*_k+ d*_k+1/k*_k+1)
!
!       Rn = (dir_Rg + sca_Rg) (1-a) + emis * ( Rat - sigma Ts**4 (t+dt) )
!
!       H  = rho Cp CH V ( Ts (t+dt) - Tas )
!
!       LE = rho Lv CH V ( qs (t+dt) - qas )
!
!      where the as subscript denotes atmospheric values at ground level
!      (and not at first half level)
!
!      The tridiagonal systel is linearized with
!
!       using      Ts**4 (t+dt) = Ts**4 (t) + 4*Ts**3 (t) * ( Ts(t+dt) - Ts(t) )
!
!       and  qs (t+dt) = Hu(t) * qsat(t) + Hu(t) dqsat/dT * ( Ts(t+dt) - Ts(t) )
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
!!                  17/10/05 (G. Pigeon) computation of storage inside the roofs
!!                  26/04/12 (G. Pigeon) add term of heating of rain (new arg PRR+XCL)
!!                     09/12 (G. Pigeon) modif of indoor conv. coef and implicit calculation
!!                     10/12 (G. Pigeon) add indoor solar heat load
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_BEM_n, ONLY : BEM_t
!
USE MODD_SURF_PAR,  ONLY : XUNDEF
USE MODD_CSTS,ONLY : XCPD, XLVTT, XSTEFAN, XCL
!
USE MODE_THERMOS
!
USE MODI_LAYER_E_BUDGET
USE MODI_LAYER_E_BUDGET_GET_COEF
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
REAL, DIMENSION(:), INTENT(INOUT) :: PQSAT_ROOF     ! q_sat(Ts)
REAL, DIMENSION(:), INTENT(IN)    :: PAC_BLD        ! aerodynamical resistance
                                                    ! inside building itself
REAL,               INTENT(IN)    :: PTSTEP         ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROOF       ! roof snow fraction
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA          ! air density
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROOF       ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROOF_WAT   ! aerodynamical conductance (for water)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD        ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PPS            ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_ROOF     ! fraction of water
REAL, DIMENSION(:), INTENT(IN)    :: PTA            ! air temperature at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PQA            ! air specific humidity
                                                    ! at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PEXNA          ! exner function
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS          ! surface exner function
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROOF   ! absorbed solar radiation
REAL, DIMENSION(:), INTENT(IN)    :: PGSNOW_ROOF    ! roof snow conduction
!                                                   ! heat fluxes at mantel
!                                                   ! base
REAL, DIMENSION(:), INTENT(IN)    :: PG_GREENROOF_ROOF ! heat conduction flux
!                                                        between greenroof and
!                                                        structural roof
REAL, DIMENSION(:), INTENT(OUT)   :: PFLX_BLD_ROOF  ! flux from bld to roof
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_ROOF      ! heat storage inside the roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_ROOF   ! absorbed infra-red rad.
REAL, DIMENSION(:), INTENT(OUT)   :: PHFREE_ROOF    ! sensible heat flux of the
                                                    ! snow free part of the roof
REAL, DIMENSION(:), INTENT(OUT)   :: PLEFREE_ROOF   ! latent heat flux of the
                                                    ! snow free part of the roof
REAL, DIMENSION(:), INTENT(OUT)   :: PIMB_ROOF      ! residual energy imbalance
                                                    ! of the roof for
REAL, DIMENSION(:), INTENT(IN)    :: PRADHT_IN      ! Indoor radiant heat transfer coefficient
                                                    ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PTS_FLOOR      ! surf. floor temp. (contact with bld air)
REAL, DIMENSION(:), INTENT(IN)    :: PTI_WALL       ! indoor wall temp.
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_ROOF_WALL ! rad. fluxes from roof to wall [W m-2(roof)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_ROOF_WIN  ! rad. fluxes from roof to win [W m-2(roof)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_ROOF_FLOOR! rad. fluxes from roof to floor [W m-2(roof)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_ROOF_MASS ! rad. fluxes from roof to mass [W m-2(roof)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_ROOF_BLD ! conv. fluxes from roof to bld [W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PRR ! rain rate [kg m-2 s-1]
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_IN_ROOF ! solar + int heat gain on roof W/m2 [roof]
!
!*      0.2    declarations of local variables
!
REAL :: ZIMPL = 1.0        ! implicit coefficient
REAL :: ZEXPL = 0.0        ! explicit coefficient
!
REAL, DIMENSION(SIZE(PTA)) :: ZDF_ROOF ! snow-free fraction
REAL, DIMENSION(SIZE(PTA),SIZE(T%XT_ROOF,2)) :: ZA,& ! lower diag.
                                              ZB,& ! main  diag.
                                              ZC,& ! upper diag.
                                              ZY   ! r.h.s.
!
REAL, DIMENSION(SIZE(PTA)) :: ZDQSAT_ROOF      ! dq_sat/dTs
REAL, DIMENSION(SIZE(PTA)) :: ZRHO_ACF_ROOF    ! conductance * rho
REAL, DIMENSION(SIZE(PTA)) :: ZRHO_ACF_ROOF_WAT! conductance * rho (for water)
REAL, DIMENSION(SIZE(PTA)) :: ZMTC_O_D_ROOF_IN ! thermal capacity times layer depth
REAL, DIMENSION(SIZE(PTA)) :: ZTS_ROOF         ! roof surface temperature at previous time step
REAL, DIMENSION(SIZE(PTA)) :: ZTRAD_ROOF       ! roof radiative surface temperature at intermediate time step
REAL, DIMENSION(SIZE(PTA)) :: ZTAER_ROOF       ! roof aerodyn. surface temperature at intermediate time step
REAL, DIMENSION(SIZE(PTA)) :: ZHEAT_RR         ! heat used too cool/heat the rain from the roof
REAL, DIMENSION(SIZE(PTA)) :: ZTI_ROOF         ! temperature of internal roof layer used for radiative exchanges
REAL, DIMENSION(SIZE(PTA)) :: ZTI_ROOF_CONV    ! temperature of internal roof layer used for convective exchanges
REAL, DIMENSION(SIZE(PTA)) :: ZCHTC_IN_ROOF      ! Indoor roof convec heat transfer coefficient
                                               ! [W K-1 m-2(bld)]
!
INTEGER :: JJ
INTEGER :: IROOF_LAYER           ! number of roof layers
INTEGER :: JLAYER                ! loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ROOF_LAYER_E_BUDGET',0,ZHOOK_HANDLE)
!
PRAD_ROOF_WALL(:) = XUNDEF
PRAD_ROOF_WIN(:)  = XUNDEF
PRAD_ROOF_FLOOR(:)= XUNDEF
PRAD_ROOF_MASS(:) = XUNDEF
PCONV_ROOF_BLD(:) = XUNDEF
!
! *Convection heat transfer coefficients [W m-2 K-1] from EP Engineering Reference
!
IROOF_LAYER = SIZE(T%XT_ROOF,2)
!
ZCHTC_IN_ROOF(:) = CHTC_DOWN_DOE(T%XT_ROOF(:,IROOF_LAYER), B%XTI_BLD(:))
DO JJ=1,SIZE(ZCHTC_IN_ROOF)
   ZCHTC_IN_ROOF(JJ) = MAX(1., ZCHTC_IN_ROOF(JJ))
ENDDO
!
 CALL LAYER_E_BUDGET_GET_COEF(T%XT_ROOF, PTSTEP, ZIMPL, T%XHC_ROOF, T%XTC_ROOF, T%XD_ROOF, &
                              ZA, ZB, ZC, ZY )
!
!
DO JJ=1,SIZE(PDN_ROOF)
  !
  ZDF_ROOF(JJ) = 1. - PDN_ROOF(JJ)
  !
  ZTS_ROOF(JJ) = T%XT_ROOF(JJ,1)
  ZTI_ROOF(JJ) = T%XT_ROOF(JJ,IROOF_LAYER)
  !
  !*      2.     Roof Ts coefficients
  !              --------------------
  !
  ZRHO_ACF_ROOF    (JJ) = PRHOA(JJ) * PAC_ROOF    (JJ)
  ZRHO_ACF_ROOF_WAT(JJ) = PRHOA(JJ) * PAC_ROOF_WAT(JJ)
  !
  IF (TOP%CBEM .EQ. 'DEF') THEN
    ZMTC_O_D_ROOF_IN(JJ) = 2. * T%XTC_ROOF(JJ,IROOF_LAYER) / T%XD_ROOF (JJ,IROOF_LAYER)
    ZMTC_O_D_ROOF_IN(JJ) = 1./(  1./ZMTC_O_D_ROOF_IN(JJ) + 1./(XCPD*PRHOA(JJ)*PAC_BLD(JJ)) ) 
  ENDIF
  !
ENDDO
!
!*      2.1    dqsat/dTs, and humidity for roofs
!              ---------------------------------
!
ZDQSAT_ROOF(:) = DQSAT(ZTS_ROOF(:),PPS(:),PQSAT_ROOF(:))
!
!*      2.2    coefficients
!              ------------
! 
DO JJ=1,SIZE(T%XT_ROOF,1)
  !
  ZB(JJ,1) = ZB(JJ,1) + ZDF_ROOF(JJ) * (1.-T%XGREENROOF(JJ)) * (                           &
                        ZIMPL * ( XCPD/PEXNS(JJ) * ZRHO_ACF_ROOF(JJ)                       &
                        + XLVTT * ZRHO_ACF_ROOF_WAT(JJ) * PDELT_ROOF(JJ) * ZDQSAT_ROOF(JJ) &
                        + XSTEFAN * T%XEMIS_ROOF(JJ) * 4.*ZTS_ROOF(JJ)**3                  &
                        + PRR(JJ) * XCL )) !! heating/cooling of rain 
  !
  ZY(JJ,1) = ZY(JJ,1) + (1.-T%XGREENROOF(JJ))                                                            &
                      * (PDN_ROOF(JJ)*PGSNOW_ROOF(JJ) + ZDF_ROOF(JJ) * ( PABS_SW_ROOF(JJ)                &
                        + XCPD * ZRHO_ACF_ROOF(JJ) * ( PTA(JJ)/PEXNA(JJ) - ZEXPL*ZTS_ROOF(JJ)/PEXNS(JJ)) &
                        + T%XEMIS_ROOF(JJ)*PLW_RAD(JJ)                                                   &                 
                        + XLVTT * ZRHO_ACF_ROOF_WAT(JJ) * PDELT_ROOF(JJ)                                 &
                          * ( PQA(JJ) - PQSAT_ROOF(JJ) + ZIMPL * ZDQSAT_ROOF(JJ) * ZTS_ROOF(JJ) )        &
                        + XSTEFAN * T%XEMIS_ROOF(JJ) * ZTS_ROOF(JJ)**4 * ( 3.*ZIMPL-ZEXPL )              &
                        + PRR(JJ) * XCL * (PTA(JJ) - ZEXPL * ZTS_ROOF(JJ)) ) ) & !! heating/cooling of rain
                        + T%XGREENROOF(JJ)*PG_GREENROOF_ROOF(JJ)
  !
  IF (TOP%CBEM=="DEF") THEN
    !
    ZB(JJ,IROOF_LAYER) = ZB(JJ,IROOF_LAYER) + ZIMPL * ZMTC_O_D_ROOF_IN(JJ)
    !
    ZY(JJ,IROOF_LAYER) = ZY(JJ,IROOF_LAYER) &
                         + ZMTC_O_D_ROOF_IN(JJ) * B%XTI_BLD(JJ) &
                         - ZEXPL * ZMTC_O_D_ROOF_IN(JJ) * T%XT_ROOF(JJ,IROOF_LAYER)
    !
  ELSEIF (TOP%CBEM=="BEM") THEN
    !
    ZB(JJ, IROOF_LAYER) = ZB(JJ,IROOF_LAYER) + ZIMPL * &
                         (ZCHTC_IN_ROOF(JJ) * 4./3. + PRADHT_IN(JJ) * &
                         (B%XF_FLOOR_MASS(JJ) + B%XF_FLOOR_WIN(JJ) + &
                          B%XF_FLOOR_WALL(JJ) + B%XF_FLOOR_ROOF(JJ) ))

    ZY(JJ,IROOF_LAYER) = ZY(JJ,IROOF_LAYER) + &
       ZCHTC_IN_ROOF(JJ) * (B%XTI_BLD(JJ) - 1./3. *  T%XT_ROOF(JJ, IROOF_LAYER)*(4*ZEXPL - 1.)) + &
       PRADHT_IN(JJ) * ( &
          B%XF_FLOOR_MASS (JJ) * (B%XT_MASS(JJ,1) - ZEXPL * T%XT_ROOF(JJ,IROOF_LAYER)) + &
          B%XF_FLOOR_WIN  (JJ) * (B%XT_WIN2  (JJ) - ZEXPL * T%XT_ROOF(JJ,IROOF_LAYER)) + &
          B%XF_FLOOR_WALL (JJ) * (PTI_WALL   (JJ) - ZEXPL * T%XT_ROOF(JJ,IROOF_LAYER)) + &
          B%XF_FLOOR_ROOF (JJ) * (PTS_FLOOR  (JJ) - ZEXPL * T%XT_ROOF(JJ,IROOF_LAYER))) + &
          PLOAD_IN_ROOF(JJ)
    !
  ENDIF
  !
ENDDO
!
!
 CALL LAYER_E_BUDGET( T%XT_ROOF, PTSTEP, ZIMPL, T%XHC_ROOF, T%XTC_ROOF, T%XD_ROOF, &
                     ZA, ZB, ZC, ZY, PDQS_ROOF )
!
!-------------------------------------------------------------------------------
!
!*     diagnostic: computation of flux between bld and internal roof layer
DO JJ=1,SIZE(T%XT_ROOF,1)
  !
  ZTI_ROOF_CONV(JJ) = 4./3. * ZIMPL * T%XT_ROOF(JJ, IROOF_LAYER) + 1./3. * ZTI_ROOF(JJ) * (4*ZEXPL -1.)
  ZTI_ROOF(JJ) = ZEXPL * ZTI_ROOF(JJ) + ZIMPL * T%XT_ROOF(JJ, IROOF_LAYER) 
  SELECT CASE(TOP%CBEM)
  CASE("DEF")
     PFLX_BLD_ROOF(JJ) = ZMTC_O_D_ROOF_IN(JJ) * (B%XTI_BLD(JJ) - ZTI_ROOF(JJ))
  CASE("BEM")
     PRAD_ROOF_WALL(JJ) = PRADHT_IN(JJ)     * (ZTI_ROOF(JJ) - PTI_WALL(JJ))
     PRAD_ROOF_WIN(JJ)  = PRADHT_IN(JJ)     * (ZTI_ROOF(JJ) - B%XT_WIN2(JJ))
     PRAD_ROOF_FLOOR(JJ)= PRADHT_IN(JJ)     * (ZTI_ROOF(JJ) - PTS_FLOOR(JJ))
     PRAD_ROOF_MASS(JJ) = PRADHT_IN(JJ)     * (ZTI_ROOF(JJ) - B%XT_MASS(JJ,1))
     PCONV_ROOF_BLD(JJ) = ZCHTC_IN_ROOF  (JJ) * (ZTI_ROOF_CONV(JJ) - B%XTI_BLD (JJ))
     PFLX_BLD_ROOF(JJ) = -(PRAD_ROOF_WALL(JJ) + PRAD_ROOF_WIN(JJ) + PRAD_ROOF_FLOOR(JJ) + &
                            PRAD_ROOF_MASS(JJ) + PCONV_ROOF_BLD(JJ))
  ENDSELECT
  
  !
  !*      8.     Infra-red radiation absorbed by roofs
  !              -------------------------------------
  !
  !* radiative surface temperature at intermediate time step
  ZTRAD_ROOF(JJ) = ( ZTS_ROOF(JJ)**4 + &
                   4.*ZIMPL*ZTS_ROOF(JJ)**3 * (T%XT_ROOF(JJ,1) - ZTS_ROOF(JJ)) )**0.25
  !
  !* absorbed LW
  PABS_LW_ROOF(JJ) = T%XEMIS_ROOF(JJ) * (PLW_RAD(JJ) - XSTEFAN * ZTRAD_ROOF(JJ)** 4)
  !
  !*      9.     Sensible heat flux between snow free roof and air
  !              -------------------------------------------------
  !
  !* aerodynamic surface temperature at the intermediate time step
  ZTAER_ROOF(JJ) = ZEXPL * ZTS_ROOF(JJ) + ZIMPL * T%XT_ROOF(JJ,1)
  PHFREE_ROOF(JJ) = ZRHO_ACF_ROOF(JJ) * XCPD * &
                   ( ZTAER_ROOF(JJ)/PEXNS(JJ) - PTA(JJ)/PEXNA(JJ) )
  !
  ZHEAT_RR(JJ) = PRR(JJ) * XCL * (ZTAER_ROOF(JJ) - PTA(JJ))
  !
  !*      10.     Latent heat flux between snow free roof and air
  !              -------------------------------------------------
  !
  PLEFREE_ROOF(JJ) = ZRHO_ACF_ROOF_WAT(JJ) * XLVTT * PDELT_ROOF(JJ) * &
                     ( PQSAT_ROOF(JJ) - PQA(JJ) +                     &
                       ZIMPL * ZDQSAT_ROOF(JJ) * (T%XT_ROOF(JJ,1) - ZTS_ROOF(JJ)) ) 
  !
  !      13.     Energy imbalance for verification
  !              ---------------------------------
  PIMB_ROOF(JJ) = PABS_SW_ROOF(JJ) + PABS_LW_ROOF(JJ) - PDQS_ROOF(JJ) &
               - ZDF_ROOF(JJ) * ( PHFREE_ROOF(JJ) + PLEFREE_ROOF(JJ)) &
               - PDN_ROOF(JJ) * PGSNOW_ROOF(JJ) + PFLX_BLD_ROOF(JJ)
  !
ENDDO
!
!*      11.     New saturated specified humidity near the roof surface
!              ------------------------------------------------------
!
PQSAT_ROOF(:) =  QSAT(T%XT_ROOF(:,1),PPS(:))
!
!-------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ROOF_LAYER_E_BUDGET',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------
!
END SUBROUTINE ROOF_LAYER_E_BUDGET
