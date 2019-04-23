!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!   ##########################################################################
    SUBROUTINE URBAN_LW_COEF(B, T, PLW_RAD, PEMIS_G, PTS_SR, PTS_G,              &  
                             PLW_WA_TO_WB, PLW_WA_TO_R, PLW_WB_TO_R,             &
                             PLW_WA_TO_NR,PLW_WB_TO_NR,                          &
                             PLW_WA_TO_G, PLW_WB_TO_G,                           &
                             PLW_WA_TO_WIN, PLW_WB_TO_WIN,                       &
                             PLW_R_TO_WA, PLW_R_TO_WB, PLW_R_TO_WIN,             &
                             PLW_G_TO_WA, PLW_G_TO_WB, PLW_G_TO_WIN,             &
                             PLW_S_TO_WA, PLW_S_TO_WB, PLW_S_TO_R,               &
                             PLW_S_TO_NR, PLW_S_TO_G, PLW_S_TO_WIN,              &
                             PLW_WIN_TO_WA, PLW_WIN_TO_WB, PLW_WIN_TO_R,         &
                             PLW_WIN_TO_NR, PLW_WIN_TO_G,                        &
                             PLW_NR_TO_WA, PLW_NR_TO_WB, PLW_NR_TO_WIN           )
!   ##########################################################################
!
!!****  *URBAN_LW_COEF*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the coefficients before each of the temperatures in the
!     radiative budgets
!         
!     
!!**  METHOD
!     ------
!
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
!!      Original    08/09/98 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_TEB_n, ONLY : TEB_t
!
USE MODD_CSTS,ONLY : XSTEFAN
USE MODD_SURF_PAR,ONLY : XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(BEM_t), INTENT(INOUT) :: B
TYPE(TEB_t), INTENT(INOUT) :: T
!
REAL, DIMENSION(:), INTENT(IN)  :: PLW_RAD  ! incoming LW radiation
REAL, DIMENSION(:), INTENT(IN)  :: PEMIS_G  ! GARDEN area emissivity
!
REAL, DIMENSION(:), INTENT(IN)  :: PTS_G    ! garden surface temperature
REAL, DIMENSION(:), INTENT(IN)  :: PTS_SR   ! snow surface temperature
!
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WA_TO_WB! L.W. interactions wall->opposite wall
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WA_TO_R  ! L.W. interactions wall->road for road balance 
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WB_TO_R  ! L.W. interactions wall->road for road balance 
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WA_TO_NR ! L.W. interactions wall->snow for snow balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WB_TO_NR ! L.W. interactions wall->snow for snow balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WA_TO_G  ! L.W. interactions wall->GARDEN areas for garden balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WB_TO_G  ! L.W. interactions wall->GARDEN areas for garden balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WA_TO_WIN! L.W. interactions wall->win for window balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WB_TO_WIN! L.W. interactions wall->win for window balance
!
REAL, DIMENSION(:), INTENT(OUT) :: PLW_R_TO_WA ! L.W. interactions road->wall for wall balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_R_TO_WB ! L.W. interactions road->wall for wall balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_R_TO_WIN  ! L.W. interactions road->win for win balance
!
REAL, DIMENSION(:), INTENT(OUT) :: PLW_G_TO_WA ! L.W. interactions GARDEN areas->wall for wall balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_G_TO_WB ! L.W. interactions GARDEN areas->wall for wall balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_G_TO_WIN  ! L.W. interactions GARDEN areas->road for window balance
!
REAL, DIMENSION(:), INTENT(OUT) :: PLW_S_TO_WA ! L.W. interactions sky->wall for wall balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_S_TO_WB ! L.W. interactions sky->wall for wall balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_S_TO_R  ! L.W. interactions sky->road for raod balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_S_TO_NR ! L.W. interactions sky->snow for snow balance 
REAL, DIMENSION(:), INTENT(OUT) :: PLW_S_TO_G  ! L.W. interactions sky->GARDEN areas for garden balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_S_TO_WIN ! L.W. interactions sky->win for window balance
!
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WIN_TO_WA ! L.W. interactions win->wall for wall balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WIN_TO_WB ! L.W. interactions win->wall for wall balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WIN_TO_R  ! L.W. interactions win->road for road balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WIN_TO_NR ! L.W. interactions win->GARDEN areas for snow balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WIN_TO_G  ! L.W. interactions win->GARDEN areas for garden balance
!
REAL, DIMENSION(:), INTENT(OUT) :: PLW_NR_TO_WA! L.W. interactions snow(road)->wall for wall balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_NR_TO_WB! L.W. interactions snow(road)->wall for wall balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_NR_TO_WIN ! L.W. interactions snow(road)->WIN areas for window balance
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(T%XBLD))  :: ZT_S        ! sky temperature
!
REAL, DIMENSION(SIZE(T%XBLD))  :: ZEMIS_WIN
REAL, DIMENSION(SIZE(T%XBLD))  :: ZF_W_W
REAL, DIMENSION(SIZE(T%XBLD))  :: ZF_R_W
REAL, DIMENSION(SIZE(T%XBLD))  :: ZF_R_WIN
REAL, DIMENSION(SIZE(T%XBLD))  :: ZF_W_R
REAL, DIMENSION(SIZE(T%XBLD))  :: ZF_W_G
REAL, DIMENSION(SIZE(T%XBLD))  :: ZF_W_WIN
REAL, DIMENSION(SIZE(T%XBLD))  :: ZF_W_NR
REAL, DIMENSION(SIZE(T%XBLD))  :: ZF_WIN_W
REAL, DIMENSION(SIZE(T%XBLD))  :: ZF_WIN_R
REAL, DIMENSION(SIZE(T%XBLD))  :: ZF_WIN_G
REAL, DIMENSION(SIZE(T%XBLD))  :: ZF_WIN_NR
REAL, DIMENSION(SIZE(T%XBLD))  :: ZF_G_W
REAL, DIMENSION(SIZE(T%XBLD))  :: ZF_G_WIN
!
REAL :: ZE1, ZE2, ZF, ZT1, ZT2, ZLW
ZLW(ZE1,ZE2,ZF,ZT1,ZT2) = 4.*XSTEFAN*ZE1*ZE2*ZF*((ZT1+ZT2)/2.)**3
!
INTEGER :: JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('URBAN_LW_COEF',0,ZHOOK_HANDLE)
!
ZT_S(:) = (PLW_RAD(:)/XSTEFAN)**0.25
!
DO JJ=1,SIZE(T%XROAD)
  !
  ZEMIS_WIN(JJ) = 0.84 !from Energy Plus Engineering Reference, p219
    ! see http://apps1.eere.energy.gov/buildings/energyplus/
  !
  ZF_W_R(JJ)     = T%XSVF_WALL(JJ) * (T%XROAD(JJ)    /(1.-T%XBLD(JJ)))
  ZF_W_G(JJ)     = T%XSVF_WALL(JJ) * (T%XGARDEN(JJ)  /(1.-T%XBLD(JJ)))
  ZF_W_NR(JJ)    = T%XSVF_WALL(JJ)    
  !  
  ZF_WIN_R(JJ)   = ZF_W_R(JJ)
  ZF_WIN_G(JJ)   = ZF_W_G(JJ)
  ZF_WIN_NR(JJ)  = ZF_W_NR(JJ)
  !
  ZF_R_W(JJ)     = (1 - T%XSVF_ROAD(JJ))*(1.-B%XGR(JJ)) * 0.5
  ZF_R_WIN(JJ)   = (1 - T%XSVF_ROAD(JJ))*B%XGR(JJ)
  !
  ZF_G_W(JJ)     = ZF_R_W(JJ) 
  ZF_G_WIN(JJ)   = ZF_R_WIN(JJ)    
  !
  ZF_WIN_W(JJ)   = (1.-2.*T%XSVF_WALL(JJ))*(1.-B%XGR(JJ))
  ZF_W_W  (JJ)   = (1.-2.*T%XSVF_WALL(JJ))*(1.-B%XGR(JJ))
  ZF_W_WIN(JJ)   = (1.-2.*T%XSVF_WALL(JJ))*B%XGR(JJ)
  !
  !
  PLW_WA_TO_WB(JJ)=ZLW(T%XEMIS_WALL(JJ),T%XEMIS_WALL(JJ),ZF_W_W(JJ),T%XT_WALL_A(JJ,1),T%XT_WALL_B(JJ,1))
  !
  PLW_WA_TO_R(JJ) = ZLW(T%XEMIS_WALL(JJ),T%XEMIS_ROAD(JJ),ZF_R_W(JJ),T%XT_WALL_A(JJ,1),T%XT_ROAD(JJ,1))
  PLW_WB_TO_R(JJ) = ZLW(T%XEMIS_WALL(JJ),T%XEMIS_ROAD(JJ),ZF_R_W(JJ),T%XT_WALL_B(JJ,1),T%XT_ROAD(JJ,1))
  PLW_R_TO_WA(JJ) = ZLW(T%XEMIS_ROAD(JJ),T%XEMIS_WALL(JJ),ZF_W_R(JJ),T%XT_ROAD(JJ,1),T%XT_WALL_A(JJ,1))
  PLW_R_TO_WB(JJ) = ZLW(T%XEMIS_ROAD(JJ),T%XEMIS_WALL(JJ),ZF_W_R(JJ),T%XT_ROAD(JJ,1),T%XT_WALL_B(JJ,1))
  !
  IF (SIZE(PTS_G)>0) THEN
    PLW_WA_TO_G(JJ) = ZLW(T%XEMIS_WALL(JJ),PEMIS_G(JJ),ZF_G_W(JJ),T%XT_WALL_A(JJ,1),PTS_G(JJ))
    PLW_WB_TO_G(JJ) = ZLW(T%XEMIS_WALL(JJ),PEMIS_G(JJ),ZF_G_W(JJ),T%XT_WALL_B(JJ,1),PTS_G(JJ))
    PLW_G_TO_WA(JJ) = ZLW(PEMIS_G(JJ),T%XEMIS_WALL(JJ),ZF_W_G(JJ),PTS_G(JJ),T%XT_WALL_A(JJ,1))
    PLW_G_TO_WB(JJ) = ZLW(PEMIS_G(JJ),T%XEMIS_WALL(JJ),ZF_W_G(JJ),PTS_G(JJ),T%XT_WALL_B(JJ,1))
  ELSE
    PLW_WA_TO_G(JJ) = 0.
    PLW_WB_TO_G(JJ) = 0.
    PLW_G_TO_WA(JJ) = 0.
    PLW_G_TO_WB(JJ) = 0.
  ENDIF
  !
  !
  IF (PTS_SR(JJ) .EQ. XUNDEF) THEN
    PLW_WA_TO_NR (JJ) = 0.
    PLW_WB_TO_NR (JJ) = 0.
    PLW_NR_TO_WA (JJ) = 0.
    PLW_NR_TO_WB (JJ) = 0.
    !
    PLW_WIN_TO_NR(JJ) = 0.
    PLW_NR_TO_WIN(JJ) = 0
    !     
    PLW_S_TO_NR  (JJ) = 0.
  ELSE
    PLW_WA_TO_NR(JJ) = ZLW(T%XEMIS_WALL(JJ),T%TSNOW_ROAD%EMIS(JJ),ZF_R_W(JJ),T%XT_WALL_A(JJ,1),PTS_SR(JJ))
    PLW_WB_TO_NR(JJ) = ZLW(T%XEMIS_WALL(JJ),T%TSNOW_ROAD%EMIS(JJ),ZF_R_W(JJ),T%XT_WALL_B(JJ,1),PTS_SR(JJ))
    PLW_NR_TO_WA(JJ) = ZLW(T%TSNOW_ROAD%EMIS(JJ),T%XEMIS_WALL(JJ),ZF_W_NR(JJ),PTS_SR(JJ),T%XT_WALL_A(JJ,1))    
    PLW_NR_TO_WB(JJ) = ZLW(T%TSNOW_ROAD%EMIS(JJ),T%XEMIS_WALL(JJ),ZF_W_NR(JJ),PTS_SR(JJ),T%XT_WALL_B(JJ,1))    
    !
    PLW_WIN_TO_NR(JJ)= ZLW(ZEMIS_WIN(JJ),T%TSNOW_ROAD%EMIS(JJ),ZF_R_WIN(JJ),B%XT_WIN1(JJ),PTS_SR(JJ))
    PLW_NR_TO_WIN(JJ)= ZLW(T%TSNOW_ROAD%EMIS(JJ),ZEMIS_WIN(JJ),ZF_WIN_NR(JJ),PTS_SR(JJ),B%XT_WIN1(JJ))
    !
    PLW_S_TO_NR(JJ) = ZLW(1.,T%TSNOW_ROAD%EMIS(JJ),T%XSVF_ROAD(JJ),ZT_S(JJ),PTS_SR(JJ))
  ENDIF
  !
  PLW_WIN_TO_R(JJ)= ZLW(ZEMIS_WIN(JJ),T%XEMIS_ROAD(JJ),ZF_R_WIN(JJ),B%XT_WIN1(JJ),T%XT_ROAD(JJ,1))
  PLW_R_TO_WIN(JJ)= ZLW(T%XEMIS_ROAD(JJ),ZEMIS_WIN(JJ),ZF_WIN_R(JJ),T%XT_ROAD(JJ,1),B%XT_WIN1(JJ))
  !
  IF (SIZE(PTS_G)>0) THEN
    PLW_WIN_TO_G(JJ)= ZLW(ZEMIS_WIN(JJ),PEMIS_G(JJ),ZF_G_WIN(JJ),B%XT_WIN1(JJ),PTS_G(JJ))
    PLW_G_TO_WIN(JJ)= ZLW(PEMIS_G(JJ),ZEMIS_WIN(JJ),ZF_WIN_G(JJ),PTS_G(JJ),B%XT_WIN1(JJ))
  ELSE
    PLW_WIN_TO_G(JJ) = 0.
    PLW_G_TO_WIN(JJ) = 0.
  ENDIF
  !  
  PLW_WIN_TO_WA(JJ) = ZLW(ZEMIS_WIN(JJ),T%XEMIS_WALL(JJ),ZF_W_WIN(JJ),B%XT_WIN1(JJ),T%XT_WALL_A(JJ,1))
  PLW_WIN_TO_WB(JJ) = ZLW(ZEMIS_WIN(JJ),T%XEMIS_WALL(JJ),ZF_W_WIN(JJ),B%XT_WIN1(JJ),T%XT_WALL_B(JJ,1))
  PLW_WA_TO_WIN(JJ) = ZLW(T%XEMIS_WALL(JJ),ZEMIS_WIN(JJ),ZF_WIN_W(JJ),T%XT_WALL_A(JJ,1),B%XT_WIN1(JJ))
  PLW_WB_TO_WIN(JJ) = ZLW(T%XEMIS_WALL(JJ),ZEMIS_WIN(JJ),ZF_WIN_W(JJ),T%XT_WALL_B(JJ,1),B%XT_WIN1(JJ))
  !
  PLW_S_TO_WIN(JJ) = ZLW(1.,ZEMIS_WIN(JJ),T%XSVF_WALL(JJ),ZT_S(JJ),B%XT_WIN1(JJ))
  !
  PLW_S_TO_WA(JJ)  = ZLW(1.,T%XEMIS_WALL(JJ),T%XSVF_WALL(JJ),ZT_S(JJ),T%XT_WALL_A(JJ,1))
  PLW_S_TO_WB(JJ)  = ZLW(1.,T%XEMIS_WALL(JJ),T%XSVF_WALL(JJ),ZT_S(JJ),T%XT_WALL_B(JJ,1))
  !
  IF (SIZE(PTS_G)>0) THEN
    PLW_S_TO_G(JJ)  = ZLW(1.,PEMIS_G(JJ),T%XSVF_ROAD(JJ),ZT_S(JJ),PTS_G(JJ))
  ENDIF
  PLW_S_TO_R(JJ)  = ZLW(1.,T%XEMIS_ROAD(JJ),T%XSVF_ROAD(JJ),ZT_S(JJ),T%XT_ROAD(JJ,1))  
  PLW_S_TO_NR(JJ) = ZLW(1.,T%TSNOW_ROAD%EMIS(JJ),T%XSVF_ROAD(JJ),ZT_S(JJ),PTS_SR(JJ))
  !
ENDDO
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('URBAN_LW_COEF',1,ZHOOK_HANDLE)
!
END SUBROUTINE URBAN_LW_COEF

