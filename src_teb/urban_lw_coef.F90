!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!   ##########################################################################
    SUBROUTINE URBAN_LW_COEF(PGR, PBLD, PLW_RAD,                                 &
                             PEMIS_R, PSVF_R, PEMIS_W, PSVF_W, PEMIS_G,          &
                             PROAD, PGARDEN, PESNOW_R,                           &
                             PTS_SR, PTS_W_A, PTS_W_B, PTS_R, PTS_G, PTS_WIN,    &  
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
!
!
REAL, DIMENSION(:), INTENT(IN)  :: PGR      ! glazing ratio
REAL, DIMENSION(:), INTENT(IN)  :: PBLD     ! building density
REAL, DIMENSION(:), INTENT(IN)  :: PLW_RAD  ! incoming LW radiation
REAL, DIMENSION(:), INTENT(IN)  :: PEMIS_R  ! road emissivity
REAL, DIMENSION(:), INTENT(IN)  :: PSVF_R   ! road sky view factor
REAL, DIMENSION(:), INTENT(IN)  :: PEMIS_W  ! wall emissivity
REAL, DIMENSION(:), INTENT(IN)  :: PSVF_W   ! wall sky view factor
REAL, DIMENSION(:), INTENT(IN)  :: PEMIS_G  ! GARDEN area emissivity
REAL, DIMENSION(:), INTENT(IN)  :: PROAD    ! road fraction
REAL, DIMENSION(:), INTENT(IN)  :: PGARDEN  ! GARDEN area fraction
REAL, DIMENSION(:), INTENT(IN)  :: PESNOW_R ! road snow emissivity
!
REAL, DIMENSION(:), INTENT(IN)  :: PTS_W_A  ! wall A surface temperature
REAL, DIMENSION(:), INTENT(IN)  :: PTS_W_B  ! wall B surface temperature
REAL, DIMENSION(:), INTENT(IN)  :: PTS_R    ! road surface tempeature
REAL, DIMENSION(:), INTENT(IN)  :: PTS_G    ! garden surface temperature
REAL, DIMENSION(:), INTENT(IN)  :: PTS_WIN  ! window surface temperature
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
REAL, DIMENSION(SIZE(PBLD))  :: ZT_S        ! sky temperature
!
REAL, DIMENSION(SIZE(PBLD))  :: ZEMIS_WIN
REAL, DIMENSION(SIZE(PBLD))  :: ZF_W_W
REAL, DIMENSION(SIZE(PBLD))  :: ZF_R_W
REAL, DIMENSION(SIZE(PBLD))  :: ZF_R_WIN
REAL, DIMENSION(SIZE(PBLD))  :: ZF_W_R
REAL, DIMENSION(SIZE(PBLD))  :: ZF_W_G
REAL, DIMENSION(SIZE(PBLD))  :: ZF_W_WIN
REAL, DIMENSION(SIZE(PBLD))  :: ZF_W_NR
REAL, DIMENSION(SIZE(PBLD))  :: ZF_WIN_W
REAL, DIMENSION(SIZE(PBLD))  :: ZF_WIN_R
REAL, DIMENSION(SIZE(PBLD))  :: ZF_WIN_G
REAL, DIMENSION(SIZE(PBLD))  :: ZF_WIN_NR
REAL, DIMENSION(SIZE(PBLD))  :: ZF_G_W
REAL, DIMENSION(SIZE(PBLD))  :: ZF_G_WIN
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
DO JJ=1,SIZE(PROAD)
  !
  ZEMIS_WIN(JJ) = 0.84 !from Energy Plus Engineering Reference, p219
    ! see http://apps1.eere.energy.gov/buildings/energyplus/
  !
  ZF_W_R(JJ)     = PSVF_W(JJ) * (PROAD(JJ)  /(1.-PBLD(JJ)))
  ZF_W_G(JJ)     = PSVF_W(JJ) * (PGARDEN(JJ)/(1.-PBLD(JJ)))
  ZF_W_NR(JJ)    = PSVF_W(JJ)    
  !  
  ZF_WIN_R(JJ)   = ZF_W_R(JJ)
  ZF_WIN_G(JJ)   = ZF_W_G(JJ)
  ZF_WIN_NR(JJ)  = ZF_W_NR(JJ)
  !
  ZF_R_W(JJ)     = (1 - PSVF_R(JJ))*(1.-PGR(JJ)) * 0.5
  ZF_R_WIN(JJ)   = (1 - PSVF_R(JJ))*PGR(JJ)
  !
  ZF_G_W(JJ)     = ZF_R_W(JJ) 
  ZF_G_WIN(JJ)   = ZF_R_WIN(JJ)    
  !
  ZF_WIN_W(JJ)   = (1.-2.*PSVF_W(JJ))*(1.-PGR(JJ))
  ZF_W_W  (JJ)   = (1.-2.*PSVF_W(JJ))*(1.-PGR(JJ))
  ZF_W_WIN(JJ)   = (1.-2.*PSVF_W(JJ))*PGR(JJ)
  !
  !
  PLW_WA_TO_WB(JJ)=ZLW(PEMIS_W(JJ),PEMIS_W(JJ),ZF_W_W(JJ),PTS_W_A(JJ),PTS_W_B(JJ))
  !
  PLW_WA_TO_R(JJ) = ZLW(PEMIS_W(JJ),PEMIS_R(JJ),ZF_R_W(JJ),PTS_W_A(JJ),PTS_R(JJ))
  PLW_WB_TO_R(JJ) = ZLW(PEMIS_W(JJ),PEMIS_R(JJ),ZF_R_W(JJ),PTS_W_B(JJ),PTS_R(JJ))
  PLW_R_TO_WA(JJ) = ZLW(PEMIS_R(JJ),PEMIS_W(JJ),ZF_W_R(JJ),PTS_R(JJ),PTS_W_A(JJ))
  PLW_R_TO_WB(JJ) = ZLW(PEMIS_R(JJ),PEMIS_W(JJ),ZF_W_R(JJ),PTS_R(JJ),PTS_W_B(JJ))
  !
  PLW_WA_TO_G(JJ) = ZLW(PEMIS_W(JJ),PEMIS_G(JJ),ZF_G_W(JJ),PTS_W_A(JJ),PTS_G(JJ))
  PLW_WB_TO_G(JJ) = ZLW(PEMIS_W(JJ),PEMIS_G(JJ),ZF_G_W(JJ),PTS_W_B(JJ),PTS_G(JJ))
  PLW_G_TO_WA(JJ) = ZLW(PEMIS_G(JJ),PEMIS_W(JJ),ZF_W_G(JJ),PTS_G(JJ),PTS_W_A(JJ))
  PLW_G_TO_WB(JJ) = ZLW(PEMIS_G(JJ),PEMIS_W(JJ),ZF_W_G(JJ),PTS_G(JJ),PTS_W_B(JJ))
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
    PLW_WA_TO_NR(JJ) = ZLW(PEMIS_W(JJ),PESNOW_R(JJ),ZF_R_W(JJ),PTS_W_A(JJ),PTS_SR(JJ))
    PLW_WB_TO_NR(JJ) = ZLW(PEMIS_W(JJ),PESNOW_R(JJ),ZF_R_W(JJ),PTS_W_B(JJ),PTS_SR(JJ))
    PLW_NR_TO_WA(JJ) = ZLW(PESNOW_R(JJ),PEMIS_W(JJ),ZF_W_NR(JJ),PTS_SR(JJ),PTS_W_A(JJ))    
    PLW_NR_TO_WB(JJ) = ZLW(PESNOW_R(JJ),PEMIS_W(JJ),ZF_W_NR(JJ),PTS_SR(JJ),PTS_W_B(JJ))    
    !
    PLW_WIN_TO_NR(JJ)= ZLW(ZEMIS_WIN(JJ),PESNOW_R(JJ),ZF_R_WIN(JJ),PTS_WIN(JJ),PTS_SR(JJ))
    PLW_NR_TO_WIN(JJ)= ZLW(PESNOW_R(JJ),ZEMIS_WIN(JJ),ZF_WIN_NR(JJ),PTS_SR(JJ),PTS_WIN(JJ))
    !
    PLW_S_TO_NR(JJ) = ZLW(1.,PESNOW_R(JJ),PSVF_R(JJ),ZT_S(JJ),PTS_SR(JJ))
  ENDIF
  !
  PLW_WIN_TO_R(JJ)= ZLW(ZEMIS_WIN(JJ),PEMIS_R(JJ),ZF_R_WIN(JJ),PTS_WIN(JJ),PTS_R(JJ))
  PLW_R_TO_WIN(JJ)= ZLW(PEMIS_R(JJ),ZEMIS_WIN(JJ),ZF_WIN_R(JJ),PTS_R(JJ),PTS_WIN(JJ))
  !  
  PLW_WIN_TO_G(JJ)= ZLW(ZEMIS_WIN(JJ),PEMIS_G(JJ),ZF_G_WIN(JJ),PTS_WIN(JJ),PTS_G(JJ))
  PLW_G_TO_WIN(JJ)= ZLW(PEMIS_G(JJ),ZEMIS_WIN(JJ),ZF_WIN_G(JJ),PTS_G(JJ),PTS_WIN(JJ))
  !  
  PLW_WIN_TO_WA(JJ) = ZLW(ZEMIS_WIN(JJ),PEMIS_W(JJ),ZF_W_WIN(JJ),PTS_WIN(JJ),PTS_W_A(JJ))
  PLW_WIN_TO_WB(JJ) = ZLW(ZEMIS_WIN(JJ),PEMIS_W(JJ),ZF_W_WIN(JJ),PTS_WIN(JJ),PTS_W_B(JJ))
  PLW_WA_TO_WIN(JJ) = ZLW(PEMIS_W(JJ),ZEMIS_WIN(JJ),ZF_WIN_W(JJ),PTS_W_A(JJ),PTS_WIN(JJ))
  PLW_WB_TO_WIN(JJ) = ZLW(PEMIS_W(JJ),ZEMIS_WIN(JJ),ZF_WIN_W(JJ),PTS_W_B(JJ),PTS_WIN(JJ))
  !
  PLW_S_TO_WIN(JJ) = ZLW(1.,ZEMIS_WIN(JJ),PSVF_W(JJ),ZT_S(JJ),PTS_WIN(JJ))
  !
  PLW_S_TO_WA(JJ)  = ZLW(1.,PEMIS_W(JJ),PSVF_W(JJ),ZT_S(JJ),PTS_W_A(JJ))
  PLW_S_TO_WB(JJ)  = ZLW(1.,PEMIS_W(JJ),PSVF_W(JJ),ZT_S(JJ),PTS_W_B(JJ))
  !
  PLW_S_TO_G(JJ)  = ZLW(1.,PEMIS_G(JJ),PSVF_R(JJ),ZT_S(JJ),PTS_G(JJ))
  PLW_S_TO_R(JJ)  = ZLW(1.,PEMIS_R(JJ),PSVF_R(JJ),ZT_S(JJ),PTS_R(JJ))  
  PLW_S_TO_NR(JJ) = ZLW(1.,PESNOW_R(JJ),PSVF_R(JJ),ZT_S(JJ),PTS_SR(JJ))
  !
ENDDO
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('URBAN_LW_COEF',1,ZHOOK_HANDLE)
!
END SUBROUTINE URBAN_LW_COEF

