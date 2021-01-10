!auto_modi:spll_urban_lw_coef.D
MODULE MODI_URBAN_LW_COEF
INTERFACE
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
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_TEB_n, ONLY : TEB_t
IMPLICIT NONE
TYPE(BEM_t), INTENT(INOUT) :: B
TYPE(TEB_t), INTENT(INOUT) :: T
REAL, DIMENSION(:), INTENT(IN)  :: PLW_RAD  ! incoming LW radiation
REAL, DIMENSION(:), INTENT(IN)  :: PEMIS_G  ! GARDEN area emissivity
REAL, DIMENSION(:), INTENT(IN)  :: PTS_G    ! garden surface temperature
REAL, DIMENSION(:), INTENT(IN)  :: PTS_SR   ! snow surface temperature
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WA_TO_WB! L.W. interactions wall->opposite wall
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WA_TO_R  ! L.W. interactions wall->road for road balance 
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WB_TO_R  ! L.W. interactions wall->road for road balance 
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WA_TO_NR ! L.W. interactions wall->snow for snow balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WB_TO_NR ! L.W. interactions wall->snow for snow balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WA_TO_G  ! L.W. interactions wall->GARDEN areas for garden balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WB_TO_G  ! L.W. interactions wall->GARDEN areas for garden balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WA_TO_WIN! L.W. interactions wall->win for window balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WB_TO_WIN! L.W. interactions wall->win for window balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_R_TO_WA ! L.W. interactions road->wall for wall balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_R_TO_WB ! L.W. interactions road->wall for wall balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_R_TO_WIN  ! L.W. interactions road->win for win balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_G_TO_WA ! L.W. interactions GARDEN areas->wall for wall balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_G_TO_WB ! L.W. interactions GARDEN areas->wall for wall balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_G_TO_WIN  ! L.W. interactions GARDEN areas->road for window balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_S_TO_WA ! L.W. interactions sky->wall for wall balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_S_TO_WB ! L.W. interactions sky->wall for wall balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_S_TO_R  ! L.W. interactions sky->road for raod balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_S_TO_NR ! L.W. interactions sky->snow for snow balance 
REAL, DIMENSION(:), INTENT(OUT) :: PLW_S_TO_G  ! L.W. interactions sky->GARDEN areas for garden balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_S_TO_WIN ! L.W. interactions sky->win for window balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WIN_TO_WA ! L.W. interactions win->wall for wall balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WIN_TO_WB ! L.W. interactions win->wall for wall balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WIN_TO_R  ! L.W. interactions win->road for road balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WIN_TO_NR ! L.W. interactions win->GARDEN areas for snow balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_WIN_TO_G  ! L.W. interactions win->GARDEN areas for garden balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_NR_TO_WA! L.W. interactions snow(road)->wall for wall balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_NR_TO_WB! L.W. interactions snow(road)->wall for wall balance
REAL, DIMENSION(:), INTENT(OUT) :: PLW_NR_TO_WIN ! L.W. interactions snow(road)->WIN areas for window balance
END SUBROUTINE URBAN_LW_COEF
END INTERFACE
END MODULE MODI_URBAN_LW_COEF
