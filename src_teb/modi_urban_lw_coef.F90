!auto_modi:spll_urban_lw_coef.D
MODULE MODI_URBAN_LW_COEF
INTERFACE
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
IMPLICIT NONE
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
REAL, DIMENSION(:), INTENT(IN)  :: PTS_W_A  ! wall A surface temperature
REAL, DIMENSION(:), INTENT(IN)  :: PTS_W_B  ! wall B surface temperature
REAL, DIMENSION(:), INTENT(IN)  :: PTS_R    ! road surface tempeature
REAL, DIMENSION(:), INTENT(IN)  :: PTS_G    ! garden surface temperature
REAL, DIMENSION(:), INTENT(IN)  :: PTS_WIN  ! window surface temperature
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
