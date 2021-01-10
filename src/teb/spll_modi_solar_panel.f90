!auto_modi:spll_solar_panel.D
MODULE MODI_SOLAR_PANEL
INTERFACE
    SUBROUTINE SOLAR_PANEL(TPN, DMT, PTSTEP, PTSUN, PRESIDENTIAL, PEMIT_LW_ROOF, &
                           PEMIT_LWDN_PANEL, PLW_RAD, PTA, PN_FLOOR, PPROD_BLD )
USE MODD_TEB_PANEL_n, ONLY : TEB_PANEL_t
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_t
IMPLICIT NONE
TYPE(TEB_PANEL_t), INTENT(INOUT) :: TPN
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: DMT
REAL,               INTENT(IN)  :: PTSTEP          ! time step  (s)
REAL, DIMENSION(:), INTENT(IN)  :: PTSUN           ! solar time (s since solar midnight)
REAL, DIMENSION(:), INTENT(IN)  :: PRESIDENTIAL    ! Buildings Residential use fraction        (-)
REAL, DIMENSION(:), INTENT(IN)  :: PEMIT_LW_ROOF   ! Upwards   LW flux from roof               (W/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PEMIT_LWDN_PANEL! Downwards LW flux from panel              (W/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PLW_RAD         ! Incoming Longwave radiation               (W/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PTA             ! Air temperature                           (K)
REAL, DIMENSION(:), INTENT(IN)  :: PN_FLOOR        ! number of floors                          (-)
REAL, DIMENSION(:), INTENT(OUT)  :: PPROD_BLD
END SUBROUTINE SOLAR_PANEL
END INTERFACE
END MODULE MODI_SOLAR_PANEL
