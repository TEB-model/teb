!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_SOLAR_PANEL
INTERFACE
    SUBROUTINE SOLAR_PANEL(PTSTEP, PTSUN, PRESIDENTIAL,                  &
                               PEMIT_LW_ROOF, PEMIT_LWDN_PANEL, PLW_RAD, &
                               PABS_SW_PANEL, PTA, PN_FLOOR, PFRAC_PANEL,&
                               PEMIS_PANEL, PALB_PANEL, PEFF_PANEL,      &
                               PABS_LW_PANEL, PH_PANEL, PRN_PANEL,       &
                               PTHER_PRODC_DAY,                          &
                               PTHER_PROD_PANEL, PPHOT_PROD_PANEL,       &
                               PPROD_PANEL,                              &
                               PTHER_PROD_BLD,   PPHOT_PROD_BLD,         &
                               PPROD_BLD                                 )
REAL,               INTENT(IN)  :: PTSTEP          ! time step  (s)
REAL, DIMENSION(:), INTENT(IN)  :: PTSUN           ! solar time (s since solar midnight)
REAL, DIMENSION(:), INTENT(IN)  :: PRESIDENTIAL    ! Buildings Residential use fraction        (-)
REAL, DIMENSION(:), INTENT(IN)  :: PEMIT_LW_ROOF   ! Upwards   LW flux from roof               (W/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PEMIT_LWDN_PANEL! Downwards LW flux from panel              (W/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PLW_RAD         ! Incoming Longwave radiation               (W/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PABS_SW_PANEL   ! Absorbed solar energy by the solar panel  (W/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PTA             ! Air temperature                           (K)
REAL, DIMENSION(:), INTENT(IN)  :: PN_FLOOR        ! number of floors                          (-)
REAL, DIMENSION(:), INTENT(IN)  :: PFRAC_PANEL     ! fraction   of solar panel on roofs        (-)
REAL, DIMENSION(:), INTENT(IN)  :: PEMIS_PANEL     ! emissivity of solar panel                 (-)
REAL, DIMENSION(:), INTENT(IN)  :: PALB_PANEL      ! albedo     of solar panel                 (-)
REAL, DIMENSION(:), INTENT(IN)  :: PEFF_PANEL      ! efficiency of solar panel                 (-)
REAL, DIMENSION(:), INTENT(OUT) :: PABS_LW_PANEL   ! Absorbed LW enerby by solar panel         (W/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PH_PANEL        ! Sensible heat released by the solar panel (W/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PRN_PANEL       ! Net radiation     of the solar panel      (W/m2)
REAL, DIMENSION(:), INTENT(INOUT)::PTHER_PRODC_DAY ! Present day integrated thermal production of energy    (J/m2 panel)
REAL, DIMENSION(:), INTENT(OUT) :: PTHER_PROD_PANEL! Thermal      Energy production of the solar panel      (W/m2 panel)
REAL, DIMENSION(:), INTENT(OUT) :: PPHOT_PROD_PANEL! Photovoltaic Energy production of the solar panel      (W/m2 panel)
REAL, DIMENSION(:), INTENT(OUT) :: PPROD_PANEL     ! Averaged     Energy production of the solar panel      (W/m2 panel)
REAL, DIMENSION(:), INTENT(OUT) :: PTHER_PROD_BLD  ! Thermal      Energy production of the solar panel      (W/m2 bld)
REAL, DIMENSION(:), INTENT(OUT) :: PPHOT_PROD_BLD  ! Photovoltaic Energy production of the solar panel      (W/m2 bld)
REAL, DIMENSION(:), INTENT(OUT) :: PPROD_BLD       ! Averaged     Energy production of the solar panel      (W/m2 bld)
END SUBROUTINE SOLAR_PANEL
END INTERFACE
END MODULE MODI_SOLAR_PANEL
