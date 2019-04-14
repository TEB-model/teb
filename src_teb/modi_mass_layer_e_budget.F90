!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_MASS_LAYER_E_BUDGET
INTERFACE
    SUBROUTINE MASS_LAYER_E_BUDGET(PT_MASS, PTSTEP, PHC_MASS, PTC_MASS, PD_MASS, &
                                    PFLX_BLD_MASS, PDQS_MASS, PIMB_MASS,           &
                                    PF_MASS_WALL, PF_MASS_WIN,       &
                                    PF_MASS_FLOOR, PRADHT_IN,           &
                                    PRAD_WALL_MASS, PRAD_ROOF_MASS,       &
                                    PRAD_WIN_MASS, PLOAD_MASS, PTI_BLD,             &
                                    PRAD_FLOOR_MASS, PCONV_MASS_BLD                  )
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_MASS     ! floor layers temperatures
REAL,                 INTENT(IN)    :: PTSTEP       ! time step
REAL, DIMENSION(:,:), INTENT(IN)    :: PHC_MASS    ! heat capacity for road layers
REAL, DIMENSION(:,:), INTENT(IN)    :: PTC_MASS    ! thermal conductivity for 
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_MASS       ! depth of road layers
REAL, DIMENSION(:),   INTENT(OUT)  :: PFLX_BLD_MASS !flux from building to floor
REAL, DIMENSION(:),   INTENT(OUT) :: PDQS_MASS !heat storage inside the floor
REAL, DIMENSION(:),   INTENT(OUT) :: PIMB_MASS !floor energy residual imbalance for verification
REAL, DIMENSION(:), INTENT(IN)    :: PF_MASS_WALL  ! View factor mass-wall
REAL, DIMENSION(:), INTENT(IN)    :: PF_MASS_WIN   ! View factor mass-window
REAL, DIMENSION(:), INTENT(IN)    :: PF_MASS_FLOOR ! View factor mass-floor
REAL, DIMENSION(:), INTENT(IN)    :: PRADHT_IN      ! Indoor radiant heat transfer coefficient
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_ROOF_MASS ! rad. fluxes from roof to floor[W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_WALL_MASS ! rad. fluxes from wall to floor[W m-2(wall)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_WIN_MASS  ! rad. fluxes from win to floor[W m-2(win)]
REAL, DIMENSION(:), INTENT(IN)    :: PTI_BLD   ! indoor air temp.
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_MASS ! solar and internal load to the floor
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_FLOOR_MASS  ! rad. fluxes from floor to mass [W m-2(floor)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_MASS_BLD  ! conv. fluxes from floor to bld [W m-2(floor)]
END SUBROUTINE MASS_LAYER_E_BUDGET
END INTERFACE
END MODULE MODI_MASS_LAYER_E_BUDGET
