!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_FLOOR_LAYER_E_BUDGET
INTERFACE
    SUBROUTINE FLOOR_LAYER_E_BUDGET(PT_FLOOR, PTSTEP, PHC_FLOOR, PTC_FLOOR, PD_FLOOR, &
                                    PFLX_BLD_FLOOR, PDQS_FLOOR, PIMB_FLOOR,           &
                                    PF_FLOOR_MASS, PF_FLOOR_WALL, PF_FLOOR_WIN,       &
                                    PF_FLOOR_ROOF, PRADHT_IN,           &
                                    PTS_MASS, PRAD_WALL_FLOOR, PRAD_ROOF_FLOOR,       &
                                    PRAD_WIN_FLOOR, PLOAD_FLOOR, PTI_BLD,             &
                                    PRAD_FLOOR_MASS, PCONV_FLOOR_BLD                  )
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_FLOOR     ! floor layers temperatures
REAL,                 INTENT(IN)    :: PTSTEP       ! time step
REAL, DIMENSION(:,:), INTENT(IN)    :: PHC_FLOOR    ! heat capacity for road layers
REAL, DIMENSION(:,:), INTENT(IN)    :: PTC_FLOOR    ! thermal conductivity for 
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_FLOOR       ! depth of road layers
REAL, DIMENSION(:),   INTENT(OUT)  :: PFLX_BLD_FLOOR !flux from building to floor
REAL, DIMENSION(:),   INTENT(OUT) :: PDQS_FLOOR !heat storage inside the floor
REAL, DIMENSION(:),   INTENT(OUT) :: PIMB_FLOOR !floor energy residual imbalance for verification
REAL, DIMENSION(:), INTENT(IN)    :: PF_FLOOR_MASS  ! View factor floor-mass
REAL, DIMENSION(:), INTENT(IN)    :: PF_FLOOR_WALL  ! View factor floor-wall
REAL, DIMENSION(:), INTENT(IN)    :: PF_FLOOR_WIN   ! View factor floor-window
REAL, DIMENSION(:), INTENT(IN)    :: PF_FLOOR_ROOF  ! View factor floor-roof
REAL, DIMENSION(:), INTENT(IN)    :: PRADHT_IN      ! Indoor radiant heat transfer coefficient
REAL, DIMENSION(:), INTENT(IN)    :: PTS_MASS  ! surf. mass temp. (contact with bld air)
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_ROOF_FLOOR ! rad. fluxes from roof to floor[W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_WALL_FLOOR ! rad. fluxes from wall to floor[W m-2(wall)]
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_WIN_FLOOR  ! rad. fluxes from win to floor[W m-2(win)]
REAL, DIMENSION(:), INTENT(IN)    :: PTI_BLD   ! indoor air temp.
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_FLOOR ! solar and internal load to the floor
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_FLOOR_MASS  ! rad. fluxes from floor to mass [W m-2(floor)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_FLOOR_BLD  ! conv. fluxes from floor to bld [W m-2(floor)]
END SUBROUTINE FLOOR_LAYER_E_BUDGET
END INTERFACE
END MODULE MODI_FLOOR_LAYER_E_BUDGET
