!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_ROOF_LAYER_E_BUDGET
INTERFACE
    SUBROUTINE ROOF_LAYER_E_BUDGET(PT_ROOF, PQSAT_ROOF, PTI_BLD, PAC_BLD, PTSTEP, &
                                   HBLD, PHC_ROOF, PTC_ROOF, PD_ROOF, PDN_ROOF,   &
                                   PRHOA, PAC_ROOF, PAC_ROOF_WAT, PLW_RAD, PPS,   &
                                   PDELT_ROOF, PTA, PQA, PEXNA, PEXNS,            &
                                   PABS_SW_ROOF, PGSNOW_ROOF, PEMIS_ROOF,         &
                                   PFLX_BLD_ROOF, PDQS_ROOF, PABS_LW_ROOF,        &
                                   PHFREE_ROOF, PLEFREE_ROOF, PIMB_ROOF,          &
                                   PFRAC_GR, PG_GREENROOF_ROOF,                   &
                                   PF_FLOOR_MASS, PF_FLOOR_WALL, PF_FLOOR_WIN,    &
                                   PF_FLOOR_ROOF, PRADHT_IN,                      &
                                   PTS_MASS, PT_WIN2, PTS_FLOOR, PTI_WALL,        &
                                   PRAD_ROOF_WALL, PRAD_ROOF_WIN, PRAD_ROOF_FLOOR,&
                                   PRAD_ROOF_MASS, PCONV_ROOF_BLD, PRR,           &
                                   PLOAD_IN_ROOF                                  )
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_ROOF      ! roof layers temperatures
REAL, DIMENSION(:), INTENT(INOUT) :: PQSAT_ROOF     ! q_sat(Ts)
REAL, DIMENSION(:), INTENT(IN)    :: PTI_BLD        ! inside building temp.
REAL, DIMENSION(:), INTENT(IN)    :: PAC_BLD        ! aerodynamical resistance
REAL,               INTENT(IN)    :: PTSTEP         ! time step
 CHARACTER(LEN=3), INTENT(IN)      :: HBLD           ! Building Energy model 'DEF' or 'BEM'  
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_ROOF       ! heat capacity for roof layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_ROOF       ! thermal conductivity for roof layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_ROOF        ! depth of roof layers
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROOF       ! roof snow fraction
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA          ! air density
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROOF       ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROOF_WAT   ! aerodynamical conductance (for water)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD        ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PPS            ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_ROOF     ! fraction of water
REAL, DIMENSION(:), INTENT(IN)    :: PTA            ! air temperature at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PQA            ! air specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PEXNA          ! exner function
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS          ! surface exner function
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROOF   ! absorbed solar radiation
REAL, DIMENSION(:), INTENT(IN)    :: PGSNOW_ROOF    ! roof snow conduction
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_ROOF     ! roof emissivity
REAL, DIMENSION(:), INTENT(IN)    :: PFRAC_GR       ! fraction of green roofs
REAL, DIMENSION(:), INTENT(IN)    :: PG_GREENROOF_ROOF ! heat conduction flux
REAL, DIMENSION(:), INTENT(OUT)   :: PFLX_BLD_ROOF  ! flux from bld to roof
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_ROOF      ! heat storage inside the roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_ROOF   ! absorbed infra-red rad.
REAL, DIMENSION(:), INTENT(OUT)   :: PHFREE_ROOF    ! sensible heat flux of the
REAL, DIMENSION(:), INTENT(OUT)   :: PLEFREE_ROOF   ! latent heat flux of the
REAL, DIMENSION(:), INTENT(OUT)   :: PIMB_ROOF      ! residual energy imbalance
REAL, DIMENSION(:), INTENT(IN)    :: PF_FLOOR_MASS  ! View factor floor-mass
REAL, DIMENSION(:), INTENT(IN)    :: PF_FLOOR_WALL  ! View factor floor-wall
REAL, DIMENSION(:), INTENT(IN)    :: PF_FLOOR_WIN   ! View factor floor-window
REAL, DIMENSION(:), INTENT(IN)    :: PF_FLOOR_ROOF  ! View factor floor-roof
REAL, DIMENSION(:), INTENT(IN)    :: PRADHT_IN      ! Indoor radiant heat transfer coefficient
REAL, DIMENSION(:), INTENT(IN)    :: PTS_MASS       ! surf. mass temp. (contact with bld air)
REAL, DIMENSION(:), INTENT(IN)    :: PT_WIN2        ! indoor wind. temp.
REAL, DIMENSION(:), INTENT(IN)    :: PTS_FLOOR      ! surf. floor temp. (contact with bld air)
REAL, DIMENSION(:), INTENT(IN)    :: PTI_WALL       ! indoor wall temp.
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_ROOF_WALL ! rad. fluxes from roof to wall [W m-2(roof)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_ROOF_WIN  ! rad. fluxes from roof to win [W m-2(roof)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_ROOF_FLOOR! rad. fluxes from roof to floor [W m-2(roof)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_ROOF_MASS ! rad. fluxes from roof to mass [W m-2(roof)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_ROOF_BLD ! conv. fluxes from roof to bld [W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PRR ! rain rate [kg m-2 s-1]
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_IN_ROOF ! solar + int heat gain on roof W/m² [roof]
END SUBROUTINE ROOF_LAYER_E_BUDGET
END INTERFACE
END MODULE MODI_ROOF_LAYER_E_BUDGET
