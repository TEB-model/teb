!auto_modi:spll_mass_layer_e_budget.D
MODULE MODI_MASS_LAYER_E_BUDGET
INTERFACE
    SUBROUTINE MASS_LAYER_E_BUDGET(PT_MASS, PTSTEP, PHC_MASS, PTC_MASS, PD_MASS, &
                                    PFLX_BLD_MASS, PDQS_MASS, PIMB_MASS,           &
                                    PF_MASS_WALL, PF_MASS_WIN,       &
                                    PF_MASS_FLOOR, PRADHT_IN,           &
                                    PRAD_WALL_MASS, PRAD_ROOF_MASS,       &
                                    PRAD_WIN_MASS, PLOAD_MASS, PTI_BLD,             &
                                    PRAD_FLOOR_MASS, PCONV_MASS_BLD                  )
IMPLICIT NONE
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
