!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_WALL_LAYER_E_BUDGET
INTERFACE
    SUBROUTINE WALL_LAYER_E_BUDGET(HBEM, PT_WALL, PTS_WALL_B, PTI_WALL_B, PTSTEP,          &
                                   PHC_WALL, PTC_WALL, PD_WALL, PDN_ROAD, PRHOA, PAC_WALL, &
                                   PAC_BLD, PTI_BLD, PLW_RAD, PPS, PEXNS,     &
                                   PABS_SW_WALL, PT_CANYON, PTS_ROAD, PTSNOW_ROAD,         &
                                   PTS_GARDEN, PTS_MASS, PTS_FLOOR, PEMIS_WALL,            &
                                   PLW_WA_TO_WB, PLW_R_TO_W, PLW_G_TO_W, PLW_NR_TO_W,      &
                                   PLW_WIN_TO_W, PLW_S_TO_W,                               &
                                   PT_WIN1,                                                &
                                   PFLX_BLD_WALL, PDQS_WALL, PT_WIN2, PABS_LW_WALL,        &
                                   PEMIT_LW_WALL, PH_WALL, PIMB_WALL, PF_WALL_MASS,        &
                                   PF_WALL_FLOOR, PF_WALL_WIN, PRADHT_IN, PRAD_ROOF_WALL,  &
                                   PRAD_WALL_WIN, PRAD_WALL_FLOOR,                         &
                                   PRAD_WALL_MASS, PCONV_WALL_BLD, PLOAD_IN_WALL ) 
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_WALL    ! wall layers temperatures
REAL, DIMENSION(:),   INTENT(IN)  :: PTS_WALL_B ! opposite wall surface temperature
REAL, DIMENSION(:),   INTENT(IN)  :: PTI_WALL_B ! opposite wall internal temperature
REAL,               INTENT(IN)    :: PTSTEP       ! time step
 CHARACTER(LEN=3), INTENT(IN)      :: HBEM         ! Building Energy model 'DEF' or 'BEM'
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_WALL     ! heat capacity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_WALL     ! thermal conductivity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_WALL      ! depth of wall layers
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROAD     ! snow-covered fraction on roads
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA        ! rho
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WALL     ! aerodynamical conductance [m/s]
REAL, DIMENSION(:), INTENT(IN)    :: PAC_BLD      ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(IN)    :: PTI_BLD      ! inside building temperature
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD      ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PPS          ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS        ! surface Exner function
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WALL ! absorbed solar radiation
REAL, DIMENSION(:), INTENT(IN)    :: PT_CANYON    ! air canyon temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_ROAD     ! road surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTSNOW_ROAD  ! road snow temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GARDEN   ! green area surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_MASS     ! surface mass temperature  [K]
REAL, DIMENSION(:), INTENT(IN)    :: PTS_FLOOR    ! floor layers temperatures [K]
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_WALL   ! wall emissivity
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_WB ! LW interactions wall  -> opposite wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_W   ! LW interactions road -> wall ; DEF formulation 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_W   ! LW interactions garden -> wall ; DEF formulation
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_W   ! LW interactions sky   -> wall 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_W  ! LW interactions road(snow) -> wall 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WIN_TO_W ! Radiative heat trasfer coeff wall-window 
REAL, DIMENSION(:), INTENT(IN)    :: PT_WIN1      ! outdoor window temperature [K]
REAL, DIMENSION(:), INTENT(OUT)   :: PFLX_BLD_WALL! flux from bld to wall
REAL, DIMENSION(:), INTENT(INOUT) :: PDQS_WALL    ! heat storage inside the wall 
REAL, DIMENSION(:), INTENT(IN)    :: PT_WIN2      ! indoor window temperature [K]
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_WALL ! absorbed infrared rad. [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PEMIT_LW_WALL  ! LW flux emitted by the wall [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PH_WALL      ! Sensible heat flux from wall to air [W/m²(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PIMB_WALL    ! wall residual energy imbalance 
REAL, DIMENSION(:),   INTENT(IN)  :: PF_WALL_FLOOR ! View factor wall-floor
REAL, DIMENSION(:),   INTENT(IN)  :: PF_WALL_MASS  ! View factor wall-mass
REAL, DIMENSION(:),   INTENT(IN)  :: PF_WALL_WIN   ! View factor wall-win
REAL, DIMENSION(:),   INTENT(IN)  :: PRADHT_IN     ! Indoor radiant heat transfer coefficient
REAL, DIMENSION(:), INTENT(IN)    :: PRAD_ROOF_WALL ! rad. fluxes from roof to wall [W m-2(roof)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WALL_WIN  ! rad. fluxes from wall to win  [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WALL_FLOOR! rad. fluxes from wall to floor [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_WALL_MASS ! rad. fluxes from wall to mass [W m-2(wall)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_WALL_BLD ! conv. fluxes from wall to bld [W m-2(wall)]
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_IN_WALL  ! LOAD from solar heat gain + rad int. gains  W/m² [Wall]
END SUBROUTINE WALL_LAYER_E_BUDGET
END INTERFACE
END MODULE MODI_WALL_LAYER_E_BUDGET
