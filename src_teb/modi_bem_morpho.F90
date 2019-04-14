!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_BEM_MORPHO
INTERFACE
SUBROUTINE BEM_MORPHO(PBLD, PWALL_O_HOR, PBLD_HEIGHT, PFLOOR_HEIGHT,            &
                      PGR, PN_FLOOR, PWALL_O_BLD, PGLAZ_O_BLD, PMASS_O_BLD,     &
                      PFLOOR_HW_RATIO,                                          &
                      PF_FLOOR_MASS, PF_FLOOR_WALL, PF_FLOOR_WIN,               &
                      PF_FLOOR_ROOF, PF_WALL_FLOOR, PF_WALL_MASS,               &
                      PF_WALL_WIN, PF_WIN_FLOOR, PF_WIN_MASS, PF_WIN_WALL,      &
                      PF_MASS_FLOOR, PF_MASS_WALL, PF_MASS_WIN, PF_WASTE_CAN,   &
                      PF_WIN_WIN)
REAL, DIMENSION(:),   INTENT(IN)  :: PBLD         ! Urban horizontal building density
REAL, DIMENSION(:),   INTENT(IN)  :: PWALL_O_HOR  ! Wall to horizontal surface ratio
REAL, DIMENSION(:),   INTENT(IN)  :: PBLD_HEIGHT  ! Average building height [m]
REAL, DIMENSION(:),   INTENT(INOUT)::PFLOOR_HEIGHT ! Building floor height [m]
REAL, DIMENSION(:),   INTENT(IN)  :: PGR           ! Glazing ratio
REAL, DIMENSION(:),   INTENT(OUT) :: PN_FLOOR ! number of floor levels
REAL, DIMENSION(:),   INTENT(OUT) :: PWALL_O_BLD  ! wall surface per ground building surface [m²(wall)/m²(bld)]
REAL, DIMENSION(:),   INTENT(OUT) :: PGLAZ_O_BLD  ! glazing surface per ground building surface [m²(glazing)/m²(bld)]
REAL, DIMENSION(:),   INTENT(OUT) :: PMASS_O_BLD  ! thermal mass surface per ground building surface [m²(mass)/m²(bld)]
REAL, DIMENSION(:),   INTENT(OUT) :: PFLOOR_HW_RATIO !Height to width ratio of the floor levels
REAL, DIMENSION(:),   INTENT(OUT) :: PF_FLOOR_MASS  !view factor from floor to mass (how much floor seen by mass)
REAL, DIMENSION(:),   INTENT(OUT) :: PF_FLOOR_WALL
REAL, DIMENSION(:),   INTENT(OUT) :: PF_FLOOR_WIN
REAL, DIMENSION(:),   INTENT(OUT) :: PF_FLOOR_ROOF
REAL, DIMENSION(:),   INTENT(OUT) :: PF_WALL_FLOOR
REAL, DIMENSION(:),   INTENT(OUT) :: PF_WALL_MASS
REAL, DIMENSION(:),   INTENT(OUT) :: PF_WALL_WIN
REAL, DIMENSION(:),   INTENT(OUT) :: PF_WIN_FLOOR
REAL, DIMENSION(:),   INTENT(OUT) :: PF_WIN_MASS
REAL, DIMENSION(:),   INTENT(OUT) :: PF_WIN_WALL
REAL, DIMENSION(:),   INTENT(OUT) :: PF_WIN_WIN
REAL, DIMENSION(:),   INTENT(OUT) :: PF_MASS_FLOOR
REAL, DIMENSION(:),   INTENT(OUT) :: PF_MASS_WALL
REAL, DIMENSION(:),   INTENT(OUT) :: PF_MASS_WIN
REAL, DIMENSION(:),   INTENT(INOUT)::PF_WASTE_CAN !fraction of waste heat released into the canyon
END SUBROUTINE BEM_MORPHO
END INTERFACE
END MODULE MODI_BEM_MORPHO
