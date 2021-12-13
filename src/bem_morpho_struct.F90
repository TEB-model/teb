!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL licence
!SFX_LIC version 2.1. See Licence_CeCILL_V2.1-en.txt and Licence_CeCILL_V2.1-fr.txt  
!SFX_LIC for details.
!#####################################################################
SUBROUTINE BEM_MORPHO_STRUCT(PBLD, PWALL_O_HOR, PBLD_HEIGHT, PFLOOR_HEIGHT,     &
                      PGR, PN_FLOOR, PWALL_O_BLD, PGLAZ_O_BLD, PMASS_O_BLD,     &
                      PFLOOR_HW_RATIO,                                          &
                      PF_FLOOR_MASS, PF_FLOOR_WALL, PF_FLOOR_WIN,               &
                      PF_FLOOR_ROOF, PF_WALL_FLOOR, PF_WALL_MASS,               &
                      PF_WALL_WIN, PF_WIN_FLOOR, PF_WIN_MASS, PF_WIN_WALL,      &
                      PF_MASS_FLOOR, PF_MASS_WALL, PF_MASS_WIN, PF_WASTE_CAN,   &
                      PF_WIN_WIN)
!#####################################################################
!
!!**** *BEM_MORPHO 
!!
!!    PURPOSE
!!    -------
!! compute morphometric parameter of the building for BEM
!! 
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      G. Pigeon   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2011
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_BEM_MORPHO
!
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_TEB_n, ONLY : TEB_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
REAL, DIMENSION(:),   INTENT(IN)  :: PBLD         ! Urban horizontal building density
REAL, DIMENSION(:),   INTENT(IN)  :: PWALL_O_HOR  ! Wall to horizontal surface ratio
REAL, DIMENSION(:),   INTENT(IN)  :: PBLD_HEIGHT  ! Average building height [m]
REAL, DIMENSION(:),   INTENT(INOUT)::PFLOOR_HEIGHT ! Building floor height [m]
REAL, DIMENSION(:),   INTENT(IN)  :: PGR           ! Glazing ratio
REAL, DIMENSION(:),   INTENT(OUT) :: PN_FLOOR ! number of floor levels
REAL, DIMENSION(:),   INTENT(OUT) :: PWALL_O_BLD  ! wall surface per ground building surface [m2(wall)/m2(bld)]
REAL, DIMENSION(:),   INTENT(OUT) :: PGLAZ_O_BLD  ! glazing surface per ground building surface [m2(glazing)/m2(bld)]
REAL, DIMENSION(:),   INTENT(OUT) :: PMASS_O_BLD  ! thermal mass surface per ground building surface [m2(mass)/m2(bld)]
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
!
!*      0.2    Declarations of local variables 
!
TYPE(BEM_t) :: B
TYPE(TEB_t) :: T
!
!------------------------------------------------------------------------------
ALLOCATE(B%XFLOOR_HEIGHT(1))
ALLOCATE(B%XGR          (1))
ALLOCATE(B%XF_WASTE_CAN (1))
ALLOCATE(B%XN_FLOOR     (1))
ALLOCATE(T%XWALL_O_BLD  (1))
ALLOCATE(B%XGLAZ_O_BLD  (1))
ALLOCATE(B%XMASS_O_BLD  (1))
ALLOCATE(B%XFLOOR_HW_RATIO  (1))
ALLOCATE(B%XF_FLOOR_MASS(1))
ALLOCATE(B%XF_FLOOR_WALL(1))
ALLOCATE(B%XF_FLOOR_WIN (1))
ALLOCATE(B%XF_FLOOR_ROOF(1))
ALLOCATE(B%XF_WALL_FLOOR(1))
ALLOCATE(B%XF_WALL_MASS (1))
ALLOCATE(B%XF_WALL_WIN  (1))
ALLOCATE(B%XF_WIN_FLOOR (1))
ALLOCATE(B%XF_WIN_MASS  (1))
ALLOCATE(B%XF_WIN_WALL  (1))
ALLOCATE(B%XF_WIN_WIN   (1))
ALLOCATE(B%XF_MASS_FLOOR(1))
ALLOCATE(B%XF_MASS_WALL (1))
ALLOCATE(B%XF_MASS_WIN  (1))
!------------------------------------------------------------------------------
!
! put input variables in structures
!
B%XFLOOR_HEIGHT = PFLOOR_HEIGHT ! Building Floor height
B%XGR           = PGR           ! Glazing ratio
B%XF_WASTE_CAN  = PF_WASTE_CAN  !fraction of waste heat released into the canyon
!
!------------------------------------------------------------------------------
!
CALL BEM_MORPHO(PBLD, PWALL_O_HOR, PBLD_HEIGHT, PWALL_O_BLD, B)
!
!------------------------------------------------------------------------------
!
! put computed variables in structures
!
 PN_FLOOR        = B%XN_FLOOR ! number of floor levels
 PGLAZ_O_BLD     = B%XGLAZ_O_BLD ! glazing surface per ground building surface [m2(glazing)/m2(bld)]
 PMASS_O_BLD     = B%XMASS_O_BLD ! thermal mass surface per ground building surface [m2(mass)/m2(bld)]
 PFLOOR_HW_RATIO = B%XFLOOR_HW_RATIO !Height to width ratio of the floor levels
 PF_FLOOR_MASS   = B%XF_FLOOR_MASS !view factor from floor to mass (how much floor seen by mass)
 PF_FLOOR_WALL   = B%XF_FLOOR_WALL
 PF_FLOOR_WIN    = B%XF_FLOOR_WIN
 PF_FLOOR_ROOF   = B%XF_FLOOR_ROOF
 PF_WALL_FLOOR   = B%XF_WALL_FLOOR
 PF_WALL_MASS    = B%XF_WALL_MASS
 PF_WALL_WIN     = B%XF_WALL_WIN
 PF_WIN_FLOOR    = B%XF_WIN_FLOOR
 PF_WIN_MASS     = B%XF_WIN_MASS
 PF_WIN_WALL     = B%XF_WIN_WALL
 PF_WIN_WIN      = B%XF_WIN_WIN
 PF_MASS_FLOOR   = B%XF_MASS_FLOOR
 PF_MASS_WALL    = B%XF_MASS_WALL
 PF_MASS_WIN     = B%XF_MASS_WIN
 PF_WASTE_CAN    = B%XF_WASTE_CAN !fraction of waste heat released into the canyon
!
!------------------------------------------------------------------------------
DEALLOCATE(B%XFLOOR_HEIGHT)
DEALLOCATE(B%XGR          )
DEALLOCATE(B%XF_WASTE_CAN )
DEALLOCATE(B%XN_FLOOR     )
DEALLOCATE(T%XWALL_O_BLD  )
DEALLOCATE(B%XGLAZ_O_BLD  )
DEALLOCATE(B%XMASS_O_BLD  )
DEALLOCATE(B%XFLOOR_HW_RATIO  )
DEALLOCATE(B%XF_FLOOR_MASS)
DEALLOCATE(B%XF_FLOOR_WALL)
DEALLOCATE(B%XF_FLOOR_WIN )
DEALLOCATE(B%XF_FLOOR_ROOF)
DEALLOCATE(B%XF_WALL_FLOOR)
DEALLOCATE(B%XF_WALL_MASS )
DEALLOCATE(B%XF_WALL_WIN  )
DEALLOCATE(B%XF_WIN_FLOOR )
DEALLOCATE(B%XF_WIN_MASS  )
DEALLOCATE(B%XF_WIN_WALL  )
DEALLOCATE(B%XF_WIN_WIN   )
DEALLOCATE(B%XF_MASS_FLOOR)
DEALLOCATE(B%XF_MASS_WALL )
DEALLOCATE(B%XF_MASS_WIN  )
!------------------------------------------------------------------------------
!
END SUBROUTINE BEM_MORPHO_STRUCT
!-------------------------------------------------------------------------------
