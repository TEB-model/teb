!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!#####################################################################
SUBROUTINE BEM_MORPHO(PBLD, PWALL_O_HOR, PBLD_HEIGHT, PFLOOR_HEIGHT,            &
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
REAL, DIMENSION(SIZE(PBLD)):: ZF_AUX2       ! Auxiliar view factor
!
!RJ: missing declaration
INTEGER :: JJ
!
!------------------------------------------------------------------------------
!
!*    1.   Verify floor height
!          -------------------
!
DO JJ=1,SIZE(PBLD)
  IF (PFLOOR_HEIGHT(JJ) < 2.5) THEN
!    WRITE(ILUOUT,*) 'WARNING: PFLOOR_HEIGHT low ',PFLOOR_HEIGHT(JJ),' grid mesh number ',JJ
    PFLOOR_HEIGHT(JJ) = 2.5
  ENDIF
  IF (PFLOOR_HEIGHT(JJ) > PBLD_HEIGHT(JJ)) THEN
!    WRITE(ILUOUT,*) 'WARNING: PFLOOR_HEIGHT higher than PBLD_HEIGHT ',PFLOOR_HEIGHT(JJ),' grid mesh number ',JJ, &
!                    'set to PBLD_HEIGHT'
    PFLOOR_HEIGHT(JJ) = PBLD_HEIGHT(JJ)
  ENDIF
!
!*    1.   Verify location of HVAC release for very extended buildings
!          -----------------------------------------------------------
!
  IF (PBLD(JJ) > 0.9999 .AND. PF_WASTE_CAN(JJ) > 0.) THEN
!    WRITE(ILUOUT,*) 'WARNING: PF_WASTE_CAN cannot be higher than 0. while PBLD is 0.9999',PF_WASTE_CAN(JJ), &
!                    ' grid mesh number ',JJ,' set to 0.'
    PF_WASTE_CAN(JJ) = 0.
  ENDIF
END DO
!
!*       1.    floor number, relative surf. and view factors caculation
!              --------------------------------------------------------
!
PN_FLOOR(:) = FLOAT(NINT(PBLD_HEIGHT(:) / PFLOOR_HEIGHT(:) ))
PWALL_O_BLD(:) = PWALL_O_HOR(:) * (1. - PGR(:)) / PBLD(:) ! [m2(wall)/m2(bld)]
PGLAZ_O_BLD(:) = PWALL_O_HOR(:) * PGR(:)        / PBLD(:) ! [m2(win)/m2(bld)]

WHERE(PN_FLOOR > 1.5)
   !* more than 1 floor level -> the floor and the roof don't see each other
   !
   !           ROOF
   !  |---------------------|
   !  |                     |
   !  |                     |
   !  |        MASS         |
   !  |---------------------|
   !  |---------------------|
   !  |        MASS         |
   !  |                     |
   !  |                     |
   !  |        FLOOR        |
   !  |---------------------|
   !
   PMASS_O_BLD(:) = 2 * (PN_FLOOR(:) - 1.)                     ! [m2(mass)/m2(bld)]
   PFLOOR_HW_RATIO(:) = PWALL_O_HOR(:) / PN_FLOOR(:) / PBLD(:) / 2.
   !
   PF_FLOOR_MASS(:) = (PFLOOR_HW_RATIO(:)**2 + 1.)**0.5  - PFLOOR_HW_RATIO(:)
   PF_FLOOR_WALL(:) = (1. - PF_FLOOR_MASS(:)) * (1. - PGR(:))
   PF_FLOOR_WIN (:) = (1. - PF_FLOOR_MASS(:)) * PGR(:)
   PF_FLOOR_ROOF(:) = 0. ! no rad interaction between floor and roof
   !
   ZF_AUX2      (:) = (1. - PF_FLOOR_MASS(:)) / PFLOOR_HW_RATIO(:)
   !
   PF_WALL_FLOOR(:) = ZF_AUX2(:) / ( 2.*PN_FLOOR(:) )
   PF_WALL_MASS (:) = PF_WALL_FLOOR(:) * (2.*PN_FLOOR(:)-2.) 
   PF_WALL_WIN  (:) = (1. - ZF_AUX2(:)) * PGR(:)
   !
   PF_WIN_FLOOR(:) = PF_WALL_FLOOR(:)
   PF_WIN_MASS (:) = PF_WALL_MASS (:)
   PF_WIN_WALL (:) = (1. - PF_WIN_FLOOR(:) * 2. - PF_WIN_MASS(:) ) * (1. - PGR(:))
   PF_WIN_WIN  (:) = (1. - PF_WIN_FLOOR(:) * 2. - PF_WIN_MASS(:) ) * PGR(:)
   !
   PF_MASS_FLOOR(:) = PF_FLOOR_MASS(:) / PMASS_O_BLD(:)
   PF_MASS_WALL (:) = PWALL_O_BLD(:) * PF_WALL_MASS(:) / PMASS_O_BLD(:)
   PF_MASS_WIN  (:) = PGLAZ_O_BLD (:) * PF_WIN_MASS(:)  / PMASS_O_BLD(:)
ELSE WHERE
   !* 1 floor level -> no mass view factors ; roof and floor see each other
   !
   !           ROOF
   !  |---------------------|
   !  |                     |
   !  |                     |
   !  |        FLOOR        |
   !  |---------------------|
   !
   PMASS_O_BLD(:) = 0.
   PFLOOR_HW_RATIO(:) = PWALL_O_HOR(:) / PBLD(:) / 2.
   !
   PF_FLOOR_ROOF(:) = (PFLOOR_HW_RATIO(:)**2 + 1.)**0.5  - PFLOOR_HW_RATIO(:)
   PF_FLOOR_MASS(:) = 0.
   PF_FLOOR_WALL(:) = (1. - PF_FLOOR_ROOF(:)) * (1. - PGR(:))
   PF_FLOOR_WIN (:) = (1. - PF_FLOOR_ROOF(:)) * PGR(:)
   !
   ZF_AUX2      (:) = (1. - PF_FLOOR_ROOF(:)) / PFLOOR_HW_RATIO(:)
   !
   PF_WALL_FLOOR(:) = ZF_AUX2(:) / 2.
   PF_WALL_MASS (:) = 0.
   PF_WALL_WIN  (:) = (1. - ZF_AUX2(:)) * PGR(:)
   !
   PF_WIN_FLOOR(:) = PF_WALL_FLOOR(:)
   PF_WIN_MASS (:) = 0.
   PF_WIN_WALL (:) = (1. - PF_WIN_FLOOR(:) * 2. - PF_WIN_MASS(:) ) * (1. - PGR(:))
   PF_WIN_WIN  (:) = (1. - PF_WIN_FLOOR(:) * 2. - PF_WIN_MASS(:) ) * PGR(:)
   !
   PF_MASS_FLOOR(:) = 0.
   PF_MASS_WALL (:) = 0.
   PF_MASS_WIN  (:) = 0.
END WHERE
!
END SUBROUTINE BEM_MORPHO
!-------------------------------------------------------------------------------
