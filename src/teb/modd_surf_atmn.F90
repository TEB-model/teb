!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ####################
      MODULE MODD_SURF_ATM_n
!     ######################
!
!!****  *MODD_SURF_ATM - declaration of surface parameters
!!
!!    PURPOSE
!!    -------
!     Declaration of surface parameters
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson and A. Boone   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_DATE_SURF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE SURF_ATM_t
!
!-----------------------------------------------------------------------------------------------------
!
! Type of each surface scheme
!
  CHARACTER(LEN=6) :: CTOWN      ! name of the urban surface scheme
  CHARACTER(LEN=6) :: CNATURE    ! name of the soil&vegetation surface scheme
  CHARACTER(LEN=6) :: CWATER     ! name of the scheme for inland water
  CHARACTER(LEN=6) :: CSEA       ! name for the ocean scheme
!
!-----------------------------------------------------------------------------------------------------
!
! Surface/Tile Fractions:
!
  REAL, POINTER, DIMENSION(:)   :: XTOWN     ! urban surface fraction of the grid box   (-)
  REAL, POINTER, DIMENSION(:)   :: XNATURE   ! natural surface fraction of the grid box (-)
  REAL, POINTER, DIMENSION(:)   :: XWATER    ! inland water fraction of the grid box    (-)
  REAL, POINTER, DIMENSION(:)   :: XSEA      ! sea/ocean fraction of the grid box       (-)
!
!-------------------------------------------------------------------------------
!
! type of initialization of vegetation: from cover types (ecoclimap) or parameters prescribed
!
  LOGICAL                        :: LECOCLIMAP ! T: parameters computed from ecoclimap
!                                              ! F: they are read in the file
!
  LOGICAL                        :: LECOSG     ! T: parameters computed from ecosg
!-------------------------------------------------------------------------------
!
! change water (not lake) to nature and/or town to rock : arrange cover properly
!
  LOGICAL                        :: LWATER_TO_NATURE ! T: Change Wetland treated as inland water into nature 
  LOGICAL                        :: LTOWN_TO_ROCK    ! T: Change Town into Rock
!
!-------------------------------------------------------------------------------
!
! include urban green areas for urbanized covers
!
  LOGICAL                        :: LGARDEN    ! T: define urban green areas
!                                              ! F: no urban green areas
!
!-----------------------------------------------------------------------------------------------------
!
! Masks and number of grid elements for each tile surface
!
! Sea/Ocean:
!
  INTEGER                               :: NSIZE_SEA    ! number of grid points by proc containing a
!                                                     ! sea surface                              (-)
  INTEGER                               :: NDIM_SEA     ! total number of grid points containing a
!                                                     ! sea surface                             (-)
  INTEGER, POINTER, DIMENSION(:)    :: NR_SEA       ! sea/ocean surface mask                  (-)
!
! Inland Water:
!
  INTEGER                               :: NSIZE_WATER  ! number of grid points containing an 
!                                                     ! inland water surface                    (-)
  INTEGER                               :: NDIM_WATER   ! total number of grid points by proc containing an
!                                                     ! inland surface
  INTEGER, POINTER, DIMENSION(:)    :: NR_WATER
!
! Town:
!
  INTEGER                               :: NSIZE_TOWN   ! number of grid points by proc containing an 
!                                                     ! urban surface                           (-)
  INTEGER                               :: NDIM_TOWN    ! total number of grid points containing an
!                                                     ! urban surface
  INTEGER, POINTER, DIMENSION(:)    :: NR_TOWN      ! urban surface mask                      (-)
!
! Natural surface:
!
  INTEGER                               :: NSIZE_NATURE ! number of grid points by proc containing a 
!                                                     ! natural surface                         (-)
  INTEGER                               :: NDIM_NATURE  ! total number of grid points containing a
!                                                     ! natural surface                         (-)
  INTEGER, POINTER, DIMENSION(:)    :: NR_NATURE    ! natural surface mask                    (-)
!
! All surfaces:
!
  INTEGER                               :: NSIZE_FULL   ! total number of grid points by proc     (-)
  INTEGER                               :: NDIM_FULL    ! total number of grid points             (-)
!
! global sizes of the physical domain
!
  INTEGER                               :: NIMAX_SURF_ll
  INTEGER                               :: NJMAX_SURF_ll
!
! local sizes of the physical domain
!
  INTEGER                               :: NIMAX_SURF_LOC
  INTEGER                               :: NJMAX_SURF_LOC
!
!-----------------------------------------------------------------------------------------------------
!
! Surface fields (only 1 horizontal dimension)
!
  REAL, POINTER, DIMENSION(:,:) :: XCOVER    ! fraction of each ecosystem for each grid box (-)
  LOGICAL, POINTER, DIMENSION(:):: LCOVER    ! GCOVER(i)=T --> ith cover field is not 0.
  REAL, POINTER, DIMENSION(:)   :: XZS       ! orography                                    (m)
!
!-------------------------------------------------------------------------------
!
  TYPE (DATE_TIME)                      :: TTIME            ! current date and time
!
  REAL                                  :: XOUT_TSTEP       ! output writing time step
!
!-----------------------------------------------------------------------------------------------------
!
! physical fields need into the restart file for ARPEGE/ALADIN run
!
  REAL, POINTER, DIMENSION(:)   :: XRAIN    ! Rainfall rate at surface               (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XSNOW    ! snowfall rate at surface               (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XZ0      ! surface roughness length for momentum  (m)
  REAL, POINTER, DIMENSION(:)   :: XZ0H     ! surface roughness length for heat      (m)
  REAL, POINTER, DIMENSION(:)   :: XQSURF   ! specific humidity at surface           (kg/kg)
!
!-----------------------------------------------------------------------------------------------------
!
!
END TYPE SURF_ATM_t
!


CONTAINS


!


!
!
SUBROUTINE SURF_ATM_INIT(YSURF_ATM)
TYPE(SURF_ATM_t), INTENT(INOUT) :: YSURF_ATM
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_SURF_ATM_N:SURF_ATM_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YSURF_ATM%XTOWN)
  NULLIFY(YSURF_ATM%XNATURE)
  NULLIFY(YSURF_ATM%XWATER)
  NULLIFY(YSURF_ATM%XSEA)
  NULLIFY(YSURF_ATM%NR_SEA)
  NULLIFY(YSURF_ATM%NR_WATER)
  NULLIFY(YSURF_ATM%NR_TOWN)
  NULLIFY(YSURF_ATM%NR_NATURE)
  NULLIFY(YSURF_ATM%XCOVER)
  NULLIFY(YSURF_ATM%LCOVER)
  NULLIFY(YSURF_ATM%XZS)
  NULLIFY(YSURF_ATM%XRAIN)
  NULLIFY(YSURF_ATM%XSNOW)
  NULLIFY(YSURF_ATM%XZ0)
  NULLIFY(YSURF_ATM%XZ0H)
  NULLIFY(YSURF_ATM%XQSURF)
YSURF_ATM%CTOWN=' '
YSURF_ATM%CNATURE=' '
YSURF_ATM%CWATER=' '
YSURF_ATM%CSEA=' '
YSURF_ATM%LECOCLIMAP=.FALSE.
YSURF_ATM%LECOSG=.FALSE.
YSURF_ATM%LWATER_TO_NATURE=.FALSE.
YSURF_ATM%LTOWN_TO_ROCK=.FALSE.
YSURF_ATM%LGARDEN=.FALSE.
YSURF_ATM%NSIZE_SEA=0
YSURF_ATM%NDIM_SEA=0
YSURF_ATM%NSIZE_WATER=0
YSURF_ATM%NDIM_WATER=0
YSURF_ATM%NSIZE_TOWN=0
YSURF_ATM%NDIM_TOWN=0
YSURF_ATM%NSIZE_NATURE=0
YSURF_ATM%NDIM_NATURE=0
YSURF_ATM%NSIZE_FULL=0
YSURF_ATM%NDIM_FULL=0
YSURF_ATM%NIMAX_SURF_ll=0
YSURF_ATM%NJMAX_SURF_ll=0
YSURF_ATM%NIMAX_SURF_LOC=0
YSURF_ATM%NJMAX_SURF_LOC=0
YSURF_ATM%XOUT_TSTEP=0.
IF (LHOOK) CALL DR_HOOK("MODD_SURF_ATM_N:SURF_ATM_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE SURF_ATM_INIT



END MODULE MODD_SURF_ATM_n

