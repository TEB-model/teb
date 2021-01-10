!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #################
      MODULE MODD_OCEAN_n
!     #################
!
!!****  *MODD_OCEAN_n - declaration of ocean varaiables 
!!                          for 1D oceanic model
!!
!!    PURPOSE
!!    -------
!     Declaration of ocean varaiables
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
!!      C. Lebeaupin   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       04/2006
!!      Modified       07/2012, P. Le Moigne : CMO1D phasing
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
TYPE OCEAN_t
!
!
!   Switche for interactive coupling with oceanic model
LOGICAL:: LMERCATOR   !set to .true. to initialize oceanic var. from Mercator
LOGICAL:: LCURRENT    !set to .true. to make initialize ocean state with current      
LOGICAL:: LPROGSST    !set to .true. to make SST evolve with tendance
INTEGER:: NTIME_COUPLING! coupling time frequency 
INTEGER:: NOCTCOUNT   !oceanic model counter
REAL :: XOCEAN_TSTEP       ! time step of the oceanic 1D model
!
! General surface: 
!
REAL, POINTER, DIMENSION(:,:) :: XSEAT  ! oceanic temperature profiles 
REAL, POINTER, DIMENSION(:,:) :: XSEAS  ! oceanic salinity profiles 
REAL, POINTER, DIMENSION(:,:) :: XSEAU  ! oceanic zonal current profiles 
REAL, POINTER, DIMENSION(:,:) :: XSEAV  ! oceanic meridian current profiles 
REAL, POINTER, DIMENSION(:,:) :: XSEAE  ! oceanic kinetic turbulent energy profiles (^(1/2))
REAL, POINTER, DIMENSION(:,:) :: XSEABATH !bathymetry indice
                                          !=1 for free sea water
                                          !=0 for sea-bed
REAL, POINTER, DIMENSION(:) ::   XSEAHMO! oceanic mixing lengths
!
REAL, POINTER, DIMENSION(:,:) :: XLE,XLK! oceanic mixing lengths
REAL, POINTER, DIMENSION(:,:) :: XKMEL,XKMELM  ! oceanic mixing coefficients
!
REAL, POINTER, DIMENSION(:) ::   XSEATEND! SST tendance
!
REAL, POINTER, DIMENSION(:,:) ::   XDTFSOL ! Temp tendancy due to solar flux
REAL, POINTER, DIMENSION(:)   ::   XDTFNSOL! -------------------- non solar flux
!
END TYPE OCEAN_t
!


CONTAINS

!




SUBROUTINE OCEAN_INIT(YOCEAN)
TYPE(OCEAN_t), INTENT(INOUT) :: YOCEAN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_OCEAN_N:OCEAN_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YOCEAN%XSEAT)
  NULLIFY(YOCEAN%XSEAS)
  NULLIFY(YOCEAN%XSEAU)
  NULLIFY(YOCEAN%XSEAV)
  NULLIFY(YOCEAN%XSEAE)
  NULLIFY(YOCEAN%XSEABATH)
  NULLIFY(YOCEAN%XSEAHMO)
  NULLIFY(YOCEAN%XLE)
  NULLIFY(YOCEAN%XLK)
  NULLIFY(YOCEAN%XKMEL)
  NULLIFY(YOCEAN%XKMELM)
  NULLIFY(YOCEAN%XSEATEND)
  NULLIFY(YOCEAN%XDTFNSOL)
  NULLIFY(YOCEAN%XDTFSOL)
YOCEAN%LMERCATOR=.FALSE.
YOCEAN%LCURRENT=.FALSE.
YOCEAN%LPROGSST=.FALSE.
YOCEAN%NTIME_COUPLING=0
YOCEAN%NOCTCOUNT=0
YOCEAN%XOCEAN_TSTEP=5*60.
IF (LHOOK) CALL DR_HOOK("MODD_OCEAN_N:OCEAN_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE OCEAN_INIT


END MODULE MODD_OCEAN_n
