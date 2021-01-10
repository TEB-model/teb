!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ####################
      MODULE MODD_DATA_COVER_n
!     ######################
!
!!****  *MODD_DATA_COVER_n* - declaration of correspondances between surface
!!                            classes and variables, for parameters that
!!                            can change as function of physical options
!!                            (GARDENs or not).
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
!!      V. Masson    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       09/2011
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE DATA_COVER_t
!
!-----------------------------------------------------------------------------------------------------
!
REAL, DIMENSION(:,:,:), POINTER :: XDATA_WEIGHT
!
REAL, DIMENSION(:),   POINTER :: XDATA_TOWN   ! artificial surfaces fraction
REAL, DIMENSION(:),   POINTER :: XDATA_NATURE ! natural and cul. fraction
REAL, DIMENSION(:),   POINTER :: XDATA_SEA    ! sea fraction
REAL, DIMENSION(:),   POINTER :: XDATA_WATER  ! inland water fraction
REAL, DIMENSION(:,:), POINTER :: XDATA_VEGTYPE! vegetation types fractions
REAL, DIMENSION(:),   POINTER :: XDATA_GARDEN ! garden fraction
REAL, DIMENSION(:),   POINTER :: XDATA_BLD    ! building fraction in
                                              ! artificial areas
REAL, DIMENSION(:),   POINTER :: XDATA_WALL_O_HOR  ! ratio of vert. surf.
!                                                  ! over hor. surf.
!
LOGICAL                           :: LGARDEN      ! T: define urban green areas
!                                                 ! F: no urban green areas
!
INTEGER :: NYEAR        ! current year for ecoclimap2
!
!-----------------------------------------------------------------------------------------------------
!
END TYPE DATA_COVER_t



CONTAINS


!




SUBROUTINE DATA_COVER_INIT(YDATA_COVER)
TYPE(DATA_COVER_t), INTENT(INOUT) :: YDATA_COVER
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DATA_COVER_N:DATA_COVER_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YDATA_COVER%XDATA_WEIGHT)
  NULLIFY(YDATA_COVER%XDATA_TOWN)
  NULLIFY(YDATA_COVER%XDATA_NATURE)
  NULLIFY(YDATA_COVER%XDATA_SEA)
  NULLIFY(YDATA_COVER%XDATA_WATER)
  NULLIFY(YDATA_COVER%XDATA_VEGTYPE)
  NULLIFY(YDATA_COVER%XDATA_GARDEN)
  NULLIFY(YDATA_COVER%XDATA_BLD)
  NULLIFY(YDATA_COVER%XDATA_WALL_O_HOR)
YDATA_COVER%LGARDEN=.FALSE.
YDATA_COVER%NYEAR=9999
IF (LHOOK) CALL DR_HOOK("MODD_DATA_COVER_N:DATA_COVER_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DATA_COVER_INIT


END MODULE MODD_DATA_COVER_n

