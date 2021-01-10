!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ################
      MODULE MODD_BEM_OPTION_n
!     ################
!
!!****  *MODD_BEM_n - declaration of parameters and option for BEM
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
!!      B. Bueno   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       10/2010
!!      G. Pigeon      06/2011 add LSHAD_DAY
!!      G. Pigeon      07/2011 add LNATVENT_NIGHT
!!      G. Pigeon      08/2011 change from MODD_BLD -> MODD_BEM
!!      G. Pigeon      10/2011 add indoor relative surf. and view factors
!!      G. Pigeon      09/2012 add TRAN_WIN
!!      G. Pigeon      10/2012 add XF_WIN_WIN
!!      V. Masson      06/2013 splits module in two
!
!*       0.   DECLARATIONS
!             ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE BEM_OPTIONS_t
! BLD scheme option
!
! Number of layers
!
  INTEGER                       :: NFLOOR_LAYER   ! number of layers in walls
  CHARACTER(LEN=6)              :: CCOOL_COIL    ! type of cooling coil
  CHARACTER(LEN=6)              :: CHEAT_COIL    ! type of heating coil
  LOGICAL                       :: LAUTOSIZE     ! Flag to activate autosize calculations
!
END TYPE BEM_OPTIONS_t
!


CONTAINS

!




SUBROUTINE BEM_OPTIONS_INIT(YBEM_OPTIONS)
TYPE(BEM_OPTIONS_t), INTENT(INOUT) :: YBEM_OPTIONS
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_BEM_N:BEM_OPTIONS_INIT",0,ZHOOK_HANDLE)
YBEM_OPTIONS%NFLOOR_LAYER = 0
YBEM_OPTIONS%CCOOL_COIL   = '      '
YBEM_OPTIONS%CHEAT_COIL   = '      '
YBEM_OPTIONS%LAUTOSIZE    = .FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_BEM_N:BEM_OPTIONS_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE BEM_OPTIONS_INIT


!----------------------------------------------------------------------------
!
END MODULE MODD_BEM_OPTION_n
