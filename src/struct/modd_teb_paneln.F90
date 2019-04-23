!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ################
      MODULE MODD_TEB_PANEL_n
!     ################
!
!!****  *MODD_TEB_PANEL_n - declaration of solar panel parameters
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       08/2013
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE TEB_PANEL_t
!
  REAL,    POINTER, DIMENSION(:) :: XEMIS_PANEL       ! emissivity of solar panels               (-)
  REAL,    POINTER, DIMENSION(:) :: XALB_PANEL        ! albedo of soolar panels                  (-)
  REAL,    POINTER, DIMENSION(:) :: XEFF_PANEL        ! Efficiency coefficient of solar panels   (-)
  REAL,    POINTER, DIMENSION(:) :: XFRAC_PANEL       ! Fraction of SOLAR panels on roofs        (-)
  REAL,    POINTER, DIMENSION(:) :: XTHER_PRODC_DAY   ! Present day integrated thermal production of energy    (J/m2 panel)
!
END TYPE TEB_PANEL_t



CONTAINS

!




SUBROUTINE TEB_PANEL_INIT(YTEB_PANEL)
TYPE(TEB_PANEL_t), INTENT(INOUT) :: YTEB_PANEL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_PANEL_N:TEB_PANEL_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YTEB_PANEL%XEMIS_PANEL)
  NULLIFY(YTEB_PANEL%XALB_PANEL)
  NULLIFY(YTEB_PANEL%XEFF_PANEL)
  NULLIFY(YTEB_PANEL%XFRAC_PANEL)
  NULLIFY(YTEB_PANEL%XTHER_PRODC_DAY)
IF (LHOOK) CALL DR_HOOK("MODD_TEB_PANEL_N:TEB_PANEL_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_PANEL_INIT


END MODULE MODD_TEB_PANEL_n
