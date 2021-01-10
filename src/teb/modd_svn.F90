!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!        ###############
         MODULE MODD_SV_n
!        ###############
!
!!****  *MODD_NSV* - declaration of scalar variables numbers
!!
!!    PURPOSE
!!    -------
!!       Arrays to store the per-model NSV_* values number (suffix _A denote an array)
!!
!!    AUTHOR
!!    ------
!!      P. Tulet   Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original  01/2004
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE SV_t
!
!###############################################################################
!
! variables updated for the current model
!
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CSV ! name of the scalar variables
  INTEGER    :: NSV_CHSBEG, NSV_CHSEND    !  index of first and last gas chemistry related scalar variable
  INTEGER    :: NBEQ                      ! number of chemical gas species in the surface scheme
  INTEGER    :: NSV_DSTBEG, NSV_DSTEND    ! index of first and last dust related scalar variable
  INTEGER    :: NDSTEQ                    ! number of dust related species in scalar variables list
  INTEGER    :: NSV_SLTBEG, NSV_SLTEND    ! index of first and last sea salt related scalar variable
  INTEGER    :: NSLTEQ                    ! number of sea salt related species in scalar variables list
  INTEGER    :: NSV_AERBEG, NSV_AEREND    ! index of first and last aerosol related scalar variabl
  INTEGER    :: NAEREQ                    ! number of aerosols variables

!
!
END TYPE SV_t



CONTAINS

!




SUBROUTINE SV_INIT(YSV)
TYPE(SV_t), INTENT(INOUT) :: YSV
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_SV_N:SV_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YSV%CSV)
YSV%NBEQ=0
YSV%NSV_CHSBEG=0
YSV%NSV_CHSEND=0
YSV%NSV_DSTBEG=0
YSV%NSV_DSTEND=0
YSV%NDSTEQ=0
YSV%NSV_SLTBEG=0
YSV%NSV_SLTEND=0
YSV%NSLTEQ=0
YSV%NSV_AERBEG=0
YSV%NSV_AEREND=0
YSV%NAEREQ=0
IF (LHOOK) CALL DR_HOOK("MODD_SV_N:SV_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE SV_INIT


END MODULE MODD_SV_n
