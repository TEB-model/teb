!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ##################
      MODULE MODD_OCEAN_REL_n
!     ##################
!
!     Author : P. Peyrille
!!    Date: 01/2012
!! 
!!  declaration of relaxation field, flag for ocean model
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
TYPE OCEAN_REL_t
!
 REAL                :: XTAU_REL   ! Time of relaxation  (s)
!
 REAL                :: XQCORR     ! correction coefficient  for surface fluxes  (BArnier et al. 1998) typcial val: 100.
!
 LOGICAL                            :: LREL_CUR   ! flag for  relxation on current
 LOGICAL                            :: LREL_TS    ! flag for relaxation on ocean temperature
 LOGICAL                            :: LFLUX_NULL ! Flag for testing zero  incoming flux at the surface 
!
 LOGICAL                            :: LFLX_CORR  ! Flag for  flux  correction 
 LOGICAL                            :: LDIAPYCNAL ! Flag for diapycnal mixing activation
!
 REAL, POINTER, DIMENSION(:,:) :: XSEAU_REL  ! ref. U current  profile (C)
 REAL, POINTER, DIMENSION(:,:) :: XSEAV_REL  ! ref. V current  profile (C)
 REAL, POINTER, DIMENSION(:,:) :: XSEAT_REL  ! ref. temperature  profile (C)
 REAL, POINTER, DIMENSION(:,:) :: XSEAS_REL  ! ---- salinity   ---------- (%)
!
END TYPE OCEAN_REL_t
!


CONTAINS
!
!


!
!
SUBROUTINE OCEAN_REL_INIT(YOCEAN_REL)
TYPE(OCEAN_REL_t), INTENT(INOUT) :: YOCEAN_REL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_OCEAN_REL_N:OCEAN_REL_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YOCEAN_REL%XSEAT_REL)
  NULLIFY(YOCEAN_REL%XSEAS_REL)
  NULLIFY(YOCEAN_REL%XSEAU_REL)
  NULLIFY(YOCEAN_REL%XSEAV_REL)
YOCEAN_REL%XTAU_REL=0.
YOCEAN_REL%XQCORR=0.
YOCEAN_REL%LREL_CUR=.FALSE.
YOCEAN_REL%LREL_TS=.FALSE.
YOCEAN_REL%LFLUX_NULL=.FALSE.
YOCEAN_REL%LFLX_CORR=.FALSE.
YOCEAN_REL%LDIAPYCNAL=.FALSE.
!
IF (LHOOK) CALL DR_HOOK("MODD_OCEAN_REL_N:OCEAN_REL_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE OCEAN_REL_INIT
!
!
END MODULE MODD_OCEAN_REL_n
