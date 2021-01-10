!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ##########################
      MODULE MODD_DATA_TSZ0_n
!     ##########################
!
!!****  *MODD_DATA_TSZ0 - declaration of TSZ0 surface parameters for TSZ0 scheme
!!
!!    PURPOSE
!!    -------
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
!!      P. Le Moigne  *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       09/2007
!
!*       0.   DECLARATIONS
!             ------------
!
!
    
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE DATA_TSZ0_t
!-------------------------------------------------------------------------------
INTEGER                     :: NTIME
REAL, POINTER, DIMENSION(:) :: XDATA_DTS
REAL, POINTER, DIMENSION(:) :: XDATA_DHUGRD       
!-------------------------------------------------------------------------------
END TYPE DATA_TSZ0_t



CONTAINS

!




SUBROUTINE DATA_TSZ0_INIT(YDATA_TSZ0)
TYPE(DATA_TSZ0_t), INTENT(INOUT) :: YDATA_TSZ0
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DATA_TSZ0_N:DATA_TSZ0_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YDATA_TSZ0%XDATA_DTS)
  NULLIFY(YDATA_TSZ0%XDATA_DHUGRD)
YDATA_TSZ0%NTIME=0
IF (LHOOK) CALL DR_HOOK("MODD_DATA_TSZ0_N:DATA_TSZ0_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DATA_TSZ0_INIT


END MODULE MODD_DATA_TSZ0_n
