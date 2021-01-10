!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ##########################
      MODULE MODD_DATA_SEAFLUX_n
!     ##########################
!
!!****  *MODD_DATA_SEAFLUX - declaration of SEAFLUX surface parameters for SEAFLUX scheme
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
USE MODD_TYPE_DATE_SURF
!      
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE DATA_SEAFLUX_t
!-------------------------------------------------------------------------------
!
TYPE (DATE_TIME), POINTER, DIMENSION(:)  :: TDATA_SST  ! date of sst field 
!
!-------------------------------------------------------------------------------
!
  REAL, POINTER, DIMENSION(:,:)  :: XDATA_SST       ! sea surface temperature for
!                                                   ! each grid mesh                   (-)
!
!-------------------------------------------------------------------------------
!
  INTEGER                        :: NTIME           ! number of time data
!                                                   ! for SST
!
!-------------------------------------------------------------------------------
!
  LOGICAL                        :: LSST_DATA       ! flag to use SST data
!
!-------------------------------------------------------------------------------
END TYPE DATA_SEAFLUX_t



CONTAINS

!




SUBROUTINE DATA_SEAFLUX_INIT(YDATA_SEAFLUX)
TYPE(DATA_SEAFLUX_t), INTENT(INOUT) :: YDATA_SEAFLUX
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DATA_SEAFLUX_N:DATA_SEAFLUX_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YDATA_SEAFLUX%XDATA_SST)
YDATA_SEAFLUX%NTIME=0
YDATA_SEAFLUX%LSST_DATA=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DATA_SEAFLUX_N:DATA_SEAFLUX_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DATA_SEAFLUX_INIT


END MODULE MODD_DATA_SEAFLUX_n
