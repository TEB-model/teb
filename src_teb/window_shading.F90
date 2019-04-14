!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     #############################################################
SUBROUTINE WINDOW_SHADING(PSHGC, PSHGC_SH, O_SHADE, PALB_WALL,      &
                          PABS_WIN, PABS_WINSH, PALB_WIN, PTRAN_WIN )
!     #############################################################
!
!
!
!!****  *WINDOW_SHADING*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the radiative properties of the window in case of shading
!     devices active
!         
!     
!!**  METHOD
!!     ------
!!     
!!    apply radiative properties coming from namelist input (SHGC_SH)  
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    MODD_CST
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!   EnergyPlus Engineering Reference V7, 2012, p. 217.
!!      
!!    AUTHOR
!!    ------
!!
!!	B. Bueno           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    11/10
!!      G. Pigeon   09/12 code reorganization to take into account propertie
!!                  from window
!-------------------------------------------------------------------------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
!
REAL,    DIMENSION(:), INTENT(IN) :: PSHGC     !Window solar heat gain coefficient
REAL,    DIMENSION(:), INTENT(IN) :: PSHGC_SH  !Window + shading solar heat gain coefficient
LOGICAL, DIMENSION(:), INTENT(IN) :: O_SHADE   !use of shadings TRUE -> shadings ;
                                               !FALSE -> no shading
REAL, DIMENSION(:), INTENT(IN)  :: PALB_WALL !albedo of the wall     
REAL, DIMENSION(:), INTENT(IN)  :: PABS_WIN  !Window absorptivity
REAL, DIMENSION(:), INTENT(OUT) :: PABS_WINSH!Window absorptivity after shading
REAL, DIMENSION(:), INTENT(OUT) :: PALB_WIN  !Albedo of the ensemble window + shading
REAL, DIMENSION(:), INTENT(INOUT) :: PTRAN_WIN !Window transmitivity
!
!local variables
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('WINDOW_SHADING',0,ZHOOK_HANDLE)
!
WHERE(O_SHADE)
  PTRAN_WIN (:) = PSHGC_SH (:)
  PABS_WINSH(:) = PTRAN_WIN(:) * PABS_WIN(:)
  PALB_WIN  (:) = PALB_WALL(:)
ELSE WHERE
  PTRAN_WIN(:) = PTRAN_WIN(:)
  PABS_WINSH(:) = PABS_WIN(:)                                
  PALB_WIN  (:) = 1. - PABS_WIN(:) - PTRAN_WIN(:)
END WHERE
!
WHERE ((PABS_WINSH(:) + PTRAN_WIN(:) + PALB_WIN) > 1.)
  PALB_WIN(:) = 1. - PABS_WINSH(:) - PTRAN_WIN(:)
END WHERE
!
IF (LHOOK) CALL DR_HOOK('WINDOW_SHADING',1,ZHOOK_HANDLE)
!
END SUBROUTINE WINDOW_SHADING
