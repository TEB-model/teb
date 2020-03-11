!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL licence
!SFX_LIC version 2.1. See Licence_CeCILL_V2.1-en.txt and Licence_CeCILL_V2.1-fr.txt  
!SFX_LIC for details.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_CIRCUMSOLAR_RAD
INTERFACE
SUBROUTINE CIRCUMSOLAR_RAD(PDIR_SW, PSCA_SW, PZENITH, PF1_o_B)
REAL, DIMENSION(:), INTENT(IN) :: PDIR_SW ! direct radiation
REAL, DIMENSION(:), INTENT(IN) :: PSCA_SW ! scattered radiation
REAL, DIMENSION(:), INTENT(IN) :: PZENITH ! solar zenith angle
REAL, DIMENSION(:), INTENT(OUT) :: PF1_o_B ! coef
END SUBROUTINE CIRCUMSOLAR_RAD
END INTERFACE
END MODULE MODI_CIRCUMSOLAR_RAD
!     ###############################################################################
SUBROUTINE CIRCUMSOLAR_RAD(PDIR_SW, PSCA_SW, PZENITH, PF1_o_B)
!     ###############################################################################
!
!!****  *CIRCUMSOLAR_RAD * 
!!
!!    PURPOSE
!!    -------
!!
! compute the fraction of the scattered solar radiation that has a
! directionnal behaviour, it is from the circumsolar region
! for an horizontal surface
!!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
! Engeneering Reference, EnergyPlus DOE, 2012
! Perez, 1990 
!!
!!    AUTHOR
!!    ------
!!     G. Pigeon 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    27/11/2012
!!---------------------------------------------------------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB      
USE MODD_CSTS ,ONLY : XPI

IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN) :: PDIR_SW ! direct radiation
REAL, DIMENSION(:), INTENT(IN) :: PSCA_SW ! scattered radiation
REAL, DIMENSION(:), INTENT(IN) :: PZENITH ! solar zenith angle
REAL, DIMENSION(:), INTENT(OUT) :: PF1_o_B ! coef
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PZENITH)) :: ZE !sky clearness factor
REAL, DIMENSION(SIZE(PZENITH)) :: ZB ! b coef
REAL, DIMENSION(SIZE(PZENITH)) :: ZF1 !factors depending on ZE and Z
REAL, DIMENSION(SIZE(PZENITH)) :: ZF11, ZF13 !factors depending on ZE
REAL, PARAMETER :: ZK = 1.041 !coef. from EP Engineering Reference 
REAL, PARAMETER :: ZBMIN = 0.087 !coef. from EP Engineering Reference 
REAL(KIND=JPRB) :: ZHOOK_HANDLE
REAL, DIMENSION(SIZE(PZENITH)) :: ZELEV
INTEGER :: JJ
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_VERT_DOE',0,ZHOOK_HANDLE)

ZELEV = XPI/2. - PZENITH

PF1_o_B = 0.
WHERE(PSCA_SW > 0)
   ZE = ((PSCA_SW + PDIR_SW) / PSCA_SW + ZK * ZELEV**3) / (1 + ZK * ZELEV**3)
ELSEWHERE
   ZE = -1
ENDWHERE

DO JJ=1,SIZE(ZELEV)
   ZB(JJ) = MAX(ZBMIN, cos(ZELEV(JJ)))
ENDDO

 CALL FIJ(ZE, ZF11, ZF13)
ZF1 = ZF11 + ZF13 * ZELEV
   
PF1_o_B = ZF1 / ZB
DO JJ=1,SIZE(PF1_o_B)
   PF1_o_B(JJ) = MIN(PF1_o_B(JJ), 1.)
   PF1_o_B(JJ) = MAX(PF1_o_B(JJ), 0.)
ENDDO

IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_VERT_DOE',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
CONTAINS
SUBROUTINE FIJ(PE, PF11, PF13)
REAL, DIMENSION(:), INTENT(IN) :: PE !sky clearness factor
REAL, DIMENSION(:), INTENT(OUT):: PF11
REAL, DIMENSION(:), INTENT(OUT):: PF13

WHERE(PE < 0 )
   PF11 = 0
   PF13 = 0
ELSEWHERE(PE < 1.065)
   PF11 = -0.0083117
   PF13 = -0.0620636
ELSEWHERE(PE >= 1.065 .AND. PE < 1.230 )
   PF11 = 0.1299457
   PF13 = -0.1513752
ELSEWHERE(PE >= 1.230 .AND. PE < 1.500 )
   PF11 = 0.3296958
   PF13 = -0.2210958
ELSEWHERE(PE >= 1.500 .AND. PE < 1.950)
   PF11 = 0.5682053
   PF13 =  -0.2951290
ELSEWHERE(PE >= 1.950 .AND. PE < 2.800 )
   PF11 = 0.8730280
   PF13 = -0.3616149
ELSEWHERE(PE >= 2.800 .AND. PE < 4.500)
   PF11 = 1.1326077
   PF13 = -0.4118494
ELSEWHERE(PE >= 4.500 .AND. PE < 6.200)
   PF11 = 1.0601591
   PF13 = -0.3589221
ELSEWHERE
   PF11 = 0.6777470
   PF13 = -0.2504286
ENDWHERE

END SUBROUTINE FIJ

END SUBROUTINE CIRCUMSOLAR_RAD
