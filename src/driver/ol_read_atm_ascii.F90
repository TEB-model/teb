!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL version 2.1
!SFX_LIC version 1. See LICENSE, Licence_CeCILL_V2.1-en.txt and Licence_CeCILL_V2.1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE OL_READ_ATM_ASCII (KFORC_STEP,                                 &
                              PTA,PQA,PWIND,PDIR_SW,PSCA_SW,PLW,PSNOW,    &
                              PRAIN,PPS,PCO2,PDIR                         )  
!**************************************************************************
!
!!    PURPOSE
!!    -------
!         Read in the ascii file the atmospheric forcing for the actual time
!         step KFORC_STEP, and for the next one.
!         The two time step are needed for the time interpolation of the
!         forcing.
!         If the end of the file  is reached, set the two step to the last
!         values.
!         Return undef value if the variable is not present
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      A. Lemonsu  *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     03/2008       
!
USE MODI_READ_SURF_ATM
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
! global variables
REAL, DIMENSION(:,:),INTENT(OUT) :: PTA
REAL, DIMENSION(:,:),INTENT(OUT) :: PQA
REAL, DIMENSION(:,:),INTENT(OUT) :: PWIND
REAL, DIMENSION(:,:),INTENT(OUT) :: PDIR_SW
REAL, DIMENSION(:,:),INTENT(OUT) :: PSCA_SW
REAL, DIMENSION(:,:),INTENT(OUT) :: PLW
REAL, DIMENSION(:,:),INTENT(OUT) :: PSNOW
REAL, DIMENSION(:,:),INTENT(OUT) :: PRAIN
REAL, DIMENSION(:,:),INTENT(OUT) :: PPS
REAL, DIMENSION(:,:),INTENT(OUT) :: PCO2
REAL, DIMENSION(:,:),INTENT(OUT) :: PDIR
INTEGER,INTENT(IN)               :: KFORC_STEP
! local variables
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! read data
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_ASCII',0,ZHOOK_HANDLE)
 CALL READ_SURF_ATM('ASCII ',PTA    (:,1:2),KFORC_STEP,2,122)
 CALL READ_SURF_ATM('ASCII ',PQA    (:,1:2),KFORC_STEP,2,123)
 CALL READ_SURF_ATM('ASCII ',PWIND  (:,1:2),KFORC_STEP,2,124)
 CALL READ_SURF_ATM('ASCII ',PLW    (:,1:2),KFORC_STEP,2,125)
 CALL READ_SURF_ATM('ASCII ',PDIR_SW(:,1:2),KFORC_STEP,2,126)
 CALL READ_SURF_ATM('ASCII ',PSCA_SW(:,1:2),KFORC_STEP,2,127)
 CALL READ_SURF_ATM('ASCII ',PRAIN  (:,1:2),KFORC_STEP,2,128)
 CALL READ_SURF_ATM('ASCII ',PSNOW  (:,1:2),KFORC_STEP,2,129)
 CALL READ_SURF_ATM('ASCII ',PPS    (:,1:2),KFORC_STEP,2,130)
 CALL READ_SURF_ATM('ASCII ',PDIR   (:,1:2),KFORC_STEP,2,131)
 CALL READ_SURF_ATM('ASCII ',PCO2   (:,1:2),KFORC_STEP,2,132)
!
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_ASCII',1,ZHOOK_HANDLE)
!
END SUBROUTINE OL_READ_ATM_ASCII
