!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL licence version 2.1.
! See the following links for details:
! https://cecill.info/licences/Licence_CeCILL_V2.1-en.txt
! https://cecill.info/licences/Licence_CeCILL_V2.1-fr.txt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_OL_READ_ATM_ASCII
INTERFACE

!**************************************************************************
SUBROUTINE OL_READ_ATM_ASCII (HSURF_FILETYPE, KFORC_STEP,                 &
                              PTA,PQA,PWIND,PDIR_SW,PSCA_SW,PLW,PSNOW,    &
                              PRAIN,PPS,PCO2,PDIR)
!
!**************************************************************************
IMPLICIT NONE
!
! global variables
INTEGER,             INTENT(IN)  :: KFORC_STEP
CHARACTER(LEN=6)    ,INTENT(IN)  :: HSURF_FILETYPE
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
!
END SUBROUTINE OL_READ_ATM_ASCII
END INTERFACE
END MODULE MODI_OL_READ_ATM_ASCII
!
!**************************************************************************
SUBROUTINE OL_READ_ATM_ASCII (HSURF_FILETYPE, KFORC_STEP,                 &
                              PTA,PQA,PWIND,PDIR_SW,PSCA_SW,PLW,PSNOW,    &
                              PRAIN,PPS,PCO2,PDIR)
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
!!	A. Lemonsu  *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     03/2008

!          
!
USE MODI_READ_SURF_ATM

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
!
INTEGER,INTENT(IN)               :: KFORC_STEP
CHARACTER(LEN=6)    ,INTENT(IN)  :: HSURF_FILETYPE

! local variables
INTEGER                          :: IRET
!

   
! read data
CALL READ_SURF_ATM('ASCII ','input/Forc_TA.txt    ',PTA    (1:2,:),KFORC_STEP,2,IRET,250)
CALL READ_SURF_ATM('ASCII ','input/Forc_QA.txt    ',PQA    (1:2,:),KFORC_STEP,2,IRET,251)
CALL READ_SURF_ATM('ASCII ','input/Forc_WIND.txt  ',PWIND  (1:2,:),KFORC_STEP,2,IRET,252)
CALL READ_SURF_ATM('ASCII ','input/Forc_LW.txt    ',PLW    (1:2,:),KFORC_STEP,2,IRET,253)
CALL READ_SURF_ATM('ASCII ','input/Forc_DIR_SW.txt',PDIR_SW(1:2,:),KFORC_STEP,2,IRET,254)
CALL READ_SURF_ATM('ASCII ','input/Forc_SCA_SW.txt',PSCA_SW(1:2,:),KFORC_STEP,2,IRET,255)
CALL READ_SURF_ATM('ASCII ','input/Forc_RAIN.txt  ',PRAIN  (1:2,:),KFORC_STEP,2,IRET,256)
CALL READ_SURF_ATM('ASCII ','input/Forc_SNOW.txt  ',PSNOW  (1:2,:),KFORC_STEP,2,IRET,257)
CALL READ_SURF_ATM('ASCII ','input/Forc_PS.txt    ',PPS    (1:2,:),KFORC_STEP,2,IRET,258)
CALL READ_SURF_ATM('ASCII ','input/Forc_DIR.txt   ',PDIR   (1:2,:),KFORC_STEP,2,IRET,259)
CALL READ_SURF_ATM('ASCII ','input/Forc_CO2.txt   ',PCO2   (1:2,:),KFORC_STEP,2,IRET,260)

END SUBROUTINE OL_READ_ATM_ASCII
