!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_OL_READ_ATM
INTERFACE

!**************************************************************************
SUBROUTINE OL_READ_ATM (HSURF_FILETYPE, HFORCING_FILETYPE, KFORC_STEP,    &
                        PTA,PQA,PWIND,PDIR_SW,PSCA_SW,PLW,PSNOW,PRAIN,PPS,&
                        PCO2,PDIR   )
!
!**************************************************************************
IMPLICIT NONE
!
! global variables
INTEGER,             INTENT(IN)  :: KFORC_STEP
CHARACTER(LEN=6)    ,INTENT(IN)  :: HSURF_FILETYPE
CHARACTER(LEN=6)    ,INTENT(IN)  :: HFORCING_FILETYPE
REAL, DIMENSION(:,:),INTENT(INOUT) :: PTA
REAL, DIMENSION(:,:),INTENT(INOUT) :: PQA
REAL, DIMENSION(:,:),INTENT(INOUT) :: PWIND
REAL, DIMENSION(:,:),INTENT(INOUT) :: PDIR_SW
REAL, DIMENSION(:,:),INTENT(INOUT) :: PSCA_SW
REAL, DIMENSION(:,:),INTENT(INOUT) :: PLW
REAL, DIMENSION(:,:),INTENT(INOUT) :: PSNOW
REAL, DIMENSION(:,:),INTENT(INOUT) :: PRAIN
REAL, DIMENSION(:,:),INTENT(INOUT) :: PPS
REAL, DIMENSION(:,:),INTENT(INOUT) :: PCO2
REAL, DIMENSION(:,:),INTENT(INOUT) :: PDIR
!
END SUBROUTINE OL_READ_ATM
END INTERFACE
END MODULE MODI_OL_READ_ATM
!
!**************************************************************************
SUBROUTINE OL_READ_ATM (HSURF_FILETYPE, HFORCING_FILETYPE, KFORC_STEP,    &
                        PTA,PQA,PWIND,PDIR_SW,PSCA_SW,PLW,PSNOW,PRAIN,PPS,&
                        PCO2,PDIR   )
!**************************************************************************
!
!!    PURPOSE
!!    -------
!         Read in the netcdf file the atmospheric forcing for the actual time
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
!!	F. Habets   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     06/2003
!!      P. Le Moigne 10/2004: set XCOUNT to 2 because of revised temporal loop in offline.f90:
!!                            time evolution is done at the end of isba time step so first 
!!                            isba computation is done on first forcing time step
!!      P. Le Moigne 10/2005: consistency checking between orographies read from forcing 
!!                            file and from initial file
!!      B. Decharme  01/2009: Optional, limitation of Qair (<= Qsat(tair))
!
USE MODI_OL_READ_ATM_ASCII 
!
USE MODE_THERMOS
!
IMPLICIT NONE
!
! global variables
REAL, DIMENSION(:,:),INTENT(INOUT) :: PTA
REAL, DIMENSION(:,:),INTENT(INOUT) :: PQA
REAL, DIMENSION(:,:),INTENT(INOUT) :: PWIND
REAL, DIMENSION(:,:),INTENT(INOUT) :: PDIR_SW
REAL, DIMENSION(:,:),INTENT(INOUT) :: PSCA_SW
REAL, DIMENSION(:,:),INTENT(INOUT) :: PLW
REAL, DIMENSION(:,:),INTENT(INOUT) :: PSNOW
REAL, DIMENSION(:,:),INTENT(INOUT) :: PRAIN
REAL, DIMENSION(:,:),INTENT(INOUT) :: PPS
REAL, DIMENSION(:,:),INTENT(INOUT) :: PCO2
REAL, DIMENSION(:,:),INTENT(INOUT) :: PDIR
!
INTEGER,INTENT(IN)               :: KFORC_STEP
CHARACTER(LEN=6)    ,INTENT(IN)  :: HSURF_FILETYPE
CHARACTER(LEN=6)    ,INTENT(IN)  :: HFORCING_FILETYPE
!
!
!set time variables
!
! read data
!
IF (HFORCING_FILETYPE == 'ASCII ') THEN
  CALL OL_READ_ATM_ASCII  (HSURF_FILETYPE, KFORC_STEP,                       &
                           PTA,PQA,PWIND,PDIR_SW,PSCA_SW,PLW,PSNOW,PRAIN,PPS,&
                           PCO2,PDIR              )
ENDIF
!
END SUBROUTINE OL_READ_ATM
