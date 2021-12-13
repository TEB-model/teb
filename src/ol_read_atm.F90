!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL version 2.1
!SFX_LIC version 1. See LICENSE, Licence_CeCILL_V2.1-en.txt and Licence_CeCILL_V2.1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE OL_READ_ATM (&
                         HSURF_FILETYPE, HFORCING_FILETYPE, KFORC_STEP,    &
                          PTA,PQA,PWIND,PDIR_SW,PSCA_SW,PLW,PSNOW,PRAIN,PPS,&
                          PCO2,PDIR,OLIMIT_QAIR                             )  
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
!!      F. Habets   *Meteo France*
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
!
!
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!USE MODD_IO_SURF_OL, ONLY : XSTART,XCOUNT,XSTRIDE,LPARTR
!         
USE MODI_OL_READ_ATM_ASCII 
!
USE MODE_THERMOS
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
! global variables
!
!
REAL, DIMENSION(:,:),INTENT(INOUT) :: PTA !K
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
INTEGER,INTENT(IN)               :: KFORC_STEP
 CHARACTER(LEN=6)    ,INTENT(IN)  :: HSURF_FILETYPE
 CHARACTER(LEN=6)    ,INTENT(IN)  :: HFORCING_FILETYPE
LOGICAL             ,INTENT(IN)  :: OLIMIT_QAIR
!
REAL, PARAMETER :: ZTAIR  = 286.0
REAL, PARAMETER :: ZPSURF = 101325.0
!
REAL, DIMENSION(SIZE(PTA,1),SIZE(PTA,2)) :: ZWORK1, ZWORK2, ZQSAT
!
INTEGER              :: JFRC, IFRC, JJ, INI, INB
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!set time variables
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM',0,ZHOOK_HANDLE)
!XSTART =KFORC_STEP
!XCOUNT =SIZE(PTA,2)
!XSTRIDE=1
!LPARTR=.TRUE.
!
! read data
!

  CALL OL_READ_ATM_ASCII  (KFORC_STEP,                       &
                           PTA,PQA,PWIND,PDIR_SW,PSCA_SW,PLW,PSNOW,PRAIN,PPS,&
                           PCO2,PDIR                                         )  

!
! Assume Qair <= Qsat_air
!
IF(OLIMIT_QAIR)THEN
!  
  INI  = SIZE(PTA,1)
  IFRC = SIZE(PTA,2)
  INB  = 0
!
  DO JFRC=1,IFRC    
     DO JJ=1,INI     
        IF(PTA(JJ,JFRC)>0.0.AND.PTA(JJ,JFRC)/=XUNDEF)THEN
           INB             = INB+1
           ZWORK1(JJ,JFRC) = PTA(JJ,JFRC)
           ZWORK2(JJ,JFRC) = PPS(JJ,JFRC)
        ELSE
           ZWORK1(JJ,JFRC) = ZTAIR
           ZWORK2(JJ,JFRC) = ZPSURF
           PTA   (JJ,JFRC) = XUNDEF
           PQA   (JJ,JFRC) = 0.0
        ENDIF
     ENDDO
  ENDDO
!  
  IF(INB==0 .AND. INI/=0)THEN
    CALL ABOR1_SFX('OL_READ_ATM: THE FORCING IS UNDEFINED')
  ENDIF
!
  ZQSAT(:,:) = QSAT(ZWORK1(:,:),ZWORK2(:,:))
!  
  PQA(:,:) = MIN(PQA(:,:),ZQSAT(:,:))
!  
ENDIF
!
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM',1,ZHOOK_HANDLE)
!
END SUBROUTINE OL_READ_ATM
