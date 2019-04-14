!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_READ_SURF_ATM
INTERFACE

!**************************************************************************
SUBROUTINE READ_SURF_ATM     (HPROGRAM, HFILE, PFIELD,                    &
                              KFORC_STEP, KNB, KRESP, KINIT               )
!**************************************************************************
IMPLICIT NONE
!
! global variables
INTEGER,             INTENT(IN)  :: KFORC_STEP
INTEGER,             INTENT(IN)  :: KNB  
INTEGER,             INTENT(IN)  :: KRESP
INTEGER,             INTENT(IN)  :: KINIT
CHARACTER(LEN=6)    ,INTENT(IN)  :: HPROGRAM
CHARACTER(LEN=15)   ,INTENT(IN)  :: HFILE
REAL, DIMENSION(:,:),INTENT(OUT) :: PFIELD
END SUBROUTINE READ_SURF_ATM
END INTERFACE
END MODULE MODI_READ_SURF_ATM
!
!**************************************************************************
SUBROUTINE READ_SURF_ATM     (HPROGRAM, HFILE, PFIELD,                    &
                              KFORC_STEP, KNB, KRESP, KINIT               )
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
IMPLICIT NONE
!

! global variables
REAL, DIMENSION(:,:),INTENT(INOUT) :: PFIELD
INTEGER,INTENT(IN)               :: KFORC_STEP
INTEGER,INTENT(IN)               :: KNB 
INTEGER,INTENT(IN)               :: KRESP
INTEGER,INTENT(IN)               :: KINIT
CHARACTER(LEN=6)    ,INTENT(IN)  :: HPROGRAM
CHARACTER(LEN=15)   ,INTENT(IN)  :: HFILE

! local variables
INTEGER                          :: I
!
!----------------------------------------------------------------------------
!
IF (HPROGRAM == 'ASCII ') THEN
!
 IF (KFORC_STEP .EQ. 1) THEN
  REWIND(KINIT)
  DO I=1,KNB
   READ(UNIT=KINIT,FMT='(50(F20.5))') PFIELD(I,:)
  ENDDO
 ELSE
  DO I=1,KNB-1
   PFIELD(I,:) = PFIELD(I+1,:)
  ENDDO     
  READ(UNIT=KINIT,FMT='(50(F20.5))') PFIELD(KNB,:)
 ENDIF
!
ENDIF

END SUBROUTINE READ_SURF_ATM
