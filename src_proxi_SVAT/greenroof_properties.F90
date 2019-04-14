!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     #########
      SUBROUTINE GREENROOF_PROPERTIES(PDIR_SW, PSCA_SW, PSW_BANDS, KSW,      &
                                     PTS, PEMIS, PALB, PTA                  )  
!     ##########################################################################
!
!!****  *GREENROOF_PROPERTIES*  
!!
!!    PURPOSE
!!    -------
!
!!
!!
!!========================================================================
!!========================================================================
!!========================================================================
!!
!! ==> YOU ARE (MORE THAN) WELCOME TO USE YOUR OWN VEGETATION DATA HERE
!!
!!========================================================================
!!========================================================================
!!========================================================================
!
!         
!!    EXTERNAL
!!    --------
!!
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      
!!    AUTHOR
!!    ------
!!
!!	V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original                           
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:),   INTENT(IN)   :: PTA                ! air temperature
REAL, DIMENSION(:,:), INTENT(IN)   :: PDIR_SW            ! direct incoming solar radiation
REAL, DIMENSION(:,:), INTENT(IN)   :: PSCA_SW            ! diffus incoming solar radiation
REAL, DIMENSION(:)  , INTENT(IN)   :: PSW_BANDS          ! mean wavelength of each shortwave band (m)
INTEGER,              INTENT(IN)   :: KSW                ! number of short-wave spectral bands
!
REAL, DIMENSION(:)  , INTENT(OUT)  :: PTS                ! radiative surface temperature
REAL, DIMENSION(:)  , INTENT(OUT)  :: PEMIS              ! green areas emissivity
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALB               ! green areas albedo
!
!-------------------------------------------------------------------------------
!
!*      0.2    Local variables
!              ---------------
!
!-------------------------------------------------------------------------------
!
!
!* albedo
PALB =  0.15
!* emissivity
PEMIS=  0.98
!* surface temperature set to air temperature
PTS  = PTA
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE GREENROOF_PROPERTIES
