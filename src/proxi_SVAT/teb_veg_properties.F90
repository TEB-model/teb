!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL licence version 2.1.
! See the following links for details:
! https://cecill.info/licences/Licence_CeCILL_V2.1-en.txt
! https://cecill.info/licences/Licence_CeCILL_V2.1-fr.txt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     #########
      SUBROUTINE TEB_VEG_PROPERTIES(PDIR_SW, PSCA_SW, PSW_BANDS, KSW,      &
                                   PTS, PEMIS, PALB, PTA                  )  
!     ##########################################################################
!
!!****  *TEB_VEG_PROPERTIES*  
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
!* albedo
PALB =  0.15
!* emissivity
PEMIS=  0.98
!* surface temperature set to air temperature
PTS  = PTA
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE TEB_VEG_PROPERTIES

