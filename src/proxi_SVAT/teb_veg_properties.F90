!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE TEB_VEG_PROPERTIES (PMASK, IO, PEK, PDIR_SW, PSCA_SW, PSW_BANDS, KSW, &
                                     PTS, PEMIS, PALB, PTA, PALBNIR_TVEG, PALBVIS_TVEG,&
                                     PALBNIR_TSOIL, PALBVIS_TSOIL      )  
!     ##########################################################################
!
!!****  *GARDEN_PROPERTIES*  
!!
!!    PURPOSE
!!    -------
!
!     Calculates grid-averaged albedo and emissivity (according to snow scheme)
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
!!      S. Belair           * Meteo-France *
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_PE_t
!
USE MODD_SURF_PAR, ONLY : XUNDEF
!
USE MODI_ISBA_PROPERTIES
USE MODI_FLAG_TEB_VEG_n

USE PROXI_SVAT_CONSTS, ONLY: VEG_ALB, VEG_EMIS
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
REAL, DIMENSION(:), INTENT(IN) :: PMASK
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_PE_t), INTENT(INOUT) :: PEK
!
REAL, DIMENSION(:,:), INTENT(IN)   :: PDIR_SW            ! direct incoming solar radiation
REAL, DIMENSION(:,:), INTENT(IN)   :: PSCA_SW            ! diffus incoming solar radiation
REAL, DIMENSION(:)  , INTENT(IN)   :: PSW_BANDS          ! mean wavelength of each shortwave band (m)
INTEGER,              INTENT(IN)   :: KSW                ! number of short-wave spectral bands
!
REAL, DIMENSION(:)  , INTENT(OUT)  :: PTS                ! radiative surface temperature
REAL, DIMENSION(:)  , INTENT(OUT)  :: PEMIS              ! green areas emissivity
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALB               ! green areas albedo
!
REAL, DIMENSION(:)  , INTENT(IN), OPTIONAL :: PTA        ! Air temperature (K)
!
REAL, DIMENSION(:)  , INTENT(OUT), OPTIONAL  :: PALBNIR_TVEG       ! nearIR  veg tot albedo
REAL, DIMENSION(:)  , INTENT(OUT), OPTIONAL  :: PALBVIS_TVEG       ! visible veg tot albedo
REAL, DIMENSION(:)  , INTENT(OUT), OPTIONAL  :: PALBNIR_TSOIL      ! nearIR  soil tot albedo
REAL, DIMENSION(:)  , INTENT(OUT), OPTIONAL  :: PALBVIS_TSOIL      ! visible soil tot albedo
!-------------------------------------------------------------------------------
!
!*      0.2    Local variables
!              ---------------
!
!-------------------------------------------------------------------------------
!
!* albedo
PALB = VEG_ALB
!* emissivity
PEMIS=  VEG_EMIS
!* surface temperature set to air temperature
PTS  = PTA
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE TEB_VEG_PROPERTIES

