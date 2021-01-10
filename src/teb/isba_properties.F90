!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE ISBA_PROPERTIES(IO, PEK, PDIR_SW, PSCA_SW, PSW_BANDS, KSW,  &
                                 PASNOW, PANOSNOW, PESNOW, PENOSNOW,         &
                                 PTSSNOW, PTSNOSNOW,                         &
                                 PALBNIR_TVEG, PALBVIS_TVEG, PALBNIR_TSOIL, PALBVIS_TSOIL )  
!     ##########################################################################
!
!!****  *ISBA_PROPERTIES*  
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
!!	S. Belair           * Meteo-France *
!!
!!    MODIFICATIONS
!!    ------------- 
!!      
!!      P. Samuelsson  02/2012  MEB
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_PE_t
!
USE MODD_TYPE_SNOW
USE MODD_SNOW_PAR   , ONLY : XEMISSN, XEMCRIN, XSNOWDMIN, &
                               XRHOSMAX_ES, XRHOSMIN_ES  
USE MODD_WATER_PAR  , ONLY : XEMISWAT
!
USE MODI_ISBA_SNOW_FRAC
USE MODI_ISBA_ALBEDO
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_PE_t), INTENT(INOUT) :: PEK
!
REAL, DIMENSION(:,:), INTENT(IN)   :: PDIR_SW            ! direct incoming solar radiation
REAL, DIMENSION(:,:), INTENT(IN)   :: PSCA_SW            ! diffus incoming solar radiation
REAL, DIMENSION(:)  , INTENT(IN)   :: PSW_BANDS          ! mean wavelength of each shortwave band (m)
INTEGER,              INTENT(IN)   :: KSW                ! number of short-wave spectral bands            
!
REAL, DIMENSION(:)  , INTENT(OUT)  :: PASNOW    ! = snow albedo
REAL, DIMENSION(:)  , INTENT(OUT)  :: PANOSNOW  ! = snow free albedo 
REAL, DIMENSION(:)  , INTENT(OUT)  :: PESNOW    ! = snow emissivity
REAL, DIMENSION(:)  , INTENT(OUT)  :: PENOSNOW  ! = snow free emissivity
REAL, DIMENSION(:)  , INTENT(OUT)  :: PTSSNOW   ! = snow radiative temperature
REAL, DIMENSION(:)  , INTENT(OUT)  :: PTSNOSNOW ! = snow free radiative temperature
!
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBNIR_TVEG       ! nearIR  veg tot albedo
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBVIS_TVEG       ! visible veg tot albedo
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBNIR_TSOIL      ! nearIR  soil tot albedo
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBVIS_TSOIL      ! visible soil tot albedo
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PDIR_SW,1)) :: ZGLOBAL_SW                 ! global incoming SW rad.
REAL, DIMENSION(SIZE(PEK%XALBNIR))   :: ZALBF
REAL, DIMENSION(SIZE(PEK%XALBNIR))   :: ZFFV
REAL, DIMENSION(SIZE(PEK%XALBNIR))   :: ZFFG
!
LOGICAL, PARAMETER :: GMEB=.FALSE.
REAL, DIMENSION(SIZE(PDIR_SW,1))   :: ZP_MEB_SCA_SW, ZALBNIR_TSNOW, ZALBVIS_TSNOW
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_PROPERTIES',0,ZHOOK_HANDLE)
!
 CALL ISBA_SNOW_FRAC(PEK%TSNOW%SCHEME, PEK%TSNOW%WSNOW, PEK%TSNOW%RHO, PEK%TSNOW%ALB, &
                     PEK%XVEG, PEK%XLAI, PEK%XZ0, &
                     PEK%XPSN, PEK%XPSNV_A, PEK%XPSNG, PEK%XPSNV )  
!
!-------------------------------------------------------------------------------
!*      2.     Compute snow-free albedo
!              ------------------------
!
!* Snow-free surface albedo for each wavelength
!
ZALBF         = 0.
ZFFV          = 0.
ZFFG          = 0.
!
 CALL ISBA_ALBEDO(PEK, IO%LTR_ML, GMEB, PDIR_SW, PSCA_SW, PSW_BANDS, KSW, &
                  ZALBF, ZFFV, ZFFG, ZGLOBAL_SW, ZP_MEB_SCA_SW,           &
                  PALBNIR_TVEG, PALBVIS_TVEG, PALBNIR_TSOIL, PALBVIS_TSOIL)

PANOSNOW(:) = PEK%XSNOWFREE_ALB(:)
!-------------------------------------------------------------------------------
!
!*      3.     Compute aggeragted albedo and emissivity
!              ----------------------------------------
!
IF(PEK%TSNOW%SCHEME == '3-L' .OR. PEK%TSNOW%SCHEME == 'CRO' .OR. IO%CISBA == 'DIF')THEN
!
! NON-SNOW covered Grid averaged albedo and emissivity for explicit snow scheme:
!
  PASNOW  (:) = PEK%TSNOW%ALB(:)
  PESNOW  (:) = PEK%TSNOW%EMIS(:)
  PENOSNOW(:) = PEK%XEMIS(:)

  PTSSNOW(:)   = PEK%TSNOW%TS(:)
  PTSNOSNOW(:) = PEK%XTG(:,1)

ELSE
!
! Grid averaged albedo and emissivity for composite snow scheme:
!
  IF(PEK%TSNOW%SCHEME =='EBA') THEN
!
    PASNOW  (:) = PEK%TSNOW%ALB(:)
    PESNOW  (:) = XEMCRIN
    PENOSNOW(:) = PEK%XEMIS(:)

    PTSSNOW  (:) = PEK%XTG(:,1)
    PTSNOSNOW(:) = PEK%XTG(:,1)

  ELSE

    PASNOW  (:) = PEK%TSNOW%ALB(:)
    PESNOW  (:) = XEMISSN
    PENOSNOW(:) = PEK%XEMIS(:)

    PTSSNOW  (:) = PEK%XTG(:,1)
    PTSNOSNOW(:) = PEK%XTG(:,1)

  ENDIF
!
ENDIF
IF (LHOOK) CALL DR_HOOK('ISBA_PROPERTIES',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ISBA_PROPERTIES
