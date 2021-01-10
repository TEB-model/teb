!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!######################
MODULE MODD_DIAG_n
!######################
!
!!****  *MODD_DIAG - declaration of diagnostics for ISBA scheme
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004
!!      Modified    01/2006 : sea flux parameterization.
!!      Modified    04/2009 : precip for/from restart file.
!!      Modified    08/2009 : BUDGETC for all tiles
!!      P. Samuelsson 10/2014 : added min max for XT2M
!!      Modified    09/2015 : M Lafaysse LSNOWDIMNC
!!      Modified    06/2016 : M Lafaysse LRESETCUMUL
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_DATE_SURF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
TYPE DIAG_OPTIONS_t
!
  REAL    :: XDIAG_TSTEP  ! time step for diagnostics writing
!
  INTEGER :: N2M          ! flag for 2 meters (and 10 meters) quantities
  LOGICAL :: LT2MMW       ! flag to perform modified weighting of 2m temperature  
  LOGICAL :: L2M_MIN_ZS   ! flag for 2 meters quantities evaluated on
!                         ! the minimum orographyy of the grid      
  LOGICAL :: LSURF_BUDGET   ! flag for surface energy budget
  LOGICAL :: LRAD_BUDGET    ! flag for radiative energy budget
!
  LOGICAL :: LCOEF        ! flag for transfer coefficients
  LOGICAL :: LSURF_VARS   ! flag for surface variables
  LOGICAL :: LFRAC        ! flag for writing fractions of each four tiles
  LOGICAL :: LDIAG_GRID   ! flag for mean grid diag
!  
  LOGICAL :: LSURF_BUDGETC       ! flag for surface cumulated energy budget
  LOGICAL :: LRESET_BUDGETC      ! flag for surface cumulated energy budget
  LOGICAL :: LREAD_BUDGETC       ! flag for surface cumulated energy budget
  LOGICAL :: LPROVAR_TO_DIAG     ! switch to write (or not) prognostic variable
                                 ! and allows puting field in diagnostics 
  LOGICAL    :: LSNOWDIMNC       ! if true create a snow layer dimension in nc files
  LOGICAL    :: LRESETCUMUL      ! reset cumulated variables at 0 at each output timestep 
  LOGICAL    :: LSELECT          ! switch to control which fields are written
                                 ! (only those whose naem appears in in text array)  
!
  TYPE(DATE_TIME):: TIME_BUDGETC
!                                  
  CHARACTER(LEN=12), POINTER, DIMENSION(:) :: CSELECT  ! Name of ouput fields if LSELECT=true
!
  LOGICAL :: LPGD          ! flag for writing of PGD files
  LOGICAL :: LPATCH_BUDGET ! flag for patch output
!
END TYPE DIAG_OPTIONS_t
!
TYPE DIAG_t
!------------------------------------------------------------------------------
!
!* variables for each patch
!
!* averaged variables
!
  REAL, POINTER, DIMENSION(:)   :: XRI       ! Bulk-Richardson number           (-)
  REAL, POINTER, DIMENSION(:)   :: XCD       ! drag coefficient for wind        (W/s2)
  REAL, POINTER, DIMENSION(:)   :: XCDN       ! neutral drag coefficient                      (-)  
  REAL, POINTER, DIMENSION(:)   :: XCH       ! drag coefficient for heat        (W/s)
  REAL, POINTER, DIMENSION(:)   :: XCE       ! drag coefficient for vapor       (W/s/K)
!
  REAL, POINTER, DIMENSION(:)   :: XHU        ! area averaged surface humidity coefficient    (-)
  REAL, POINTER, DIMENSION(:)   :: XHUG       ! baresoil surface humidity coefficient         (-)
  REAL, POINTER, DIMENSION(:)   :: XHV        ! Halstead coefficient                          (-)  
!
  REAL, POINTER, DIMENSION(:)   :: XRN       ! net radiation at surface         (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XH        ! sensible heat flux               (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLE       ! total latent heat flux           (W/m2) 
  REAL, POINTER, DIMENSION(:)   :: XLEI      ! sublimation latent heat flux     (W/m2) 
  REAL, POINTER, DIMENSION(:)   :: XGFLUX    ! net soil-vegetation flux         (W/m2)
!
  REAL, POINTER, DIMENSION(:)   :: XEVAP    ! total evaporation                (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XSUBL    ! sublimation                      (kg/m2/s)
!
  REAL, POINTER, DIMENSION(:)   :: XTS       ! surface temperature              (K)
  REAL, POINTER, DIMENSION(:)   :: XTSRAD    ! surface temperature              (K)
  REAL, POINTER, DIMENSION(:)   :: XALBT          ! Total Albedo  
  REAL, POINTER, DIMENSION(:)   :: XSWE     ! snow water equivalent (kg/m2)
!  
  REAL, POINTER, DIMENSION(:)   :: XT2M      ! temperature at 2 meters          (K)
  REAL, POINTER, DIMENSION(:)   :: XT2M_MIN  ! Minimum temperature at 2 meters          (K)
  REAL, POINTER, DIMENSION(:)   :: XT2M_MAX  ! Maximum temperature at 2 meters          (K)
  REAL, POINTER, DIMENSION(:)   :: XQ2M      ! humidity    at 2 meters          (kg/kg)
  REAL, POINTER, DIMENSION(:)   :: XHU2M     ! relative humidity at 2 meters    (-)
  REAL, POINTER, DIMENSION(:)   :: XHU2M_MIN ! Minimum relative humidity at 2 meters    (-)
  REAL, POINTER, DIMENSION(:)   :: XHU2M_MAX ! Maximum relative humidity at 2 meters    (-)
  REAL, POINTER, DIMENSION(:)   :: XQS       ! humidity at surface              (kg/kg)
  REAL, POINTER, DIMENSION(:)   :: XZON10M   ! zonal wind at 10 meters          (m/s)
  REAL, POINTER, DIMENSION(:)   :: XMER10M   ! meridian wind at 10 meters       (m/s)
  REAL, POINTER, DIMENSION(:)   :: XWIND10M  ! wind at 10 meters                (m/s)
  REAL, POINTER, DIMENSION(:)   :: XWIND10M_MAX  ! Maximum wind at 10 meters    (m/s)
!
  REAL, POINTER, DIMENSION(:)   :: XSFCO2    ! CO2 flux                         (m/s*kg_CO2/kg_air)
!
  REAL, POINTER, DIMENSION(:,:) :: XSWBD     ! downward short wave radiation by spectral band   (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XSWBU     ! upward short wave radiation by spectral band (W/m2)
!
  REAL, POINTER, DIMENSION(:)   :: XLWD      ! downward long wave radiation     (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLWU      ! upward long wave radiation       (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XSWD      ! downward short wave radiation    (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XSWU      ! upward short wave radiation      (W/m2)
!
  REAL, POINTER, DIMENSION(:)   :: XFMU      ! horizontal momentum flux zonal   (m2/s2)
  REAL, POINTER, DIMENSION(:)   :: XFMV      ! horizontal momentum flux meridian (m2/s2) 
  !                                                
  REAL, POINTER, DIMENSION(:)   :: XZ0       ! roughness length for momentum
                                                 ! for vegetation and snow    (m)
  REAL, POINTER, DIMENSION(:)   :: XZ0H      ! roughness length for heat
                                                 ! for vegetation and snow    (m)
  REAL, POINTER, DIMENSION(:)   :: XZ0EFF    ! effective roughness length for heat
                                                 ! for vegetation and snow    (m)
!
  REAL, POINTER, DIMENSION(:)   :: XT2M_MIN_ZS ! air temperature at 2 meters   (K)
  REAL, POINTER, DIMENSION(:)   :: XQ2M_MIN_ZS ! air humidity at 2 meters      (kg/kg)
  REAL, POINTER, DIMENSION(:)   :: XHU2M_MIN_ZS! air relative humidity at 2 m  (-)



  REAL, POINTER, DIMENSION(:)   :: XPS          ! air pressure at the surface      (Pa)
  REAL, POINTER, DIMENSION(:)   :: XRHOA        ! air density  at the surface      (kg/m3)

  REAL, POINTER, DIMENSION(:)   :: XSSO_FMU     ! zonal friction    (with SSO)     (Pa)
  REAL, POINTER, DIMENSION(:)   :: XSSO_FMV     ! meridian friction (with SSO)     (Pa)
!

  REAL, POINTER, DIMENSION(:)   :: XUREF   ! reference height for momentum    (m)
  REAL, POINTER, DIMENSION(:)   :: XZREF   ! reference height for heat        (m)
  REAL, POINTER, DIMENSION(:)   :: XTRAD   ! radiative temperature at t       (K)
  REAL, POINTER, DIMENSION(:)   :: XEMIS   ! surface emissivity at t          (-)

!------------------------------------------------------------------------------
!
END TYPE DIAG_t
!
TYPE DIAG_NP_t
!
TYPE(DIAG_t), POINTER :: AL(:)=>NULL() 
!
END TYPE DIAG_NP_t
!
CONTAINS
!
SUBROUTINE DIAG_OPTIONS_INIT(DGO)
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DGO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_N:DIAG_OPTIONS_INIT",0,ZHOOK_HANDLE)
  NULLIFY(DGO%CSELECT)
DGO%XDIAG_TSTEP=0.
DGO%N2M=0
DGO%LT2MMW=.FALSE.
DGO%L2M_MIN_ZS=.FALSE.
DGO%LSURF_BUDGET=.FALSE.
DGO%LRAD_BUDGET=.FALSE.
DGO%LCOEF=.FALSE.
DGO%LSURF_VARS=.FALSE.
DGO%LFRAC=.FALSE.
DGO%LDIAG_GRID=.FALSE.
DGO%LPGD=.FALSE.
DGO%LPATCH_BUDGET=.FALSE.
DGO%LSURF_BUDGETC=.FALSE.
DGO%LRESET_BUDGETC=.FALSE.
DGO%LREAD_BUDGETC=.FALSE.
DGO%LPROVAR_TO_DIAG=.FALSE.
DGO%LSNOWDIMNC=.FALSE.
DGO%LRESETCUMUL=.FALSE.
DGO%LSELECT=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_N:DIAG_OPTIONS_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_OPTIONS_INIT
!
SUBROUTINE DIAG_NP_INIT(ND,KPATCH)
TYPE(DIAG_NP_t), INTENT(INOUT) :: ND 
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_N:DIAG_NP_INIT",0,ZHOOK_HANDLE)
IF (.NOT.ASSOCIATED(ND%AL)) THEN
  ALLOCATE(ND%AL(KPATCH))
  DO JP=1,KPATCH
    CALL DIAG_INIT(ND%AL(JP))
  ENDDO
ELSE
  DO JP=1,KPATCH
    CALL DIAG_INIT(ND%AL(JP))
  ENDDO
  DEALLOCATE(ND%AL)
ENDIF
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_N:DIAG_NP_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_NP_INIT
!
!
SUBROUTINE DIAG_INIT(D)
TYPE(DIAG_t), INTENT(INOUT) :: D
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_N:DIAG_INIT",0,ZHOOK_HANDLE)
  NULLIFY(D%XRI)
  NULLIFY(D%XCD)
  NULLIFY(D%XCDN)
  NULLIFY(D%XCH)
  NULLIFY(D%XCE)
  NULLIFY(D%XHU)
  NULLIFY(D%XHUG)
  NULLIFY(D%XHV)  
  NULLIFY(D%XRN)
  NULLIFY(D%XH)
  NULLIFY(D%XLE)
  NULLIFY(D%XLEI)
  NULLIFY(D%XGFLUX)
  NULLIFY(D%XEVAP)
  NULLIFY(D%XSUBL)  
  NULLIFY(D%XTS)
  NULLIFY(D%XTSRAD)
  NULLIFY(D%XALBT)  
  NULLIFY(D%XSWE)  
  NULLIFY(D%XT2M)
  NULLIFY(D%XT2M_MIN)
  NULLIFY(D%XT2M_MAX)
  NULLIFY(D%XQ2M)
  NULLIFY(D%XHU2M)
  NULLIFY(D%XQS)
  NULLIFY(D%XZON10M)
  NULLIFY(D%XMER10M)
  NULLIFY(D%XWIND10M)
  NULLIFY(D%XWIND10M_MAX)  
  NULLIFY(D%XLWD)
  NULLIFY(D%XLWU)
  NULLIFY(D%XSWD)
  NULLIFY(D%XSWU)
  NULLIFY(D%XSWBD)
  NULLIFY(D%XSWBU)
  NULLIFY(D%XFMU)
  NULLIFY(D%XFMV)
  NULLIFY(D%XZ0)
  NULLIFY(D%XZ0H)
  NULLIFY(D%XZ0EFF)
  NULLIFY(D%XT2M_MIN_ZS)
  NULLIFY(D%XQ2M_MIN_ZS)
  NULLIFY(D%XHU2M_MIN_ZS)
  NULLIFY(D%XPS)
  NULLIFY(D%XRHOA)
  NULLIFY(D%XSSO_FMU)
  NULLIFY(D%XSSO_FMV)
  NULLIFY(D%XUREF)
  NULLIFY(D%XZREF)
  NULLIFY(D%XTRAD)
  NULLIFY(D%XEMIS)
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_N:DIAG_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_INIT


END MODULE MODD_DIAG_n
