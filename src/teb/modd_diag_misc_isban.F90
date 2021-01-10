!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ############################
      MODULE MODD_DIAG_MISC_ISBA_n
!     ############################
!
!!****  *MODD_DIAG_MISC_ISBA - declaration of packed surface parameters for ISBA scheme
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
!!      P. Le Moigne   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       07/10/04
!!      A.L. Gibelin 04/2009 : Add respiration diagnostics
!!      A.L. Gibelin 05/2009 : Add carbon spinup
!!      A.L. Gibelin 07/2009 : Suppress RDK and transform GPP as a diagnostic
!!       B. Decharme 05/2012 : Carbon fluxes in diag_evap
!!       B. Decharme 07/2012 : New diag for DIF under LSURF_MISC_DIF key
!!                               F2 stress
!!                               Root zone swi, wg and wgi
!!                               swi, wg and wgi comparable to ISBA-FR-DG2 and DG3 layers
!!       B. Decharme 10/2012 : New diag for DIF 
!!                               active layer thickness over permafrost area
!!                               frozen layer thickness over non-permafrost area
!!       M. Lafaysse 09/2015 : new Crocus-MEPRA outputs
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE DIAG_MISC_ISBA_t
!------------------------------------------------------------------------------
!
  LOGICAL :: LSURF_MISC_BUDGET   ! flag for miscellaneous terms of isba scheme
  LOGICAL :: LSURF_DIAG_ALBEDO   ! flag to write out diagnostic albedo
  LOGICAL :: LSURF_MISC_DIF      ! flag for miscellaneous terms of isba-dif scheme
  LOGICAL :: LPROSNOW            ! flag for Crocus-MEPRA outputs
  LOGICAL :: LVOLUMETRIC_SNOWLIQ ! volumetric snow liquid water content (kg m-3)  
!
!* variables for each patch
!
  REAL, POINTER, DIMENSION(:) :: XHV       ! Halstead coefficient
  REAL, POINTER, DIMENSION(:) :: XLAI      ! leaf average index  
!      
  REAL, POINTER, DIMENSION(:,:) :: XSWI        ! Soil wetness index
  REAL, POINTER, DIMENSION(:,:) :: XTSWI       ! Total soil wetness index
  REAL, POINTER, DIMENSION(:)   :: XSOIL_SWI     ! Soil wetness index
  REAL, POINTER, DIMENSION(:)   :: XSOIL_TSWI    ! Total Soil wetness index
  REAL, POINTER, DIMENSION(:)   :: XSOIL_TWG     ! Soil water content (liquid+ice) (kg.m-2)
  REAL, POINTER, DIMENSION(:)   :: XSOIL_TWGI    ! Soil ice content (kg.m-2)
  REAL, POINTER, DIMENSION(:)   :: XSOIL_WG     ! Soil water content (liquid+ice) (m3.m-3)
  REAL, POINTER, DIMENSION(:)   :: XSOIL_WGI    ! Soil ice content (m3.m-3)  
!     
  REAL, POINTER, DIMENSION(:) :: XFRD2_TSWI      ! ISBA-FR-DG2 comparable soil wetness index (DIF option)
  REAL, POINTER, DIMENSION(:) :: XFRD2_TWG       ! ISBA-FR-DG2 comparable soil water content (liquid+ice) (DIF option)
  REAL, POINTER, DIMENSION(:) :: XFRD2_TWGI      ! ISBA-FR-DG2 comparable soil ice content (DIF option)  
  REAL, POINTER, DIMENSION(:) :: XFRD3_TSWI      ! ISBA-FR-Deep comparable soil wetness index (DIF option)
  REAL, POINTER, DIMENSION(:) :: XFRD3_TWG       ! ISBA-FR-Deep comparable soil water content (liquid+ice) (DIF option)
  REAL, POINTER, DIMENSION(:) :: XFRD3_TWGI      ! ISBA-FR-Deep comparable soil ice content (DIF option)
!
  REAL, POINTER, DIMENSION(:)   :: XALT        ! Active layer thickness in permafrost area
  REAL, POINTER, DIMENSION(:)   :: XFLT        ! Frozen layer thickness in non-permmafrost area
!
  REAL, POINTER, DIMENSION(:) :: XRNSNOW    ! net radiative flux from snow (ISBA-ES:3-L)    (W/m2)
  REAL, POINTER, DIMENSION(:) :: XHSNOW     ! sensible heat flux from snow (ISBA-ES:3-L)    (W/m2)
  REAL, POINTER, DIMENSION(:) :: XHPSNOW    ! heat release from rainfall (ISBA-ES:3-L)      (W/m2)
  REAL, POINTER, DIMENSION(:) :: XGFLUXSNOW ! net surface energy flux into snowpack      
!                                               ! (ISBA-ES:3-L)                                 (W/m2)
  REAL, POINTER, DIMENSION(:) :: XUSTARSNOW ! friction velocity  over snow 
!                                               ! (ISBA-ES:3-L)                                 (m/s)
  REAL, POINTER, DIMENSION(:) :: XGRNDFLUX  ! soil/snow interface heat flux (ISBA-ES:3-L)   (W/m2)
  REAL, POINTER, DIMENSION(:) :: XSRSFC     ! snowfall over snowpack (ISBA-ES:3-L)          (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XRRSFC     ! rainfall over snowpack (ISBA-ES:3-L)          (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XCDSNOW    ! snow drag coefficient (ISBA-ES:3-L)           (-)
  REAL, POINTER, DIMENSION(:) :: XCHSNOW    ! heat turbulent transfer coefficient 
!                                               ! (ISBA-ES:3-L)                                 (-)
  REAL, POINTER, DIMENSION(:,:) :: XSNOWDZ     ! snow layer thicknesses                        (m)
  REAL, POINTER, DIMENSION(:,:) :: XSNOWDEND   ! dendricity (Crocus)
  REAL, POINTER, DIMENSION(:,:) :: XSNOWSPHER  ! sphericity (Crocus)
  REAL, POINTER, DIMENSION(:,:) :: XSNOWSIZE   ! grain size (Crocus)
  REAL, POINTER, DIMENSION(:,:) :: XSNOWSSA    ! snow specific surface area (Crocus)
  REAL, POINTER, DIMENSION(:,:) :: XSNOWTYPEMEPRA ! snow grain type (Crocus-MEPRA)
  REAL, POINTER, DIMENSION(:,:) :: XSNOWRAM    ! snow ram resistance (Crocus-MEPRA)
  REAL, POINTER, DIMENSION(:,:) :: XSNOWSHEAR  ! snow shear resistance (Crocus-MEPRA)  
!
  REAL, POINTER, DIMENSION(:) :: XSNOWHMASS ! heat content change due to mass
!                                           ! changes in snowpack: for budget
  REAL, POINTER, DIMENSION(:,:) :: XSNOWLIQ    ! snow liquid water profile (ISBA-ES:3-L)
  REAL, POINTER, DIMENSION(:,:) :: XSNOWTEMP   ! snow temperature profile  (ISBA-ES:3-L)
!     
  REAL, POINTER, DIMENSION(:) :: XTWSNOW       ! Total snow reservoir
  REAL, POINTER, DIMENSION(:) :: XTDSNOW       ! Total snow height
  REAL, POINTER, DIMENSION(:) :: XTTSNOW       ! Total snow temperature
!
  REAL, POINTER, DIMENSION(:) :: XSNDPT_1DY, XSNDPT_3DY, XSNDPT_5DY, XSNDPT_7DY ! fresh snow depth in 1, 3, 5, 7 days
  REAL, POINTER, DIMENSION(:) :: XSNSWE_1DY, XSNSWE_3DY, XSNSWE_5DY, XSNSWE_7DY! fresh snow water equivalent in 1, 3, 5, 7 days
  REAL, POINTER, DIMENSION(:) :: XSNRAM_SONDE ! penetration depth of the ram resistance sensor (2 DaN)
  REAL, POINTER, DIMENSION(:) :: XSN_WETTHCKN ! Thickness of wet snow at the top of the snowpack
  REAL, POINTER, DIMENSION(:) :: XSN_REFRZNTHCKN  ! Thickness of refrozen snow at the top of the snowpack

  REAL, POINTER, DIMENSION(:) :: XPSNG         ! Snow fraction over ground, diag at time t
  REAL, POINTER, DIMENSION(:) :: XPSNV         ! Snow fraction over vegetation, diag at time t
  REAL, POINTER, DIMENSION(:) :: XPSN          ! Total Snow fraction, diag at time t
!    
  REAL, POINTER, DIMENSION(:) :: XFSAT         ! Topmodel/dt92 saturated fraction
!
  REAL, POINTER, DIMENSION(:) :: XFFG          ! Flood fraction over ground, diag at time t
  REAL, POINTER, DIMENSION(:) :: XFFV          ! Flood fraction over vegetation, diag at time t
  REAL, POINTER, DIMENSION(:) :: XFF           ! Total Flood fraction, diag at time t
!
  REAL, POINTER, DIMENSION(:) :: XSEUIL        ! Irrigation threshold
!
  REAL, POINTER, DIMENSION(:) :: XFAPAR        ! Fapar of vegetation
  REAL, POINTER, DIMENSION(:) :: XFAPIR        ! Fapir of vegetation
  REAL, POINTER, DIMENSION(:) :: XDFAPARC      ! Fapar of vegetation (cumul)
  REAL, POINTER, DIMENSION(:) :: XDFAPIRC      ! Fapir of vegetation (cumul)
  REAL, POINTER, DIMENSION(:) :: XFAPAR_BS     ! Fapar of bare soil
  REAL, POINTER, DIMENSION(:) :: XFAPIR_BS     ! Fapir of bare soil
  REAL, POINTER, DIMENSION(:) :: XDLAI_EFFC    ! Effective LAI (cumul)
!
  REAL, POINTER, DIMENSION(:) :: XCG        ! heat capacity of the ground
  REAL, POINTER, DIMENSION(:) :: XC1        ! coefficients for the moisure
  REAL, POINTER, DIMENSION(:) :: XC2        ! equation.
  REAL, POINTER, DIMENSION(:) :: XWGEQ      ! equilibrium volumetric water
!                                               ! content
  REAL, POINTER, DIMENSION(:) :: XCT        ! area-averaged heat capacity
  REAL, POINTER, DIMENSION(:) :: XRS        ! stomatal resistance                            (s/m)
!
!------------------------------------------------------------------------------
!
END TYPE DIAG_MISC_ISBA_t
!
TYPE DIAG_MISC_ISBA_NP_t
!
TYPE(DIAG_MISC_ISBA_t), POINTER :: AL(:)=>NULL() 
!
END TYPE DIAG_MISC_ISBA_NP_t
!
CONTAINS
!
SUBROUTINE DIAG_MISC_ISBA_NP_INIT(YNDIAG_MISC_ISBA,KPATCH)
TYPE(DIAG_MISC_ISBA_NP_t), INTENT(INOUT) :: YNDIAG_MISC_ISBA 
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_ISBA_N:DIAG_MISC_ISBA_NP_INIT",0,ZHOOK_HANDLE)
IF (.NOT.ASSOCIATED(YNDIAG_MISC_ISBA%AL)) THEN
  ALLOCATE(YNDIAG_MISC_ISBA%AL(KPATCH))
  DO JP=1,KPATCH
    CALL DIAG_MISC_ISBA_INIT(YNDIAG_MISC_ISBA%AL(JP))
  ENDDO
ELSE
  DO JP=1,KPATCH
    CALL DIAG_MISC_ISBA_INIT(YNDIAG_MISC_ISBA%AL(JP))
  ENDDO
  DEALLOCATE(YNDIAG_MISC_ISBA%AL)        
ENDIF
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_ISBA_N:DIAG_MISC_ISBA_NP_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_MISC_ISBA_NP_INIT
!
SUBROUTINE DIAG_MISC_ISBA_INIT(DMI)
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DMI
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_ISBA_N:DIAG_MISC_ISBA_INIT",0,ZHOOK_HANDLE)
  NULLIFY(DMI%XHV)
  NULLIFY(DMI%XLAI)  
  NULLIFY(DMI%XSWI)
  NULLIFY(DMI%XTSWI)
  NULLIFY(DMI%XSOIL_SWI)
  NULLIFY(DMI%XSOIL_TSWI)
  NULLIFY(DMI%XSOIL_TWG)
  NULLIFY(DMI%XSOIL_TWGI)
  NULLIFY(DMI%XSOIL_WG)
  NULLIFY(DMI%XSOIL_WGI)
  NULLIFY(DMI%XFRD2_TWG)
  NULLIFY(DMI%XFRD2_TWGI)
  NULLIFY(DMI%XFRD3_TSWI)
  NULLIFY(DMI%XFRD3_TWG)
  NULLIFY(DMI%XFRD3_TWGI)    
  NULLIFY(DMI%XALT)
  NULLIFY(DMI%XFLT)
  NULLIFY(DMI%XRNSNOW)
  NULLIFY(DMI%XHSNOW)
  NULLIFY(DMI%XHPSNOW)
  NULLIFY(DMI%XGFLUXSNOW)
  NULLIFY(DMI%XUSTARSNOW)
  NULLIFY(DMI%XGRNDFLUX)
  NULLIFY(DMI%XSRSFC)
  NULLIFY(DMI%XRRSFC)
  NULLIFY(DMI%XCDSNOW)
  NULLIFY(DMI%XCHSNOW)
  NULLIFY(DMI%XSNOWDZ)
  NULLIFY(DMI%XSNOWDEND)
  NULLIFY(DMI%XSNOWSPHER)
  NULLIFY(DMI%XSNOWSIZE)
  NULLIFY(DMI%XSNOWTYPEMEPRA)
  NULLIFY(DMI%XSNOWSSA)
  NULLIFY(DMI%XSNOWRAM)
  NULLIFY(DMI%XSNOWSHEAR)
  NULLIFY(DMI%XSNDPT_1DY)
  NULLIFY(DMI%XSNDPT_3DY)
  NULLIFY(DMI%XSNDPT_5DY)
  NULLIFY(DMI%XSNDPT_7DY) 
  NULLIFY(DMI%XSNSWE_1DY)
  NULLIFY(DMI%XSNSWE_3DY)
  NULLIFY(DMI%XSNSWE_5DY)
  NULLIFY(DMI%XSNSWE_7DY)
  NULLIFY(DMI%XSNRAM_SONDE)
  NULLIFY(DMI%XSN_REFRZNTHCKN)
  NULLIFY(DMI%XSN_WETTHCKN)  
  NULLIFY(DMI%XSNOWHMASS)  
  NULLIFY(DMI%XSNOWLIQ)
  NULLIFY(DMI%XSNOWTEMP)
  NULLIFY(DMI%XTWSNOW)
  NULLIFY(DMI%XTDSNOW)
  NULLIFY(DMI%XTTSNOW)
  NULLIFY(DMI%XPSNG)
  NULLIFY(DMI%XPSNV)
  NULLIFY(DMI%XPSN)
  NULLIFY(DMI%XFFG)
  NULLIFY(DMI%XFFV)
  NULLIFY(DMI%XFF)
  NULLIFY(DMI%XSEUIL)
  NULLIFY(DMI%XFAPAR)
  NULLIFY(DMI%XFAPIR)
  NULLIFY(DMI%XDFAPARC)
  NULLIFY(DMI%XDFAPIRC)
  NULLIFY(DMI%XFAPAR_BS)
  NULLIFY(DMI%XFAPIR_BS)
  NULLIFY(DMI%XDLAI_EFFC) 
  NULLIFY(DMI%XCG)
  NULLIFY(DMI%XC1)
  NULLIFY(DMI%XC2)
  NULLIFY(DMI%XWGEQ)
  NULLIFY(DMI%XCT)
  NULLIFY(DMI%XRS)  
DMI%LSURF_MISC_BUDGET=.FALSE.
DMI%LSURF_DIAG_ALBEDO=.FALSE.
DMI%LSURF_MISC_DIF=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_ISBA_N:DIAG_MISC_ISBA_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_MISC_ISBA_INIT


END MODULE MODD_DIAG_MISC_ISBA_n
