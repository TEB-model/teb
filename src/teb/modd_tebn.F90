!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ################
      MODULE MODD_TEB_n
!     ################
!
!!****  *MODD_TEB_n - declaration of surface parameters for urban surface
!!
!!    PURPOSE
!!    -------
!     Declaration of surface parameters
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
!!      A. Lemonsu      07/2012         Key for urban hydrology
!!      V. Masson       06/2013         splits module in two
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_SNOW
USE MODD_TYPE_DATE_SURF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!--------------------------------------------------------------------------

TYPE TEB_t
! TEB scheme option
!
! Geometric Parameters:
!
  REAL, POINTER, DIMENSION(:)   :: XROAD_DIR     ! Road direction (deg from North, clockwise)
  REAL, POINTER, DIMENSION(:)   :: XGARDEN       ! fraction of veg in the streets   (-)
  REAL, POINTER, DIMENSION(:)   :: XGREENROOF    ! fraction of greenroofs on roofs  (-)
  REAL, POINTER, DIMENSION(:)   :: XBLD          ! fraction of buildings            (-)
  REAL, POINTER, DIMENSION(:)   :: XROAD         ! fraction of roads                (-)
  REAL, POINTER, DIMENSION(:)   :: XCAN_HW_RATIO ! canyon    h/W                    (-)
  REAL, POINTER, DIMENSION(:)   :: XBLD_HEIGHT   ! buildings height 'h'             (m)
  REAL, POINTER, DIMENSION(:)   :: XWALL_O_HOR   ! wall surf. / hor. surf.          (-)
  REAL, POINTER, DIMENSION(:)   :: XROAD_O_GRND  ! road surf. / (road + garden surf.) (-)
  REAL, POINTER, DIMENSION(:)   :: XGARDEN_O_GRND! gard. surf. / (road + garden surf.)(-)
  REAL, POINTER, DIMENSION(:)   :: XWALL_O_GRND  ! wall surf. / (road + garden surf.) (-)
  REAL, POINTER, DIMENSION(:)   :: XWALL_O_BLD   ! wall surf. / bld surf. (-)
  REAL, POINTER, DIMENSION(:)   :: XZ0_TOWN      ! roughness length for momentum    (m)
  REAL, POINTER, DIMENSION(:)   :: XSVF_ROAD     ! road sky view factor             (-)
  REAL, POINTER, DIMENSION(:)   :: XSVF_GARDEN   ! green area sky view factor       (-)
  REAL, POINTER, DIMENSION(:)   :: XSVF_WALL     ! wall sky view factor             (-)
!
! Roof parameters
!
  REAL, POINTER, DIMENSION(:)   :: XALB_ROOF     ! roof albedo                      (-)
  REAL, POINTER, DIMENSION(:)   :: XEMIS_ROOF    ! roof emissivity                  (-)
  REAL, POINTER, DIMENSION(:,:) :: XHC_ROOF      ! roof layers heat capacity        (J/K/m3)
  REAL, POINTER, DIMENSION(:,:) :: XTC_ROOF      ! roof layers thermal conductivity (W/K/m)
  REAL, POINTER, DIMENSION(:,:) :: XD_ROOF       ! depth of roof layers             (m)
  REAL, POINTER, DIMENSION(:)   :: XROUGH_ROOF   ! roof roughness coef
!
!
! Road parameters
!
  REAL, POINTER, DIMENSION(:)   :: XALB_ROAD     ! road albedo                      (-)
  REAL, POINTER, DIMENSION(:)   :: XEMIS_ROAD    ! road emissivity                  (-)
  REAL, POINTER, DIMENSION(:,:) :: XHC_ROAD      ! road layers heat capacity        (J/K/m3)
  REAL, POINTER, DIMENSION(:,:) :: XTC_ROAD      ! road layers thermal conductivity (W/K/m)
  REAL, POINTER, DIMENSION(:,:) :: XD_ROAD       ! depth of road layers             (m)
!
! Wall parameters
!
  REAL, POINTER, DIMENSION(:)   :: XALB_WALL     ! wall albedo                      (-)
  REAL, POINTER, DIMENSION(:)   :: XEMIS_WALL    ! wall emissivity                  (-)
  REAL, POINTER, DIMENSION(:,:) :: XHC_WALL      ! wall layers heat capacity        (J/K/m3)
  REAL, POINTER, DIMENSION(:,:) :: XTC_WALL      ! wall layers thermal conductivity (W/K/m)
  REAL, POINTER, DIMENSION(:,:) :: XD_WALL       ! depth of wall layers             (m)
  REAL, POINTER, DIMENSION(:)   :: XROUGH_WALL   ! wall roughness coef
!
! Building's use type
!
  REAL, POINTER, DIMENSION(:)   :: XRESIDENTIAL  ! fraction of Residential use      (-)
  REAL                          :: XDT_RES       ! target temperature change when unoccupied (K) (residential buildings)
  REAL                          :: XDT_OFF       ! target temperature change when unoccupied (K) (offices buildings)
  
!
! anthropogenic fluxes
!
  REAL, POINTER, DIMENSION(:)   :: XH_TRAFFIC    ! anthropogenic sensible
!                                                  ! heat fluxes due to traffic       (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLE_TRAFFIC   ! anthropogenic latent
!                                                  ! heat fluxes due to traffic       (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XH_INDUSTRY   ! anthropogenic sensible                   
!                                                  ! heat fluxes due to factories     (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLE_INDUSTRY  ! anthropogenic latent
!                                                  ! heat fluxes due to factories     (W/m2)
!
! temperatures for boundary conditions
!
  REAL, POINTER, DIMENSION(:)   :: XTI_ROAD      ! road interior temperature        (K)
!
! Prognostic variables:
!
  REAL, POINTER, DIMENSION(:)   :: XWS_ROOF      ! roof water reservoir             (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XWS_ROAD      ! road water reservoir             (kg/m2)
  REAL, POINTER, DIMENSION(:,:) :: XT_ROOF       ! roof layer temperatures          (K)
  REAL, POINTER, DIMENSION(:,:) :: XT_ROAD       ! road layer temperatures          (K)
  REAL, POINTER, DIMENSION(:,:) :: XT_WALL_A     ! wall layer temperatures          (K)
  REAL, POINTER, DIMENSION(:,:) :: XT_WALL_B     ! wall layer temperatures          (K)
!
  REAL, POINTER, DIMENSION(:)   :: XAC_ROOF      ! roof aerodynamic conductance     ()
  REAL, POINTER, DIMENSION(:)   :: XAC_ROAD      ! road aerodynamic conductance     ()
  REAL, POINTER, DIMENSION(:)   :: XAC_WALL      ! wall aerodynamic conductance     ()
  REAL, POINTER, DIMENSION(:)   :: XAC_TOP       ! top  aerodynamic conductance     ()
  REAL, POINTER, DIMENSION(:)   :: XAC_ROOF_WAT  ! water aerodynamic conductance    ()
  REAL, POINTER, DIMENSION(:)   :: XAC_ROAD_WAT  ! water aerodynamic conductance    ()
!
  REAL, POINTER, DIMENSION(:)   :: XQSAT_ROOF    ! humidity of saturation for roofs (kg/kg)
  REAL, POINTER, DIMENSION(:)   :: XQSAT_ROAD    ! humidity of saturation for roads (kg/kg)
!
  REAL, POINTER, DIMENSION(:)   :: XDELT_ROOF    ! humidity of saturation for roofs (-)
  REAL, POINTER, DIMENSION(:)   :: XDELT_ROAD    ! humidity of saturation for roads (-)
!
! Semi-prognostic variables:
!
  REAL, POINTER, DIMENSION(:)   :: XT_CANYON     ! canyon air temperature           (K)
  REAL, POINTER, DIMENSION(:)   :: XQ_CANYON     ! canyon air specific humidity     (kg/kg)
!
!
! Prognostic snow:
!
  TYPE(SURF_SNOW)                 :: TSNOW_ROOF      ! snow state on roofs: 
!                                                  ! scheme type/option               (-)
!                                                  ! number of layers                 (-)
!                                                  ! snow (& liq. water) content      (kg/m2)
!                                                  ! heat content                     (J/m2)
!                                                  ! temperature                      (K)
!                                                  ! density                          (kg m-3)
  TYPE(SURF_SNOW)                 :: TSNOW_ROAD      ! snow state on roads: 
!                                                  ! scheme type/option               (-)
!                                                  ! number of layers                 (-)
!                                                  ! snow (& liq. water) content      (kg/m2)
!                                                  ! heat content                     (J/m2)
!                                                  ! temperature                      (K)
!                                                  ! density                          (kg m-3)
!                                                  ! density                          (kg m-3)
  TYPE(SURF_SNOW)                 :: TSNOW_GARDEN    ! snow state on green areas:
!                                                  ! scheme type/option               (-)
!                                                  ! number of layers                 (-)
!                                                  ! snow (& liq. water) content      (kg/m2)
!                                                  ! heat content                     (J/m2)
!                                                  ! temperature                      (K)
!                                                  ! density                          (kg m-3)
!
END TYPE TEB_t

TYPE TEB_NP_t
  !
  TYPE(TEB_t), POINTER :: AL(:) => NULL()
  !
END TYPE TEB_NP_t
!
CONTAINS
!----------------------------------------------------------------------------
!
SUBROUTINE TEB_INIT(YTEB)
TYPE(TEB_t), INTENT(INOUT) :: YTEB
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_N:TEB_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YTEB%XROAD_DIR)
  NULLIFY(YTEB%XGARDEN)
  NULLIFY(YTEB%XGREENROOF)
  NULLIFY(YTEB%XBLD)
  NULLIFY(YTEB%XROAD)
  NULLIFY(YTEB%XCAN_HW_RATIO)
  NULLIFY(YTEB%XBLD_HEIGHT)
  NULLIFY(YTEB%XWALL_O_HOR)
  NULLIFY(YTEB%XROAD_O_GRND)
  NULLIFY(YTEB%XGARDEN_O_GRND)
  NULLIFY(YTEB%XWALL_O_GRND)
  NULLIFY(YTEB%XWALL_O_BLD)
  NULLIFY(YTEB%XZ0_TOWN)
  NULLIFY(YTEB%XSVF_ROAD)
  NULLIFY(YTEB%XSVF_GARDEN)
  NULLIFY(YTEB%XSVF_WALL)
  NULLIFY(YTEB%XALB_ROOF)
  NULLIFY(YTEB%XEMIS_ROOF)
  NULLIFY(YTEB%XHC_ROOF)
  NULLIFY(YTEB%XTC_ROOF)
  NULLIFY(YTEB%XD_ROOF)
  NULLIFY(YTEB%XALB_ROAD)
  NULLIFY(YTEB%XEMIS_ROAD)
  NULLIFY(YTEB%XHC_ROAD)
  NULLIFY(YTEB%XTC_ROAD)
  NULLIFY(YTEB%XD_ROAD)
  NULLIFY(YTEB%XALB_WALL)
  NULLIFY(YTEB%XEMIS_WALL)
  NULLIFY(YTEB%XHC_WALL)
  NULLIFY(YTEB%XTC_WALL)
  NULLIFY(YTEB%XD_WALL)
  NULLIFY(YTEB%XH_TRAFFIC)
  NULLIFY(YTEB%XLE_TRAFFIC)
  NULLIFY(YTEB%XH_INDUSTRY)
  NULLIFY(YTEB%XLE_INDUSTRY)
  NULLIFY(YTEB%XTI_ROAD)
  NULLIFY(YTEB%XWS_ROOF)
  NULLIFY(YTEB%XWS_ROAD)
  NULLIFY(YTEB%XT_ROOF)
  NULLIFY(YTEB%XT_ROAD)
  NULLIFY(YTEB%XT_WALL_A)
  NULLIFY(YTEB%XT_WALL_B)
  NULLIFY(YTEB%XAC_ROOF)
  NULLIFY(YTEB%XAC_ROAD)
  NULLIFY(YTEB%XAC_WALL)
  NULLIFY(YTEB%XAC_TOP)
  NULLIFY(YTEB%XAC_ROOF_WAT)
  NULLIFY(YTEB%XAC_ROAD_WAT)
  NULLIFY(YTEB%XQSAT_ROOF)
  NULLIFY(YTEB%XQSAT_ROAD)
  NULLIFY(YTEB%XDELT_ROOF)
  NULLIFY(YTEB%XDELT_ROAD)
  NULLIFY(YTEB%XT_CANYON)
  NULLIFY(YTEB%XQ_CANYON)
  NULLIFY(YTEB%XROUGH_ROOF)
  NULLIFY(YTEB%XROUGH_WALL)
  NULLIFY(YTEB%XRESIDENTIAL)
YTEB%XDT_RES=0.
YTEB%XDT_OFF=0.
IF (LHOOK) CALL DR_HOOK("MODD_TEB_N:TEB_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_INIT
!
SUBROUTINE TEB_NP_INIT(YNTEB,KPATCH)
TYPE(TEB_NP_t), INTENT(INOUT) :: YNTEB
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_N:TEB_NP_INIT",0,ZHOOK_HANDLE)
IF (.NOT.ASSOCIATED(YNTEB%AL)) THEN
  ALLOCATE(YNTEB%AL(KPATCH))
  DO JP=1,KPATCH
    CALL TEB_INIT(YNTEB%AL(JP))
  ENDDO
ELSE
  DO JP=1,KPATCH
    CALL TEB_INIT(YNTEB%AL(JP))
  ENDDO
  DEALLOCATE(YNTEB%AL)        
ENDIF
IF (LHOOK) CALL DR_HOOK("MODD_TEB_N:TEB_NP_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_NP_INIT


END MODULE MODD_TEB_n
