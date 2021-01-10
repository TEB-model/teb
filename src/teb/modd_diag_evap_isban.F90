!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!########################
MODULE MODD_DIAG_EVAP_ISBA_n
!########################
!
!!****  *MODD_DIAG_NAT - declaration of packed surface parameters for NAT scheme
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
!!      Original       07/11/03
!!      P. Samuelsson  04/2012   MEB
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE DIAG_EVAP_ISBA_t
!------------------------------------------------------------------------------
!
  LOGICAL :: LSURF_EVAP_BUDGET   ! flag for all terms of evaporation
  LOGICAL :: LWATER_BUDGET       ! flag for isba water budget including input  
                                 ! fluxes (rain and snow) and reservoir tendencies
!
!* variables for each patch
!
  REAL, POINTER, DIMENSION(:) :: XLEG          ! latent heat of evaporation over the ground   (W/m2)
  REAL, POINTER, DIMENSION(:) :: XLEGI         ! surface soil ice sublimation                 (W/m2)
  REAL, POINTER, DIMENSION(:) :: XLEV          ! latent heat of evaporation over vegetation   (W/m2)
  REAL, POINTER, DIMENSION(:) :: XLES          ! latent heat of sublimation over the snow     (W/m2)
  REAL, POINTER, DIMENSION(:) :: XLESL         ! latent heat of evaporation over the snow     (W/m2)
  REAL, POINTER, DIMENSION(:) :: XLER          ! evaporation from canopy water interception   (W/m2)
  REAL, POINTER, DIMENSION(:) :: XLETR         ! evapotranspiration of the vegetation         (W/m2)
  REAL, POINTER, DIMENSION(:) :: XUSTAR       ! friction velocity from snow-free 
!                                               ! surface (ISBA-ES:3-L)                         (m/s  
  REAL, POINTER, DIMENSION(:) :: XSNDRIFT      ! blowing snow sublimation (ES or Crocus)      (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XRESTORE   ! surface energy budget restore term            (W/m2)  
  REAL, POINTER, DIMENSION(:) :: XDRAIN        ! soil drainage flux                           (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XQSB          ! lateral subsurface flux (dif option)         (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XRUNOFF       ! sub-grid and supersaturation runoff          (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XHORT         ! sub-grid Horton runoff from the SGH scheme   (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XRRVEG        !  precipitation intercepted by the vegetation (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XMELT         ! snow melt                                    (kg/m2/s)  
  REAL, POINTER, DIMENSION(:) :: XMELTADV   ! advective energy from snow melt water 
!                                               ! (ISBA-ES:3-L)                                 (W/m2)
  REAL, POINTER, DIMENSION(:) :: XIFLOOD       ! Floodplains infiltration                     (kg/m2/s)      
  REAL, POINTER, DIMENSION(:) :: XPFLOOD       ! Precipitation intercepted by the floodplains (kg/m2/s)      
  REAL, POINTER, DIMENSION(:) :: XLE_FLOOD     ! Floodplains evapotration                     (W/m2)      
  REAL, POINTER, DIMENSION(:) :: XLEI_FLOOD    ! Floodplains evapotration                     (W/m2) 
  REAL, POINTER, DIMENSION(:) :: XICEFLUX    ! ice calving flux                             (kg/m2)  
  REAL, POINTER, DIMENSION(:) :: XDRIP         ! dripping from the vegetation reservoir       (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XIRRIG_FLUX   ! irrigation rate (as soil input)              (kg/m2/s)
!  
  REAL, POINTER, DIMENSION(:) :: XGPP          ! Gross Primary Production                     (kgCO2/m2/s)
  REAL, POINTER, DIMENSION(:) :: XRESP_AUTO    ! Autotrophic respiration                      (kgCO2/m2/s)
  REAL, POINTER, DIMENSION(:) :: XRESP_ECO     ! Ecosystem respiration                        (kgCO2/m2/s)
!
!
  REAL, POINTER, DIMENSION(:) :: XLELITTER     ! MEB: interception evaporation from understory vegetation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLELITTERI    ! MEB: interception evaporation from understory vegetation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XDRIPLIT      ! 
  REAL, POINTER, DIMENSION(:) :: XRRLIT        ! 
!
  REAL, POINTER, DIMENSION(:) :: XLEV_CV        ! MEB: total evapotranspiration from vegetation canopy overstory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLES_CV        ! XLESC MEB: total snow sublimation from vegetation canopy overstory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLETR_CV       ! MEB: transpiration from overstory canopy vegetation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLER_CV        ! MEB: interception evaporation from overstory canopy vegetation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLE_CV         ! XLE_V_C MEB: latent heat flux from vegetation canopy overstory [W/m2]  
  REAL, POINTER, DIMENSION(:) :: XH_CV          ! H_V_C MEB: sensible heat flux from vegetation canopy overstory [W/m2] 
  REAL, POINTER, DIMENSION(:) :: XMELT_CV       ! MEB: snow melt rate from the overstory snow reservoir [kg/m2/s]
  REAL, POINTER, DIMENSION(:) :: XFRZ_CV        ! MEB: snow refreeze rate from the overstory snow reservoir [kg/m2/s]  
!
  REAL, POINTER, DIMENSION(:) :: XLETR_GV       ! MEB: transpiration from understory vegetation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLER_GV        ! MEB: interception evaporation from understory vegetation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLE_GV         ! LEG_C MEB: latent heat flux from understory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XH_GV          ! H_G_C MEB: sensible heat flux from understory [W/m2]  
!
  REAL, POINTER, DIMENSION(:) :: XLE_GN         ! LE_N_C MEB: latent heat flux from the snow on the ground [W/m2]
                                                !      NOTE total latent heat flux from the snowpack
                                                !      possibly includes a contribution from snow covering the canopy  
  REAL, POINTER, DIMENSION(:) :: XEVAP_GN       ! EVAP_N_C MEB: Total evap from snow on the ground to canopy air space  [kg/m2/s]
  REAL, POINTER, DIMENSION(:) :: XH_GN          ! H_N_C MEB: sensible heat flux from the snow on the ground [W/m2]
                                                !      NOTE total sensible heat flux from the snowpack
                                                !      possibly includes a contribution from snow covering the canopy   
  REAL, POINTER, DIMENSION(:) :: XSR_GN         ! MEB: snow unloading rate from the overstory reservoir [kg/m2/s]
  REAL, POINTER, DIMENSION(:) :: XSWDOWN_GN     ! MEB: total shortwave radiation transmitted through the canopy
                                                !      reaching the snowpack/ground understory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLWDOWN_GN     ! MEB: total shortwave radiation transmitted through and emitted by the canopy
                                                !      reaching the snowpack/ground understory (explicit part) [W/m2]                                               
!
  REAL, POINTER, DIMENSION(:) :: XEVAP_G        ! EVAP_G_C MEB: Total evap from ground to canopy air space [kg/m2/s]
  REAL, POINTER, DIMENSION(:) :: XLE_CA         ! LE_C_A MEB: latent heat flux from canopy air space to the atmosphere [W/m2] 
                                                !      NOTE total latent heat flux to the atmosphere also possibly 
                                                !      includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:) :: XH_CA          ! H_C_A MEB: sensible heat flux from canopy air space to the atmosphere [W/m2] 
                                                !      NOTE total sensible heat flux to the atmosphere also possibly 
                                                !      includes a contribution from snow covering the canopy
                                       
!
  REAL, POINTER, DIMENSION(:) :: XSWUP         ! MEB: net *total* (surface) upwelling shortwave radiation to atmosphere [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLWUP         ! MEB: net *total* (surface) upwelling longwave radiation to atmosphere [W/m2]     

  REAL, POINTER, DIMENSION(:) :: XSWNET_V      ! MEB: net vegetation canopy shortwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XSWNET_G      ! MEB: net ground shortwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XSWNET_N      ! MEB: net snow shortwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XSWNET_NS     ! MEB: net snow shortwave radiation for *surface* layer 
                                               !     (i.e. net snow shortwave radiation less absorbed radiation) [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLWNET_V      ! MEB: net vegetation canopy longwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLWNET_G      ! MEB: net ground longwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLWNET_N      ! MEB: net snow longwave radiation [W/m2]
!
!
  REAL, POINTER, DIMENSION(:) :: XDWG          ! liquid soil moisture time tendencies         (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XDWGI         ! solid soil moisture time tendencies          (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XDWR          ! canopy water time tendencies                 (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XDSWE         ! snow water equivalent time tendencies        (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XWATBUD       ! ISBA water budget                            (kg/m2/s)
!
  REAL, POINTER, DIMENSION(:)   :: XRAINFALL     ! input rainfall rate for LWATER_BUDGET        (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XSNOWFALL     ! input snowfall rate for LWATER_BUDGET        (kg/m2/s)
! 
! Snowfree diagnostics
!
  REAL, POINTER, DIMENSION(:) :: XRN_SN_FR
  REAL, POINTER, DIMENSION(:) :: XH_SN_FR
  REAL, POINTER, DIMENSION(:) :: XLEI_SN_FR
  REAL, POINTER, DIMENSION(:) :: XLE_SN_FR
  REAL, POINTER, DIMENSION(:) :: XGFLUX_SN_FR
!  
  REAL, POINTER, DIMENSION(:) :: XLEG_SN_FR
  REAL, POINTER, DIMENSION(:) :: XLEGI_SN_FR
  REAL, POINTER, DIMENSION(:) :: XLEV_SN_FR
  REAL, POINTER, DIMENSION(:) :: XLETR_SN_FR
  REAL, POINTER, DIMENSION(:) :: XUSTAR_SN_FR
  REAL, POINTER, DIMENSION(:) :: XLER_SN_FR
!
!------------------------------------------------------------------------------
!

END TYPE DIAG_EVAP_ISBA_t
!
TYPE DIAG_EVAP_ISBA_NP_t
!
TYPE(DIAG_EVAP_ISBA_t), POINTER :: AL(:)=>NULL()
!
END TYPE DIAG_EVAP_ISBA_NP_t
!
CONTAINS
!
SUBROUTINE DIAG_EVAP_ISBA_NP_INIT(YNDIAG_EVAP_ISBA,KPATCH)
TYPE(DIAG_EVAP_ISBA_NP_t), INTENT(INOUT) :: YNDIAG_EVAP_ISBA
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_EVAP_ISBA_N:DIAG_EVAP_ISBA_NP_INIT",0,ZHOOK_HANDLE)
IF (.NOT.ASSOCIATED(YNDIAG_EVAP_ISBA%AL)) THEN
  ALLOCATE(YNDIAG_EVAP_ISBA%AL(KPATCH))
  DO JP=1,KPATCH
    CALL DIAG_EVAP_ISBA_INIT(YNDIAG_EVAP_ISBA%AL(JP))
  ENDDO
ELSE
  DO JP=1,KPATCH
    CALL DIAG_EVAP_ISBA_INIT(YNDIAG_EVAP_ISBA%AL(JP))
  ENDDO
  DEALLOCATE(YNDIAG_EVAP_ISBA%AL)
ENDIF  
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_EVAP_ISBA_N:DIAG_EVAP_ISBA_NP_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_EVAP_ISBA_NP_INIT
!
SUBROUTINE DIAG_EVAP_ISBA_INIT(DE)
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DE
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_EVAP_ISBA_N:DIAG_EVAP_ISBA_INIT",0,ZHOOK_HANDLE)
  NULLIFY(DE%XLEG)
  NULLIFY(DE%XLEGI)
  NULLIFY(DE%XLEV)
  NULLIFY(DE%XLES)
  NULLIFY(DE%XLESL)
  NULLIFY(DE%XLER)
  NULLIFY(DE%XLETR)
  NULLIFY(DE%XUSTAR)
  NULLIFY(DE%XSNDRIFT)
  NULLIFY(DE%XRESTORE)
  NULLIFY(DE%XDRAIN)
  NULLIFY(DE%XQSB)
  NULLIFY(DE%XRUNOFF)
  NULLIFY(DE%XHORT)
  NULLIFY(DE%XRRVEG)
  NULLIFY(DE%XMELT)
  NULLIFY(DE%XMELTADV)
  NULLIFY(DE%XIFLOOD)
  NULLIFY(DE%XPFLOOD)
  NULLIFY(DE%XLE_FLOOD)
  NULLIFY(DE%XLEI_FLOOD)
!
  NULLIFY(DE%XICEFLUX)  
!
  NULLIFY(DE%XLELITTER)
  NULLIFY(DE%XLELITTERI)
  NULLIFY(DE%XDRIPLIT)
  NULLIFY(DE%XRRLIT)

  NULLIFY(DE%XLEV_CV)
  NULLIFY(DE%XLES_CV)
  NULLIFY(DE%XLETR_CV)
  NULLIFY(DE%XLER_CV)
  NULLIFY(DE%XLE_CV)
  NULLIFY(DE%XH_CV)
  NULLIFY(DE%XMELT_CV)
  NULLIFY(DE%XFRZ_CV)  

  NULLIFY(DE%XLETR_GV)
  NULLIFY(DE%XLER_GV)
  NULLIFY(DE%XLE_GV)  
  NULLIFY(DE%XH_GV)  

  NULLIFY(DE%XLE_GN)  
  NULLIFY(DE%XEVAP_GN)
  NULLIFY(DE%XH_GN)  
  NULLIFY(DE%XSR_GN)  
  NULLIFY(DE%XSWDOWN_GN)
  NULLIFY(DE%XLWDOWN_GN)  

  NULLIFY(DE%XEVAP_G)  
  NULLIFY(DE%XLE_CA)
  NULLIFY(DE%XH_CA)
  
  NULLIFY(DE%XSWUP)
  NULLIFY(DE%XLWUP)
  
  NULLIFY(DE%XSWNET_V)
  NULLIFY(DE%XSWNET_G)
  NULLIFY(DE%XSWNET_N)
  NULLIFY(DE%XSWNET_NS)
  NULLIFY(DE%XLWNET_V)
  NULLIFY(DE%XLWNET_G)
  NULLIFY(DE%XLWNET_N)
!
  NULLIFY(DE%XDRIP)
  NULLIFY(DE%XIRRIG_FLUX)
  NULLIFY(DE%XGPP)
  NULLIFY(DE%XRESP_AUTO)
  NULLIFY(DE%XRESP_ECO)  
  NULLIFY(DE%XDWG)
  NULLIFY(DE%XDWGI)
  NULLIFY(DE%XDWR)
  NULLIFY(DE%XDSWE)
  NULLIFY(DE%XWATBUD)  
!
  NULLIFY(DE%XRAINFALL)
  NULLIFY(DE%XSNOWFALL)
!
  NULLIFY(DE%XRN_SN_FR)
  NULLIFY(DE%XH_SN_FR)
  NULLIFY(DE%XLEI_SN_FR)
  NULLIFY(DE%XLE_SN_FR)
  NULLIFY(DE%XGFLUX_SN_FR)
  NULLIFY(DE%XLEG_SN_FR)
  NULLIFY(DE%XLEGI_SN_FR)
  NULLIFY(DE%XLEV_SN_FR)
  NULLIFY(DE%XLETR_SN_FR)
  NULLIFY(DE%XUSTAR_SN_FR)
  NULLIFY(DE%XLER_SN_FR)
  
DE%LSURF_EVAP_BUDGET=.FALSE.
DE%LWATER_BUDGET=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_EVAP_ISBA_N:DIAG_EVAP_ISBA_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_EVAP_ISBA_INIT
!
END MODULE MODD_DIAG_EVAP_ISBA_n
