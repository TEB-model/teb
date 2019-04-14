!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ################
      MODULE MODD_BEM_n
!     ################
!
!!****  *MODD_BEM_n - declaration of parameters and option for BEM
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
!!      B. Bueno   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       10/2010
!!      G. Pigeon      06/2011 add LSHAD_DAY
!!      G. Pigeon      07/2011 add LNATVENT_NIGHT
!!      G. Pigeon      08/2011 change from MODD_BLD -> MODD_BEM
!!      G. Pigeon      10/2011 add indoor relative surf. and view factors
!!      G. Pigeon      09/2012 add TRAN_WIN
!!      G. Pigeon      10/2012 add XF_WIN_WIN
!!      V. Masson      06/2013 splits module in two
!
!*       0.   DECLARATIONS
!             ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!--------------------------------------------------------------------------
!
TYPE BEM_t
!
! Floor parameters
!
  REAL, POINTER, DIMENSION(:,:) :: XHC_FLOOR     ! floor layers heat capacity        (J/K/m3)
  REAL, POINTER, DIMENSION(:,:) :: XTC_FLOOR     ! floor layers thermal conductivity (W/K/m)
  REAL, POINTER, DIMENSION(:,:) :: XD_FLOOR      ! depth of floor layers             (m)
!
! HVAC parameters
!
  REAL, POINTER, DIMENSION(:)   :: XTCOOL_TARGET ! cooling setpoint of indoor air
  REAL, POINTER, DIMENSION(:)   :: XTHEAT_TARGET ! heating setpoint of indoor air
  REAL, POINTER, DIMENSION(:)   :: XF_WASTE_CAN  ! fraction of waste heat released into the canyon
  REAL, POINTER, DIMENSION(:)   :: XEFF_HEAT     ! efficiency of the heating system
!
! Indoor parameters
!
  REAL, POINTER, DIMENSION(:)   :: XTI_BLD       ! building interior temperature    (K)
  REAL, POINTER, DIMENSION(:,:) :: XT_FLOOR      ! floor layer temperatures         (K)
  REAL, POINTER, DIMENSION(:,:) :: XT_MASS       ! Air cooled building internal th. mass temperature (K)
!
  REAL, POINTER, DIMENSION(:)   :: XQIN          ! internal heat gains [W m-2(floor)]
  REAL, POINTER, DIMENSION(:)   :: XQIN_FRAD     ! radiant fraction of internal heat gains
  REAL, POINTER, DIMENSION(:)   :: XSHGC         ! solar heat gain coef. of windows
  REAL, POINTER, DIMENSION(:)   :: XSHGC_SH      ! solar heat gain coef. of windows + shading
  REAL, POINTER, DIMENSION(:)   :: XU_WIN        ! window U-factor [K m W-2]
  REAL, POINTER, DIMENSION(:)   :: XTRAN_WIN     ! window transmittance (-)
  REAL, POINTER, DIMENSION(:)   :: XGR           ! glazing ratio
  REAL, POINTER, DIMENSION(:)   :: XFLOOR_HEIGHT ! building floor height [m]
  REAL, POINTER, DIMENSION(:)   :: XINF          ! infiltration/ventilation flow rate [AC/H]
!
! New parameters
!
  REAL, POINTER, DIMENSION(:)   :: XF_WATER_COND  ! fraction of evaporation for condensers (cooling system)
  REAL, POINTER, DIMENSION(:)   :: XAUX_MAX      ! Auxiliar variable for autosize calcs
  REAL, POINTER, DIMENSION(:)   :: XQIN_FLAT     ! Latent franction of internal heat gains
  REAL, POINTER, DIMENSION(:)   :: XHR_TARGET    ! Relative humidity setpoint
  REAL, POINTER, DIMENSION(:)   :: XT_WIN2       ! Indoor window temperature [K]
  REAL, POINTER, DIMENSION(:)   :: XQI_BLD       ! Indoor air specific humidity [kg kg-1]
  REAL, POINTER, DIMENSION(:)   :: XV_VENT       ! Ventilation flow rate [AC/H]
  REAL, POINTER, DIMENSION(:)   :: XCAP_SYS_HEAT ! Capacity of the heating system 
                                                 ! [W m-2(bld)]
  REAL, POINTER, DIMENSION(:)   :: XCAP_SYS_RAT  ! Rated capacity of the cooling system
                                                 ! [W m-2(bld)]
  REAL, POINTER, DIMENSION(:)   :: XT_ADP        ! Apparatus dewpoint temperature of the
                                                 ! cooling coil [K]
  REAL, POINTER, DIMENSION(:)   :: XM_SYS_RAT    ! Rated HVAC mass flow rate 
                                                 ! [kg s-1 m-2(bld)]
  REAL, POINTER, DIMENSION(:)   :: XCOP_RAT      ! Rated COP of the cooling system
  REAL, POINTER, DIMENSION(:)   :: XT_WIN1       ! outdoor window temperature [K]
  REAL, POINTER, DIMENSION(:)   :: XALB_WIN      ! window albedo
  REAL, POINTER, DIMENSION(:)   :: XABS_WIN      ! window absortance
  REAL, POINTER, DIMENSION(:)   :: XT_SIZE_MAX   ! Maximum outdoor air temperature for
                                                 ! HVAC sizing [K]
  REAL, POINTER, DIMENSION(:)   :: XT_SIZE_MIN   ! Minimum outdoor air temperature for
                                                 ! HVAC sizing [K]
  REAL, POINTER, DIMENSION(:)   :: XUGG_WIN      ! Window glass-to-glass U-factor [K m W-2]
  LOGICAL, POINTER, DIMENSION(:):: LSHADE        ! flag to activate shading devices -> LOGICAL in the code
  REAL,    POINTER, DIMENSION(:):: XSHADE        ! flag to activate shading devices -> REAL for i/o 0. or 1.
  CHARACTER(LEN=4), POINTER, DIMENSION(:) :: CNATVENT ! flag to activate natural ventilation 'NONE', 'MANU', 'AUTO'
  REAL,    POINTER, DIMENSION(:):: XNATVENT      ! flag to describe surventilation system for i/o 
                                                 ! 0 for NONE, 1 for MANU and 2 for AUTO
  LOGICAL, POINTER, DIMENSION(:):: LSHAD_DAY     !Has shading been necessary this day ?
  LOGICAL, POINTER, DIMENSION(:):: LNATVENT_NIGHT !Has nocturnal surventilation been necessary and possible this night ?
  !
  !indoor relative surfaces and view factors
  REAL, POINTER, DIMENSION(:) :: XN_FLOOR        ! Number of floors     
  REAL, POINTER, DIMENSION(:) :: XGLAZ_O_BLD    ! Window area [m2_win/m2_bld]
  REAL, POINTER, DIMENSION(:) :: XMASS_O_BLD    ! Mass area [m2_mass/m2_bld]
  REAL, POINTER, DIMENSION(:) :: XFLOOR_HW_RATIO ! H/W ratio of 1 floor level
  REAL, POINTER, DIMENSION(:) :: XF_FLOOR_MASS   ! View factor floor-mass
  REAL, POINTER, DIMENSION(:) :: XF_FLOOR_WALL   ! View factor floor-wall
  REAL, POINTER, DIMENSION(:) :: XF_FLOOR_WIN    ! View factor floor-window
  REAL, POINTER, DIMENSION(:) :: XF_FLOOR_ROOF   ! View factor floor-roof
  REAL, POINTER, DIMENSION(:) :: XF_WALL_FLOOR   ! View factor wall-floor
  REAL, POINTER, DIMENSION(:) :: XF_WALL_MASS    ! View factor wall-mass
  REAL, POINTER, DIMENSION(:) :: XF_WALL_WIN     ! View factor wall-win
  REAL, POINTER, DIMENSION(:) :: XF_WIN_FLOOR    ! View factor win-floor
  REAL, POINTER, DIMENSION(:) :: XF_WIN_MASS     ! View factor win-mass
  REAL, POINTER, DIMENSION(:) :: XF_WIN_WALL     ! View factor win-wall
  REAL, POINTER, DIMENSION(:) :: XF_WIN_WIN      ! indoor View factor win-win
  REAL, POINTER, DIMENSION(:) :: XF_MASS_FLOOR   ! View factor mass-floor
  REAL, POINTER, DIMENSION(:) :: XF_MASS_WALL    ! View factor mass-wall
  REAL, POINTER, DIMENSION(:) :: XF_MASS_WIN     ! View factor mass-window


! 
END TYPE BEM_t
!
TYPE BEM_NP_t
  !
  TYPE(BEM_t), POINTER :: AL(:) => NULL()
  !
END TYPE BEM_NP_t
!
CONTAINS

!----------------------------------------------------------------------------
!
SUBROUTINE BEM_INIT(YBEM)
TYPE(BEM_t), INTENT(INOUT) :: YBEM
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_BEM_N:BEM_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YBEM%XF_WATER_COND)
  NULLIFY(YBEM%XHC_FLOOR)
  NULLIFY(YBEM%XTC_FLOOR)
  NULLIFY(YBEM%XD_FLOOR)
  NULLIFY(YBEM%XTCOOL_TARGET)
  NULLIFY(YBEM%XTHEAT_TARGET)
  NULLIFY(YBEM%XTI_BLD)
  NULLIFY(YBEM%XT_FLOOR)
  NULLIFY(YBEM%XT_MASS)
  NULLIFY(YBEM%XQIN)
  NULLIFY(YBEM%XQIN_FRAD)
  NULLIFY(YBEM%XSHGC)
  NULLIFY(YBEM%XSHGC_SH)
  NULLIFY(YBEM%XU_WIN)
  NULLIFY(YBEM%XTRAN_WIN)
  NULLIFY(YBEM%XGR)
  NULLIFY(YBEM%XFLOOR_HEIGHT)
  NULLIFY(YBEM%XEFF_HEAT)
  NULLIFY(YBEM%XINF)
  NULLIFY(YBEM%XF_WASTE_CAN)
  NULLIFY(YBEM%XAUX_MAX)
  NULLIFY(YBEM%XQIN_FLAT)
  NULLIFY(YBEM%XHR_TARGET)
  NULLIFY(YBEM%XT_WIN2)
  NULLIFY(YBEM%XQI_BLD)
  NULLIFY(YBEM%XV_VENT)
  NULLIFY(YBEM%XCAP_SYS_HEAT)
  NULLIFY(YBEM%XCAP_SYS_RAT)
  NULLIFY(YBEM%XT_ADP)
  NULLIFY(YBEM%XM_SYS_RAT)
  NULLIFY(YBEM%XCOP_RAT)
  NULLIFY(YBEM%XT_WIN1)
  NULLIFY(YBEM%XALB_WIN)
  NULLIFY(YBEM%XABS_WIN)
  NULLIFY(YBEM%XT_SIZE_MAX)
  NULLIFY(YBEM%XT_SIZE_MIN)
  NULLIFY(YBEM%XUGG_WIN)
  NULLIFY(YBEM%LSHAD_DAY)
  NULLIFY(YBEM%LNATVENT_NIGHT)
  NULLIFY(YBEM%LSHADE)
  NULLIFY(YBEM%XSHADE)
  NULLIFY(YBEM%CNATVENT)
  NULLIFY(YBEM%XNATVENT)
  NULLIFY(YBEM%XN_FLOOR)
  NULLIFY(YBEM%XGLAZ_O_BLD)
  NULLIFY(YBEM%XMASS_O_BLD)
  NULLIFY(YBEM%XFLOOR_HW_RATIO)
  NULLIFY(YBEM%XF_FLOOR_MASS)
  NULLIFY(YBEM%XF_FLOOR_WALL)
  NULLIFY(YBEM%XF_FLOOR_WIN)
  NULLIFY(YBEM%XF_FLOOR_ROOF)
  NULLIFY(YBEM%XF_WALL_FLOOR)
  NULLIFY(YBEM%XF_WALL_MASS)
  NULLIFY(YBEM%XF_WALL_WIN)
  NULLIFY(YBEM%XF_WIN_FLOOR)
  NULLIFY(YBEM%XF_WIN_MASS)
  NULLIFY(YBEM%XF_WIN_WALL)
  NULLIFY(YBEM%XF_WIN_WIN)
  NULLIFY(YBEM%XF_MASS_FLOOR)
  NULLIFY(YBEM%XF_MASS_WALL)
  NULLIFY(YBEM%XF_MASS_WIN)
IF (LHOOK) CALL DR_HOOK("MODD_BEM_N:BEM_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE BEM_INIT
!
SUBROUTINE BEM_NP_INIT(YNBEM,KPATCH)
TYPE(BEM_NP_t), INTENT(INOUT) :: YNBEM
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_BEM_N:BEM_NP_INIT",0,ZHOOK_HANDLE)
IF (.NOT.ASSOCIATED(YNBEM%AL)) THEN
  ALLOCATE(YNBEM%AL(KPATCH))
  DO JP=1,KPATCH
    CALL BEM_INIT(YNBEM%AL(JP))
  ENDDO
ELSE
  DO JP=1,KPATCH
    CALL BEM_INIT(YNBEM%AL(JP))
  ENDDO
  DEALLOCATE(YNBEM%AL)  
ENDIF 
IF (LHOOK) CALL DR_HOOK("MODD_BEM_N:BEM_NP_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE BEM_NP_INIT
!
!
END MODULE MODD_BEM_n
