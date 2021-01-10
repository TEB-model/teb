!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ################
      MODULE MODD_DATA_BEM_n
!     ################
!
!!****  *MODD_DATA_BEM_n - declaration of surface parameters for BEM
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
!!      G. Pigeon   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       08/2011
!
!*       0.   DECLARATIONS
!             ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE DATA_BEM_t
  !
  LOGICAL :: LDATA_TCOOL_TARGET 
  LOGICAL :: LDATA_THEAT_TARGET   
  LOGICAL :: LDATA_F_WASTE_CAN   
  LOGICAL :: LDATA_EFF_HEAT   
  LOGICAL :: LDATA_HC_FLOOR     
  LOGICAL :: LDATA_TC_FLOOR     
  LOGICAL :: LDATA_D_FLOOR 
  LOGICAL :: LDATA_QIN          
  LOGICAL :: LDATA_QIN_FRAD     
  LOGICAL :: LDATA_SHGC         
  LOGICAL :: LDATA_U_WIN        
  LOGICAL :: LDATA_GR           
  LOGICAL :: LDATA_SHGC_SH
  LOGICAL :: LDATA_FLOOR_HEIGHT 
  LOGICAL :: LDATA_INF  
  LOGICAL :: LDATA_F_WATER_COND  
  LOGICAL :: LDATA_SHADE
  LOGICAL :: LDATA_NATVENT
  LOGICAL :: LDATA_QIN_FLAT 
  LOGICAL :: LDATA_HR_TARGET
  LOGICAL :: LDATA_V_VENT   
  LOGICAL :: LDATA_CAP_SYS_HEAT 
  LOGICAL :: LDATA_CAP_SYS_RAT  
  LOGICAL :: LDATA_T_ADP        
  LOGICAL :: LDATA_M_SYS_RAT    
  LOGICAL :: LDATA_COP_RAT      
  LOGICAL :: LDATA_T_SIZE_MAX    
  LOGICAL :: LDATA_T_SIZE_MIN  
!
! Number of layers in the specification of thermal characteristics
!
  INTEGER                    :: NPAR_FLOOR_LAYER   ! number of layers in roofs
!
!
! Parameters:
!
  REAL, POINTER, DIMENSION(:)   :: XPAR_TCOOL_TARGET ! cooling setpoint of indoor air
  REAL, POINTER, DIMENSION(:)   :: XPAR_THEAT_TARGET ! heating setpoint of indoor air
  REAL, POINTER, DIMENSION(:)   :: XPAR_F_WASTE_CAN  ! fraction of waste heat into the canyon
  REAL, POINTER, DIMENSION(:)   :: XPAR_EFF_HEAT     ! efficiency of the heating system
  REAL, POINTER, DIMENSION(:,:) :: XPAR_HC_FLOOR     ! heat capacity of floor layers [J m-3 K-1]
  REAL, POINTER, DIMENSION(:,:) :: XPAR_TC_FLOOR     ! thermal conductivity of floor layers [W m-1 K-1]
  REAL, POINTER, DIMENSION(:,:) :: XPAR_D_FLOOR      ! thickness of floor layers [m]
  REAL, POINTER, DIMENSION(:)   :: XPAR_QIN          ! internal heat gains [W m-2(floor)]
  REAL, POINTER, DIMENSION(:)   :: XPAR_QIN_FRAD     ! radiant fraction of internal heat gains
  REAL, POINTER, DIMENSION(:)   :: XPAR_SHGC         ! solar transmitance of windows
  REAL, POINTER, DIMENSION(:)   :: XPAR_U_WIN        ! glazing thermal resistance [K m W-2]
  REAL, POINTER, DIMENSION(:)   :: XPAR_GR           ! glazing ratio
  REAL, POINTER, DIMENSION(:)   :: XPAR_SHGC_SH      ! solar transmitance of windows + shading
  REAL, POINTER, DIMENSION(:)   :: XPAR_FLOOR_HEIGHT ! building floor height [m]
  REAL, POINTER, DIMENSION(:)   :: XPAR_INF          ! infiltration/ventilation flow rate [AC/H]
!
  REAL, POINTER, DIMENSION(:)   :: XPAR_F_WATER_COND ! fraction of evaporation for condensers
  REAL, POINTER, DIMENSION(:)   :: XPAR_SHADE        ! Flag to activate shading devices 0->inactivated , 1->activated
  REAL, POINTER, DIMENSION(:)   :: XPAR_NATVENT      ! Flag to describe bld surventilation solution 0-> NONE ; 1 -> MANU ; 2-> AUTO
  REAL, POINTER, DIMENSION(:)   :: XPAR_QIN_FLAT     ! Latent franction of internal heat gains
  REAL, POINTER, DIMENSION(:)   :: XPAR_HR_TARGET    ! Relative humidity setpoint
  REAL, POINTER, DIMENSION(:)   :: XPAR_V_VENT       ! Ventilation flow rate [AC/H]
  REAL, POINTER, DIMENSION(:)   :: XPAR_CAP_SYS_HEAT ! Capacity of the heating system 
                                                     ! [W m-2(bld)]
  REAL, POINTER, DIMENSION(:)   :: XPAR_CAP_SYS_RAT  ! Rated capacity of the cooling system
                                                     ! [W m-2(bld)]
  REAL, POINTER, DIMENSION(:)   :: XPAR_T_ADP        ! Apparatus dewpoint temperature of the
                                                     ! cooling coil [K]
  REAL, POINTER, DIMENSION(:)   :: XPAR_M_SYS_RAT    ! Rated HVAC mass flow rate 
                                                     ! [kg s-1 m-2(bld)]
  REAL, POINTER, DIMENSION(:)   :: XPAR_COP_RAT      ! Rated COP of the cooling system
  REAL, POINTER, DIMENSION(:)   :: XPAR_T_SIZE_MAX   ! 
  REAL, POINTER, DIMENSION(:)   :: XPAR_T_SIZE_MIN   ! 
!
END TYPE DATA_BEM_t



CONTAINS

!
!


!
!
SUBROUTINE DATA_BEM_INIT(YDATA_BEM)
TYPE(DATA_BEM_t), INTENT(INOUT) :: YDATA_BEM
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DATA_BEM_N:DATA_BEM_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YDATA_BEM%XPAR_F_WASTE_CAN)  
  NULLIFY(YDATA_BEM%XPAR_EFF_HEAT)
  NULLIFY(YDATA_BEM%XPAR_HC_FLOOR)
  NULLIFY(YDATA_BEM%XPAR_TC_FLOOR)
  NULLIFY(YDATA_BEM%XPAR_D_FLOOR)
  NULLIFY(YDATA_BEM%XPAR_QIN)
  NULLIFY(YDATA_BEM%XPAR_QIN_FRAD)
  NULLIFY(YDATA_BEM%XPAR_SHGC) 
  NULLIFY(YDATA_BEM%XPAR_U_WIN)
  NULLIFY(YDATA_BEM%XPAR_GR)
  NULLIFY(YDATA_BEM%XPAR_FLOOR_HEIGHT)
  NULLIFY(YDATA_BEM%XPAR_INF)
  NULLIFY(YDATA_BEM%XPAR_F_WATER_COND)
  NULLIFY(YDATA_BEM%XPAR_QIN_FLAT)
  NULLIFY(YDATA_BEM%XPAR_V_VENT)
  NULLIFY(YDATA_BEM%XPAR_CAP_SYS_HEAT)
  NULLIFY(YDATA_BEM%XPAR_CAP_SYS_RAT)
  NULLIFY(YDATA_BEM%XPAR_T_ADP)
  NULLIFY(YDATA_BEM%XPAR_M_SYS_RAT)
  NULLIFY(YDATA_BEM%XPAR_COP_RAT)
  NULLIFY(YDATA_BEM%XPAR_T_SIZE_MAX)
  NULLIFY(YDATA_BEM%XPAR_T_SIZE_MIN)
  NULLIFY(YDATA_BEM%XPAR_SHADE)
  NULLIFY(YDATA_BEM%XPAR_NATVENT)
YDATA_BEM%LDATA_F_WASTE_CAN  = .FALSE.   
YDATA_BEM%LDATA_EFF_HEAT     = .FALSE.   
YDATA_BEM%LDATA_HC_FLOOR     = .FALSE.     
YDATA_BEM%LDATA_TC_FLOOR     = .FALSE.     
YDATA_BEM%LDATA_D_FLOOR      = .FALSE. 
YDATA_BEM%LDATA_QIN          = .FALSE.          
YDATA_BEM%LDATA_QIN_FRAD     = .FALSE.     
YDATA_BEM%LDATA_SHGC         = .FALSE.         
YDATA_BEM%LDATA_U_WIN        = .FALSE.        
YDATA_BEM%LDATA_GR           = .FALSE.           
YDATA_BEM%LDATA_SHGC_SH      = .FALSE.
YDATA_BEM%LDATA_FLOOR_HEIGHT = .FALSE. 
YDATA_BEM%LDATA_INF          = .FALSE.  
YDATA_BEM%LDATA_F_WATER_COND = .FALSE.  
YDATA_BEM%LDATA_SHADE        = .FALSE.
YDATA_BEM%LDATA_NATVENT      = .FALSE.
YDATA_BEM%LDATA_V_VENT       = .FALSE.   
YDATA_BEM%LDATA_QIN_FLAT     = .FALSE.     
YDATA_BEM%LDATA_CAP_SYS_HEAT = .FALSE. 
YDATA_BEM%LDATA_CAP_SYS_RAT  = .FALSE.  
YDATA_BEM%LDATA_T_ADP        = .FALSE.        
YDATA_BEM%LDATA_M_SYS_RAT    = .FALSE.    
YDATA_BEM%LDATA_COP_RAT      = .FALSE.      
YDATA_BEM%LDATA_T_SIZE_MAX   = .FALSE.    
YDATA_BEM%LDATA_T_SIZE_MIN   = .FALSE.
YDATA_BEM%NPAR_FLOOR_LAYER   =0

IF (LHOOK) CALL DR_HOOK("MODD_DATA_BEM_N:DATA_BEM_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DATA_BEM_INIT
!
!
END MODULE MODD_DATA_BEM_n
