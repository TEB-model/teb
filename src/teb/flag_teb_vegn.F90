!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE FLAG_TEB_VEG_n (PEK, IO, PMASK, KFLAG)
!     ##################################
!
!!****  *FLAG_TEB_VEG_n* - routine to flag ISBA variables where gardens are
!!                            not present
!!                         
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2011
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
!
USE MODD_ISBA_n, ONLY : ISBA_PE_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
!
USE MODD_CO2V_PAR,       ONLY : XANFMINIT, XCONDCTMIN
!                                
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
USE MODI_FLAG_GR_SNOW
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
TYPE(ISBA_PE_t), INTENT(INOUT) :: PEK
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
REAL, DIMENSION(:), INTENT(IN) :: PMASK
!
INTEGER, INTENT(IN) :: KFLAG ! 1 : to put physical values to run ISBA afterwards
!                            ! 2 : to flag with XUNDEF value for points wihtout garden
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL :: ZWR, ZTG, ZWG, ZRESA, ZANFM, ZDEF
INTEGER :: JL1, JL2 ! loop counter on layers
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!
!* 1D physical dimension
!
IF (LHOOK) CALL DR_HOOK('FLAG_TEB_VEG_N',0,ZHOOK_HANDLE)
!
ZWR = XUNDEF
!
IF (KFLAG==1) THEN
  ZTG   = 300.
  ZWG   = 0.5
  ZRESA = 100.
  ZANFM = XANFMINIT
  ZDEF  = 0.
ELSEIF (KFLAG==2) THEN
  ZTG   = XUNDEF
  ZWG   = XUNDEF
  ZRESA = XUNDEF
  ZANFM = XUNDEF
  ZDEF  = XUNDEF
ENDIF
!
!-------------------------------------------------------------------------------
!     
  !
  DO JL1=1,IO%NGROUND_LAYER
    WHERE (PMASK(:)==0.) 
      PEK%XTG (:,JL1) = ZTG
      PEK%XWG (:,JL1) = ZWG
      PEK%XWGI(:,JL1) = ZDEF
    END WHERE
  END DO
  !
  WHERE (PMASK(:)==0.) 
    PEK%XWR  (:) = ZWR
    PEK%XRESA(:) = ZRESA
  END WHERE
  !
  IF (IO%CPHOTO/='NON') THEN
    !
    WHERE (PMASK(:)==0.)
      PEK%XANFM (:) = ZANFM              
      PEK%XAN   (:) = ZDEF
      PEK%XANDAY(:) = ZDEF
      PEK%XLE   (:) = ZDEF
    END WHERE
    !
    IF (IO%CPHOTO=='NIT' .OR. IO%CPHOTO=='NCB') THEN
      !
      WHERE (PMASK(:)==0.) PEK%XLAI(:) = ZDEF
      !
    ELSE IF (IO%CPHOTO=='AST') THEN
      !
      DO JL1=1,SIZE(PEK%XBIOMASS,2)
        WHERE (PMASK(:)==0.)
          PEK%XBIOMASS     (:,JL1) = ZDEF
          PEK%XRESP_BIOMASS(:,JL1) = ZDEF
        END WHERE
      END DO
      !
    END IF
    !
  ENDIF
  !
!
!-------------------------------------------------------------------------------
!
!* Flag snow characteristics
!
 CALL FLAG_GR_SNOW(KFLAG,PMASK(:)==0.,PEK%TSNOW)
!
!
!* snow-free characteristics
!
IF (KFLAG==1) THEN
  WHERE (PMASK(:)==0.) 
    PEK%XSNOWFREE_ALB(:)      = 0.2
    PEK%XSNOWFREE_ALB_VEG(:)  = 0.2
    PEK%XSNOWFREE_ALB_SOIL(:) = 0.2
  END WHERE
ELSEIF (KFLAG==2) THEN
  WHERE (PMASK==0.) 
    PEK%XSNOWFREE_ALB(:)      = XUNDEF
    PEK%XSNOWFREE_ALB_VEG(:)  = XUNDEF
    PEK%XSNOWFREE_ALB_SOIL(:) = XUNDEF
  END WHERE
END IF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('FLAG_TEB_VEG_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE FLAG_TEB_VEG_n
