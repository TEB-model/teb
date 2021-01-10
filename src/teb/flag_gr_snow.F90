!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE FLAG_GR_SNOW(KFLAG,OMASK,TPSNOW)  
!     ##########################################################
!
!!****  *FLAG_GR_SNOW* - routine to flag snow surface fields
!!
!!    PURPOSE
!!    -------
!       Initialize snow surface fields.
!
!!**  METHOD
!!    ------
!!    
!!    
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
!!
!!    AUTHOR
!!    ------
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       10/2011
!!      P. Samuelsson  07/2014 Added snow albedos
!-----------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
USE MODD_TYPE_SNOW
!
USE MODD_SURF_PAR, ONLY : XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
INTEGER, INTENT(IN) :: KFLAG ! 1 : to put physical values to run ISBA afterwards
!                            ! 2 : to flag with XUNDEF value for points wihtout
LOGICAL, DIMENSION(:), INTENT(IN) :: OMASK ! T: points where snow values 
!                                          !    must be flagged
TYPE(SURF_SNOW), INTENT(INOUT) :: TPSNOW   ! snow characteristics
!
!*       0.2   declarations of local variables
!
REAL            :: ZVAL
INTEGER         :: JLAYER
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FLAG_GR_SNOW',0,ZHOOK_HANDLE)
!
IF (KFLAG==1) THEN
  ZVAL = 0.
ELSEIF (KFLAG==2) THEN 
  ZVAL = XUNDEF
ENDIF
!
  DO JLAYER = 1,TPSNOW%NLAYER
    !
    WHERE(OMASK(:)) TPSNOW%WSNOW(:,JLAYER) = ZVAL
    !
    IF (KFLAG==1) THEN 
      !
      WHERE(OMASK(:)) TPSNOW%RHO  (:,JLAYER) = XUNDEF
      !
      IF (SIZE(TPSNOW%TEMP ) >0) THEN
        WHERE(OMASK(:)) 
          TPSNOW%TEMP (:,JLAYER) = XUNDEF
          TPSNOW%HEAT (:,JLAYER) = XUNDEF
        END WHERE
      ENDIF
      !
      IF (SIZE(TPSNOW%T    ) >0) WHERE(OMASK(:)) TPSNOW%T(:,JLAYER) = XUNDEF
      !
      IF (SIZE(TPSNOW%GRAN1) >0) THEN
        WHERE(OMASK(:)) 
          TPSNOW%GRAN1(:,JLAYER) = XUNDEF
          TPSNOW%GRAN2(:,JLAYER) = XUNDEF
          TPSNOW%HIST (:,JLAYER) = XUNDEF
          TPSNOW%AGE  (:,JLAYER) = XUNDEF
        END WHERE
      END IF
      !
    ENDIF
    !
  ENDDO
  !
  IF (KFLAG==1) THEN
    !
    WHERE(OMASK(:)) TPSNOW%ALB    (:) = XUNDEF
    WHERE(OMASK(:)) TPSNOW%ALBVIS (:) = XUNDEF
    WHERE(OMASK(:)) TPSNOW%ALBNIR (:) = XUNDEF
    WHERE(OMASK(:)) TPSNOW%ALBFIR (:) = XUNDEF
    !
    IF (SIZE(TPSNOW%EMIS ) >0) THEN
      WHERE(OMASK(:))
        TPSNOW%EMIS (:) = XUNDEF
        TPSNOW%TS   (:) = XUNDEF
      END WHERE
    END IF
    !
  ENDIF
!
IF (LHOOK) CALL DR_HOOK('FLAG_GR_SNOW',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE FLAG_GR_SNOW
