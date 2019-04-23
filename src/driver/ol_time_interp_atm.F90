!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_OL_TIME_INTERP_ATM
INTERFACE

!**************************************************************************
SUBROUTINE OL_TIME_INTERP_ATM (KSURF_STEP,KNB_ATM,                               &
                               PTA,PQA,PWIND,PDIR_SW,PSCA_SW,PLW,PSNOW,PRAIN,PPS,&
                               PCO2,PDIR             )
!
!**************************************************************************
IMPLICIT NONE
!
! global variables
REAL, DIMENSION(:,:),INTENT(IN) :: PTA,PQA,PWIND,PDIR_SW,PSCA_SW,PLW,PSNOW,PRAIN,PPS,PCO2,PDIR
INTEGER,INTENT(IN)   :: KSURF_STEP,KNB_ATM
!
END SUBROUTINE OL_TIME_INTERP_ATM
END INTERFACE
END MODULE MODI_OL_TIME_INTERP_ATM
!
!**************************************************************************
SUBROUTINE OL_TIME_INTERP_ATM (KSURF_STEP,KNB_ATM,                               &
                               PTA,PQA,PWIND,PDIR_SW,PSCA_SW,PLW,PSNOW,PRAIN,PPS,&
                               PCO2,PDIR )
!**************************************************************************
!
!!    PURPOSE
!!    -------
!        Time interpolation of the atmospheric forcing
!        So far, it is a simple linear interpolation.
!        More complex interpolation may be added, especially for the atmospheric
!        radiation (Option to use).
!        Output are in the module 
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
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
!!	F. Habets   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/2003

!          
!
USE MODD_CSTS,       ONLY : XPI, XRD, XRV, XG
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_FORC_ATM,  ONLY: XTA       ,&! air temperature forcing               (K)
                          XQA       ,&! air specific humidity forcing         (kg/m3)
                          XRHOA     ,&! air density forcing                   (kg/m3)
                          XZS       ,&! orography                             (m)
                          XU        ,&! zonal wind                            (m/s)
                          XV        ,&! meridian wind                         (m/s)
                          XDIR_SW   ,&! direct  solar radiation (on horizontal surf.)
                          XSCA_SW   ,&! diffuse solar radiation (on horizontal surf.)
                          XLW       ,&! longwave radiation (on horizontal surf.)
                          XPS       ,&! pressure at atmospheric model surface (Pa)
                          XPA       ,&! pressure at forcing level             (Pa)
                          XRHOA     ,&! density at forcing level              (kg/m3)
                          XCO2      ,&! CO2 concentration in the air          (kg/kg)
                          XSNOW     ,&! snow precipitation                    (kg/m2/s)
                          XRAIN     ,&! liquid precipitation                  (kg/m2/s)
                          XZREF       ! height of T,q forcing                 (m)

!
IMPLICIT NONE
!

! global variables
REAL, DIMENSION(:,:),INTENT(IN) :: PTA,PQA,PWIND,PDIR_SW,PSCA_SW,PLW,PSNOW,PRAIN,PPS,PCO2,PDIR
INTEGER,INTENT(IN)   :: KSURF_STEP,KNB_ATM

! local variables
REAL, DIMENSION(SIZE(PTA,2)) :: ZDELTA_TA,ZDELTA_QA,                   &
                                ZDELTA_DIR_SW,ZDELTA_SCA_SW,ZDELTA_LW, &
                                ZDELTA_SNOW,ZDELTA_RAIN,ZDELTA_PS,     &
                                ZDELTA_CO2,ZDELTA_U, ZDELTA_V,         &
                                ZU1, ZV1, ZU2, ZV2
INTEGER :: ILUOUT
!========================================================================
!
ILUOUT = 10
!
! Compute wind components
! zonal wind
ZU1 = PWIND  (1,:) * SIN(PDIR(1,:)*XPI/180.)
ZV1 = PWIND  (1,:) * COS(PDIR(1,:)*XPI/180.)
ZU2 = PWIND  (2,:) * SIN(PDIR(2,:)*XPI/180.)
ZV2 = PWIND  (2,:) * COS(PDIR(2,:)*XPI/180.)
!
! Compute variation from atmospheric time step J and J+1
!
ZDELTA_TA    (:)=(PTA    (2,:)-PTA    (1,:))/(1.*KNB_ATM)
ZDELTA_QA    (:)=(PQA    (2,:)-PQA    (1,:))/(1.*KNB_ATM)
ZDELTA_U     (:)=(ZU2      (:)-ZU1      (:))/(1.*KNB_ATM)
ZDELTA_V     (:)=(ZV2      (:)-ZV1      (:))/(1.*KNB_ATM)
ZDELTA_RAIN  (:)=(PRAIN  (2,:)-PRAIN  (1,:))/(1.*KNB_ATM)
ZDELTA_SNOW  (:)=(PSNOW  (2,:)-PSNOW  (1,:))/(1.*KNB_ATM)
ZDELTA_DIR_SW(:)=(PDIR_SW(2,:)-PDIR_SW(1,:))/(1.*KNB_ATM)
ZDELTA_SCA_SW(:)=(PSCA_SW(2,:)-PSCA_SW(1,:))/(1.*KNB_ATM)
ZDELTA_LW    (:)=(PLW    (2,:)-PLW    (1,:))/(1.*KNB_ATM)
ZDELTA_PS    (:)=(PPS    (2,:)-PPS    (1,:))/(1.*KNB_ATM)
ZDELTA_CO2   (:)=(PCO2   (2,:)-PCO2   (1,:))/(1.*KNB_ATM)
!
XRAIN (:)= PRAIN(1,:)
XSNOW (:)= PSNOW(1,:)
!
XTA    (:)= PTA    (1,:)+(KSURF_STEP*1.-1.)*ZDELTA_TA    (:)
XQA    (:)= PQA    (1,:)+(KSURF_STEP*1.-1.)*ZDELTA_QA    (:)
XU     (:)= ZU1      (:)+(KSURF_STEP*1.-1.)*ZDELTA_U     (:)
XV     (:)= ZV1      (:)+(KSURF_STEP*1.-1.)*ZDELTA_V     (:)
XLW    (:)= PLW    (1,:)+(KSURF_STEP*1.-1.)*ZDELTA_LW    (:)
XPS    (:)= PPS    (1,:)+(KSURF_STEP*1.-1.)*ZDELTA_PS    (:)
XCO2   (:)= PCO2   (1,:)+(KSURF_STEP*1.-1.)*ZDELTA_CO2   (:)
!
! Only 1 band
XDIR_SW(:,1)= PDIR_SW(1,:)+(KSURF_STEP*1.-1.)*ZDELTA_DIR_SW(:)
XSCA_SW(:,1)= PSCA_SW(1,:)+(KSURF_STEP*1.-1.)*ZDELTA_SCA_SW(:)
!
! air density
XRHOA (:) = XPS (:) / XRD /  XTA (:) / ( 1.+((XRV/XRD)-1.)*XQA (:) )
XRHOA(:) = XPS(:) / ( XTA(:)*XRD * ( 1.+((XRV/XRD)-1.)*XQA(:) ) + XZREF(:)*XG )
!
! humidity in kg/m3
XQA(:) = XQA(:) * XRHOA(:)
!
! Check No value data
!---------------------
! Error cases
!
IF ((MINVAL(XTA)  .EQ.XUNDEF).OR.(MINVAL(XQA).EQ.XUNDEF).OR.&
    (MINVAL(XU).EQ.XUNDEF).OR.(MINVAL(XRAIN).EQ.XUNDEF).OR.&
    (MINVAL(XSNOW).EQ.XUNDEF)) THEN
    WRITE(ILUOUT,*)'MINVAL(XTA),MINVAL(XQA),MINVAL(XU),MINVAL(XRAIN),MINVAL(XSNOW)'
    WRITE(ILUOUT,*)MINVAL(XTA),MINVAL(XQA),MINVAL(XU),MINVAL(XRAIN),MINVAL(XSNOW)
    CALL ABOR1_SFX('OL_TIME_INTERP_ATM: UNDEFINED VALUE IN ATMOSPHERIC FORCING')
ENDIF
!
IF ((MINVAL(XDIR_SW).EQ.XUNDEF).AND.(MINVAL(XSCA_SW).EQ.XUNDEF)) THEN
    WRITE(ILUOUT,*)'MINVAL(XSCA_SW),MINVAL(XDIR_SW)'
    WRITE(ILUOUT,*)MINVAL(XSCA_SW),MINVAL(XDIR_SW)
    CALL ABOR1_SFX('OL_TIME_INTERP_ATM: UNDEFINED VALUE IN ATMOSPHERIC FORCING')
ENDIF
!
IF ((MINVAL(XPS).EQ.XUNDEF).AND.(MINVAL(XZS).EQ.XUNDEF)) THEN
    WRITE(ILUOUT,*)'MINVAL(XPS),MINVAL(XZS)'
    WRITE(ILUOUT,*)MINVAL(XPS),MINVAL(XZS)
    CALL ABOR1_SFX('OL_TIME_INTERP_ATM: UNDEFINED VALUE IN ATMOSPHERIC FORCING')
ENDIF
!
IF (MINVAL(XDIR_SW).EQ.XUNDEF) XDIR_SW(:,:)=0. ! No direct solar radiation
IF (MINVAL(XSCA_SW).EQ.XUNDEF) XSCA_SW(:,:)=0. ! No diffuse solar radiation
IF (MINVAL(XPS)    .EQ.XUNDEF) THEN            ! No surface Pressure 
   WRITE(ILUOUT,*)' OL_TIME_INTERP_ATM: SURFACE PRESSURE COMPUTED FROM ZS'
   XPS(:)  = 101325*(1-0.0065 * XZS(:)/288.15)**5.31
ENDIF
!
!* forcing level pressure from hydrostatism
XPA = XPS - XRHOA * XG * XZREF
!
END SUBROUTINE OL_TIME_INTERP_ATM

