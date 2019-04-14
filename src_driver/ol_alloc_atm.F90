!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ############################
    MODULE MODI_OL_ALLOC_ATM
!   ############################
!
INTERFACE
!
SUBROUTINE OL_ALLOC_ATM(KNI,KBANDS,KSCAL)
!
!*      0.1    declarations of arguments
!
INTEGER,            INTENT(IN)  :: KNI       ! grid dimension
INTEGER,            INTENT(IN)  :: KBANDS    ! grid dimension
INTEGER,            INTENT(IN)  :: KSCAL     ! grid dimension
!
END SUBROUTINE OL_ALLOC_ATM
!
!
END INTERFACE
!
!
END MODULE MODI_OL_ALLOC_ATM
!     #################################################################################
SUBROUTINE OL_ALLOC_ATM(KNI,KBANDS,KSCAL)
!     #################################################################################
!
USE MODD_SURF_PAR,       ONLY : XUNDEF

USE MODD_FORC_ATM,  ONLY: CSV       ,&! name of all scalar variables
                          XDIR_ALB  ,&! direct albedo for each band
                          XSCA_ALB  ,&! diffuse albedo for each band
                          XEMIS     ,&! emissivity
                          XTSRAD    ,&! radiative temperature
                          XTSUN     ,&! solar time                    (s from midnight)
                          XZS       ,&! orography                             (m)
                          XZREF     ,&! height of T,q forcing                 (m)
                          XUREF     ,&! height of wind forcing                (m)
                          XTA       ,&! air temperature forcing               (K)
                          XQA       ,&! air specific humidity forcing         (kg/m3)
                          XRHOA     ,&! air density forcing                   (kg/m3)
                          XSV       ,&! scalar variables
                          XU        ,&! zonal wind                            (m/s)
                          XV        ,&! meridian wind                         (m/s)
                          XDIR_SW   ,&! direct  solar radiation (on horizontal surf.)
                          XSCA_SW   ,&! diffuse solar radiation (on horizontal surf.)
                          XSW_BANDS ,&! mean wavelength of each shortwave band (m)
                          XZENITH   ,&! zenithal angle at t      (radian from the vertical)
                          XZENITH2  ,&! zenithal angle at t+1    (radian from the vertical)
                          XAZIM     ,&! azimuthal angle      (radian from North, clockwise)
                          XLW       ,&! longwave radiation (on horizontal surf.)
                          XPS       ,&! pressure at atmospheric model surface (Pa)
                          XPA       ,&! pressure at forcing level             (Pa)
                          XCO2      ,&! CO2 concentration in the air          (kg/kg)
                          XSNOW     ,&! snow precipitation                    (kg/m2/s)
                          XRAIN     ,&! liquid precipitation                  (kg/m2/s)
                          XSFTH     ,&! flux of heat                          (W/m2)
                          XSFTQ     ,&! flux of water vapor                   (kg/m2/s)
                          XSFU      ,&! zonal momentum flux                   (pa)
                          XSFV      ,&! meridian momentum flux                (pa)
                          XSFCO2    ,&! flux of CO2                           (kg/m2/s)
                          XSFTS     ,&! flux of scalar var.                   (kg/m2/s)
                          XPEW_A_COEF ,&! implicit coefficients
                          XPEW_B_COEF ,&! needed if HCOUPLING='I'
                          XPET_A_COEF ,&
                          XPEQ_A_COEF ,&
                          XPET_B_COEF ,&
                          XPEQ_B_COEF

!
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
INTEGER,            INTENT(IN)  :: KNI       ! grid dimension
INTEGER,            INTENT(IN)  :: KBANDS    ! grid dimension
INTEGER,            INTENT(IN)  :: KSCAL     ! grid dimension
!
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
ALLOCATE(CSV       (KSCAL)     )! name of all scalar variables
ALLOCATE(XDIR_ALB  (KNI,KBANDS)) ! direct albedo for each band
ALLOCATE(XSCA_ALB  (KNI,KBANDS)) ! diffuse albedo for each band
ALLOCATE(XEMIS     (KNI)       )! emissivity
ALLOCATE(XTSRAD    (KNI)       )! radiative temperature
ALLOCATE(XTSUN     (KNI)       )! solar time                    (s from midnight)
ALLOCATE(XZREF     (KNI)       )! height of T,q forcing                 (m)
ALLOCATE(XUREF     (KNI)       )! height of wind forcing                (m)
ALLOCATE(XTA       (KNI)       )! air temperature forcing               (K)
ALLOCATE(XQA       (KNI)       )! air specific humidity forcing         (kg/m3)
ALLOCATE(XZS       (KNI)       )! orography                             (m)
ALLOCATE(XRHOA     (KNI)       )! air density forcing                   (kg/m3)
ALLOCATE(XSV       (KNI,KSCAL) ) ! scalar variables
ALLOCATE(XU        (KNI)       )! zonal wind                            (m/s)
ALLOCATE(XV        (KNI)       )! meridian wind                         (m/s)
ALLOCATE(XDIR_SW   (KNI,KBANDS))    ! direct  solar radiation (on horizontal surf.)
ALLOCATE(XSCA_SW   (KNI,KBANDS))    ! diffuse solar radiation (on horizontal surf.)
ALLOCATE(XSW_BANDS (KBANDS)    ) ! mean wavelength of each shortwave band (m)
ALLOCATE(XZENITH   (KNI)       )! zenithal angle       (radian from the vertical)
ALLOCATE(XZENITH2  (KNI)       )! zenithal angle at t+1(radian from the vertical)
ALLOCATE(XAZIM     (KNI)       )! azimuthal angle      (radian from North, clockwise)
ALLOCATE(XLW       (KNI)       )! longwave radiation (on horizontal surf.)
ALLOCATE(XPS       (KNI)       )! pressure at atmospheric model surface (Pa)
ALLOCATE(XPA       (KNI)       )! pressure at forcing level             (Pa)
ALLOCATE(XCO2      (KNI)       )! CO2 concentration in the air          (kg/kg)
ALLOCATE(XSNOW     (KNI)       )! snow precipitation                    (kg/m2/s)
ALLOCATE(XRAIN     (KNI)       )! liquid precipitation                  (kg/m2/s)
ALLOCATE(XSFTH     (KNI)       )! flux of heat                          (W/m2)
ALLOCATE(XSFTQ     (KNI)       )! flux of water vapor                   (kg/m2/s)
ALLOCATE(XSFU      (KNI)       )! zonal momentum flux                   (pa)
ALLOCATE(XSFV      (KNI)       )! meridian momentum flux                (pa)
ALLOCATE(XSFCO2    (KNI)       )! flux of CO2                           (kg/m2/s)
ALLOCATE(XSFTS     (KNI,KSCAL) ) ! flux of scalar var.                  (kg/m2/s)
ALLOCATE(XPEW_A_COEF (KNI)     )! implicit coefficients
ALLOCATE(XPEW_B_COEF (KNI)     )! needed if HCOUPLING='I'
ALLOCATE(XPET_A_COEF (KNI)     )
ALLOCATE(XPEQ_A_COEF (KNI)     )
ALLOCATE(XPET_B_COEF (KNI)     )
ALLOCATE(XPEQ_B_COEF (KNI)     )

CSV       (:)  ='UNDEF '! name of all scalar variables
XDIR_ALB  (:,:)=XUNDEF ! direct albedo for each band
XSCA_ALB  (:,:)=XUNDEF ! diffuse albedo for each band
XEMIS     (:)=XUNDEF ! emissivity
XTSRAD    (:)=XUNDEF ! radiative temperature
XTSUN     (:)=XUNDEF ! solar time                    (s from midnight)
XZREF     (:)=XUNDEF ! height of T,q forcing                 (m)
XUREF     (:)=XUNDEF ! height of wind forcing                (m)
XTA       (:)=XUNDEF ! air temperature forcing               (K)
XQA       (:)=XUNDEF ! air specific humidity forcing         (kg/m3)
XRHOA     (:)=XUNDEF ! air density forcing                   (kg/m3)
XSV       (:,:)=XUNDEF ! scalar variables
XU        (:)=XUNDEF ! zonal wind                            (m/s)
XV        (:)=XUNDEF ! meridian wind                         (m/s)
XDIR_SW   (:,:)=XUNDEF ! direct  solar radiation (on horizontal surf.)
XSCA_SW   (:,:)=XUNDEF ! diffuse solar radiation (on horizontal surf.)
XSW_BANDS (:)=XUNDEF ! mean wavelength of each shortwave band (m)
XZENITH   (:)=XUNDEF ! zenithal angle at t      (radian from the vertical)
XZENITH2  (:)=XUNDEF ! zenithal angle at t+1      (radian from the vertical)
XAZIM     (:)=XUNDEF ! azimuthal angle      (radian from North, clockwise)
XLW       (:)=XUNDEF ! longwave radiation (on horizontal surf.)
XPS       (:)=XUNDEF ! pressure at atmospheric model surface (Pa)
XPA       (:)=XUNDEF ! pressure at forcing level             (Pa)
XZS       (:)=XUNDEF ! atmospheric model orography           (m)
XCO2      (:)=XUNDEF ! CO2 concentration in the air          (kg/kg)
XSNOW     (:)=XUNDEF ! snow precipitation                    (kg/m2/s)
XRAIN     (:)=XUNDEF ! liquid precipitation                  (kg/m2/s)
XSFTH     (:)=XUNDEF ! flux of heat                          (W/m2)
XSFTQ     (:)=XUNDEF ! flux of water vapor                   (kg/m2/s)
XSFU      (:)=XUNDEF ! zonal momentum flux                   (pa)
XSFV      (:)=XUNDEF ! meridian momentum flux                (pa)
XSFCO2    (:)=XUNDEF ! flux of CO2                           (kg/m2/s)
XSFTS     (:,:)=XUNDEF ! flux of scalar var.                 (kg/m2/s)
XPEW_A_COEF (:)=XUNDEF ! implicit coefficients
XPEW_B_COEF (:)=XUNDEF ! needed if HCOUPLING='I'
XPET_A_COEF (:)=XUNDEF
XPEQ_A_COEF (:)=XUNDEF
XPET_B_COEF (:)=XUNDEF
XPEQ_B_COEF (:)=XUNDEF

END SUBROUTINE OL_ALLOC_ATM
