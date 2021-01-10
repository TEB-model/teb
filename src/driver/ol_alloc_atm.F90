!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE OL_ALLOC_ATM(KNI,KBANDS,KSCAL)
!     #################################################################################
!
!!
!!    MODIFICATIONS
!!    -------------
!     05/2013  B. Decharme : New coupling variables (for AGCM)
!-------------------------------------------------------------------------------
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
                            XPEQ_B_COEF ,&
                            XTSURF    ,&
                            XZ0       ,&
                            XZ0H      ,&
                            XQSURF
!
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
INTEGER,            INTENT(IN)  :: KNI       ! grid dimension
INTEGER,            INTENT(IN)  :: KBANDS    ! grid dimension
INTEGER,            INTENT(IN)  :: KSCAL     ! grid dimension
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('OL_ALLOC_ATM',0,ZHOOK_HANDLE)
IF (.NOT.ALLOCATED(CSV)) ALLOCATE(CSV       (KSCAL)     )! name of all scalar variables
IF (.NOT.ALLOCATED(XDIR_ALB)) ALLOCATE(XDIR_ALB  (KNI,KBANDS)) ! direct albedo for each band
IF (.NOT.ALLOCATED(XSCA_ALB)) ALLOCATE(XSCA_ALB  (KNI,KBANDS)) ! diffuse albedo for each band
IF (.NOT.ALLOCATED(XEMIS)) ALLOCATE(XEMIS     (KNI)       )! emissivity
IF (.NOT.ALLOCATED(XTSRAD)) ALLOCATE(XTSRAD    (KNI)       )! radiative temperature
IF (.NOT.ALLOCATED(XTSUN)) ALLOCATE(XTSUN     (KNI)       )! solar time                    (s from midnight)
IF (.NOT.ALLOCATED(XZREF)) ALLOCATE(XZREF     (KNI)       )! height of T,q forcing                 (m)
IF (.NOT.ALLOCATED(XUREF)) ALLOCATE(XUREF     (KNI)       )! height of wind forcing                (m)
IF (.NOT.ALLOCATED(XTA)) ALLOCATE(XTA       (KNI)       )! air temperature forcing               (K)
IF (.NOT.ALLOCATED(XQA)) ALLOCATE(XQA       (KNI)       )! air specific humidity forcing         (kg/m3)
IF (.NOT.ALLOCATED(XZS)) ALLOCATE(XZS       (KNI)       )! orography                             (m)
IF (.NOT.ALLOCATED(XRHOA)) ALLOCATE(XRHOA     (KNI)       )! air density forcing                   (kg/m3)
IF (.NOT.ALLOCATED(XSV)) ALLOCATE(XSV       (KNI,KSCAL) ) ! scalar variables
IF (.NOT.ALLOCATED(XU)) ALLOCATE(XU        (KNI)       )! zonal wind                            (m/s)
IF (.NOT.ALLOCATED(XV)) ALLOCATE(XV        (KNI)       )! meridian wind                         (m/s)
IF (.NOT.ALLOCATED(XDIR_SW)) ALLOCATE(XDIR_SW   (KNI,KBANDS))    ! direct  solar radiation (on horizontal surf.)
IF (.NOT.ALLOCATED(XSCA_SW)) ALLOCATE(XSCA_SW   (KNI,KBANDS))    ! diffuse solar radiation (on horizontal surf.)
IF (.NOT.ALLOCATED(XSW_BANDS)) ALLOCATE(XSW_BANDS (KBANDS)    ) ! mean wavelength of each shortwave band (m)
IF (.NOT.ALLOCATED(XZENITH)) ALLOCATE(XZENITH   (KNI)       )! zenithal angle       (radian from the vertical)
IF (.NOT.ALLOCATED(XZENITH2)) ALLOCATE(XZENITH2  (KNI)       )! zenithal angle at t+1(radian from the vertical)
IF (.NOT.ALLOCATED(XAZIM)) ALLOCATE(XAZIM     (KNI)       )! azimuthal angle      (radian from North, clockwise)
IF (.NOT.ALLOCATED(XLW)) ALLOCATE(XLW       (KNI)       )! longwave radiation (on horizontal surf.)
IF (.NOT.ALLOCATED(XPS)) ALLOCATE(XPS       (KNI)       )! pressure at atmospheric model surface (Pa)
IF (.NOT.ALLOCATED(XPA)) ALLOCATE(XPA       (KNI)       )! pressure at forcing level             (Pa)
IF (.NOT.ALLOCATED(XCO2)) ALLOCATE(XCO2      (KNI)       )! CO2 concentration in the air          (kg/kg)
IF (.NOT.ALLOCATED(XSNOW)) ALLOCATE(XSNOW     (KNI)       )! snow precipitation                    (kg/m2/s)
IF (.NOT.ALLOCATED(XRAIN)) ALLOCATE(XRAIN     (KNI)       )! liquid precipitation                  (kg/m2/s)
IF (.NOT.ALLOCATED(XSFTH)) ALLOCATE(XSFTH     (KNI)       )! flux of heat                          (W/m2)
IF (.NOT.ALLOCATED(XSFTQ)) ALLOCATE(XSFTQ     (KNI)       )! flux of water vapor                   (kg/m2/s)
IF (.NOT.ALLOCATED(XSFU)) ALLOCATE(XSFU      (KNI)       )! zonal momentum flux                   (pa)
IF (.NOT.ALLOCATED(XSFV)) ALLOCATE(XSFV      (KNI)       )! meridian momentum flux                (pa)
IF (.NOT.ALLOCATED(XSFCO2)) ALLOCATE(XSFCO2    (KNI)       )! flux of CO2                           (kg/m2/s)
IF (.NOT.ALLOCATED(XSFTS)) ALLOCATE(XSFTS     (KNI,KSCAL) ) ! flux of scalar var.                  (kg/m2/s)
IF (.NOT.ALLOCATED(XPEW_A_COEF)) ALLOCATE(XPEW_A_COEF (KNI)     )! implicit coefficients
IF (.NOT.ALLOCATED(XPEW_B_COEF)) ALLOCATE(XPEW_B_COEF (KNI)     )! needed if HCOUPLING='I'
IF (.NOT.ALLOCATED(XPET_A_COEF)) ALLOCATE(XPET_A_COEF (KNI)     )
IF (.NOT.ALLOCATED(XPEQ_A_COEF)) ALLOCATE(XPEQ_A_COEF (KNI)     )
IF (.NOT.ALLOCATED(XPET_B_COEF)) ALLOCATE(XPET_B_COEF (KNI)     )
IF (.NOT.ALLOCATED(XPEQ_B_COEF)) ALLOCATE(XPEQ_B_COEF (KNI)     )
IF (.NOT.ALLOCATED(XTSURF)) ALLOCATE(XTSURF(KNI)     )
IF (.NOT.ALLOCATED(XZ0)   ) ALLOCATE(XZ0   (KNI)     )
IF (.NOT.ALLOCATED(XZ0H)  ) ALLOCATE(XZ0H  (KNI)     )
IF (.NOT.ALLOCATED(XQSURF)) ALLOCATE(XQSURF(KNI)     )
!
IF (SIZE(CSV)>=1) CSV(1) = '#CO   '
IF (SIZE(CSV)>=2) CSV(2) = '#O3   '     
IF (SIZE(CSV)>=3) CSV(3) = '#H2O2 ' 
IF (SIZE(CSV)>=4) CSV(4) = '#NO   '  
IF (SIZE(CSV)>=5) CSV(5) = '#NO2  ' 
IF (SIZE(CSV)>=6) CSV(6) = '#NO3  '
IF (SIZE(CSV)>=7) CSV(7) = '#N2O5 '
IF (SIZE(CSV)>=8) CSV(8) = '#HONO '
IF (SIZE(CSV)>=9) CSV(9) = '#HNO3 '
IF (SIZE(CSV)>=10) CSV(10) = '#HNO4 '
IF (SIZE(CSV)>=11) CSV(11) = '#NH3  '
IF (SIZE(CSV)>=12) CSV(12) = '#SO2  '
IF (SIZE(CSV)>=13) CSV(13) = '#SULF '
IF (SIZE(CSV)>=14) CSV(14) = '#OH   '
IF (SIZE(CSV)>=15) CSV(15) = '#HO2  '
IF (SIZE(CSV)>=16) CSV(16) = '#CH4  '
IF (SIZE(CSV)>=17) CSV(17) = '#ETH  '
IF (SIZE(CSV)>=18) CSV(18) = '#ALKA '
IF (SIZE(CSV)>=19) CSV(19) = '#ALKE '
IF (SIZE(CSV)>=20) CSV(20) = '#BIO  '
IF (SIZE(CSV)>=21) CSV(21) = '#ARO  '
IF (SIZE(CSV)>=22) CSV(22) = '#HCHO '
IF (SIZE(CSV)>=23) CSV(23) = '#ALD  '
IF (SIZE(CSV)>=24) CSV(24) = '#KET  '
IF (SIZE(CSV)>=25) CSV(25) = '#CARBO'
IF (SIZE(CSV)>=26) CSV(26) = '#ONIT '
IF (SIZE(CSV)>=27) CSV(27) = '#PAN  '
IF (SIZE(CSV)>=28) CSV(28) = '#OP1  '
IF (SIZE(CSV)>=29) CSV(29) = '#OP2  '
IF (SIZE(CSV)>=30) CSV(30) = '#ORA  '
IF (SIZE(CSV)>=31) CSV(31) = '#ORA2 '
IF (SIZE(CSV)>=32) CSV(32) = '#MO2  '
IF (SIZE(CSV)>=33) CSV(33) = '#ALKAP'
IF (SIZE(CSV)>=34) CSV(34) = '#ALKEP'
IF (SIZE(CSV)>=35) CSV(35) = '#BIOP '
IF (SIZE(CSV)>=36) CSV(36) = '#PHO  '
IF (SIZE(CSV)>=37) CSV(37) = '#ADD  '
IF (SIZE(CSV)>=38) CSV(38) = '#AROP '
IF (SIZE(CSV)>=39) CSV(39) = '#CARBO'
IF (SIZE(CSV)>=40) CSV(40) = '#OLN  '
IF (SIZE(CSV)>=41) CSV(41) = '#XO2  '
IF (SIZE(CSV)>=42) CSV(42) = '@M0I '
IF (SIZE(CSV)>=43) CSV(43) = '@M0J '
IF (SIZE(CSV)>=44) CSV(44) = '@M6I '
IF (SIZE(CSV)>=45) CSV(45) = '@M6J '
IF (SIZE(CSV)>=46) CSV(46) = '@H2OI'
IF (SIZE(CSV)>=47) CSV(47) = '@H2OJ'
IF (SIZE(CSV)>=48) CSV(48) = '@SO4I'
IF (SIZE(CSV)>=49) CSV(49) = '@SO4J'
IF (SIZE(CSV)>=50) CSV(50) = '@NO3I'
IF (SIZE(CSV)>=51) CSV(51) = '@NO3J'
IF (SIZE(CSV)>=52) CSV(52) = '@NH3I'
IF (SIZE(CSV)>=53) CSV(53) = '@NH3J'
IF (SIZE(CSV)>=54) CSV(54) = '@OCI'
IF (SIZE(CSV)>=55) CSV(55) = '@OCJ'
IF (SIZE(CSV)>=56) CSV(56) = '@BCI'
IF (SIZE(CSV)>=57) CSV(57) = '@BCJ'
IF (SIZE(CSV)>=58) CSV(58) = '@DSTI'
IF (SIZE(CSV)>=59) CSV(59) = '@DSTJ'
!
!CSV       (:)  ='UNDEF '! name of all scalar variables
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
XTSURF    (:)=XUNDEF ! effective temperature                  (K)
XZ0       (:)=XUNDEF ! surface roughness length for momentum  (m)
XZ0H      (:)=XUNDEF ! surface roughness length for heat      (m)
XQSURF    (:)=XUNDEF ! specific humidity at surface           (kg/kg)
!
IF (LHOOK) CALL DR_HOOK('OL_ALLOC_ATM',1,ZHOOK_HANDLE)

END SUBROUTINE OL_ALLOC_ATM
