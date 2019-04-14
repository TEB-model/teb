!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ######spl
      MODULE MODE_THERMOS
!     ####################
!
!!****  *MODE_THERMO* -
!!
!!    PURPOSE
!!    -------
!      
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!       NONE          
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    28/08/94 
!!      Modified    01/2006 : sea flux parameterization.
!--------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!-------------------------------------------------------------------------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
INTERFACE PSAT
  MODULE PROCEDURE PSAT_0D
  MODULE PROCEDURE PSAT_1D
END INTERFACE
INTERFACE QSAT
  MODULE PROCEDURE QSATW_0D         
  MODULE PROCEDURE QSATW_1D
  MODULE PROCEDURE QSATW_2D
END INTERFACE
INTERFACE QSAT_SEAWATER
  MODULE PROCEDURE QSATSEAW_1D
END INTERFACE
INTERFACE DQSAT
  MODULE PROCEDURE DQSATW_O_DT_1D
END INTERFACE
INTERFACE QSATI
  MODULE PROCEDURE QSATI_1D
  MODULE PROCEDURE QSATI_2D
END INTERFACE
INTERFACE DQSATI
  MODULE PROCEDURE DQSATI_O_DT_1D
END INTERFACE
CONTAINS
!-------------------------------------------------------------------------------
!     ######################################
      FUNCTION PSAT_0D(PT) RESULT(PPSAT)
!     ######################################
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CSTS
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, INTENT(IN)                :: PT     ! Temperature (Kelvin)
REAL                            :: PPSAT  ! saturation vapor 
                                          ! specific humidity
                                          ! with respect to
                                          ! water (kg/kg)
REAL(KIND=JPRB) :: ZHOOK_HANDLE                                          
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_THERMOS:PSAT_0D',0,ZHOOK_HANDLE)
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
PPSAT = EXP( XALPW - XBETAW/PT - XGAMW*LOG(PT)  )
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_THERMOS:PSAT_0D',1,ZHOOK_HANDLE)
!
END FUNCTION PSAT_0D
!-------------------------------------------------------------------------------
!     ######################################
      FUNCTION PSAT_1D(PT) RESULT(PPSAT)
!     ######################################
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CSTS
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, DIMENSION(:), INTENT(IN)                :: PT     ! Temperature (Kelvin)
REAL, DIMENSION(SIZE(PT))                     :: PPSAT  ! saturation vapor pressure (Pa)
!
INTEGER                         :: JJ !loop index
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_THERMOS:PSAT_1D',0,ZHOOK_HANDLE)
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
!cdir nodep
DO JJ=1,SIZE(PT)
  PPSAT(JJ) = EXP( XALPW - XBETAW/PT(JJ) - XGAMW*LOG(PT(JJ))  )
ENDDO
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_THERMOS:PSAT_1D',1,ZHOOK_HANDLE)
!
END FUNCTION PSAT_1D
!-------------------------------------------------------------------------------
!     ######################################
      FUNCTION QSATW_0D(PT,PP) RESULT(PQSAT)
!     ######################################
!
!!****  *QSATW * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor 
!     pressure from temperature 
!      
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor 
!!    pressure of the triple point es(Tt) (XESTT), i.e  
!!     
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!  
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt) 
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!  
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function  
!!      
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH 
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CSTS
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, INTENT(IN)                :: PT     ! Temperature (Kelvin)
REAL, INTENT(IN)                :: PP     ! Pressure (Pa)
REAL                            :: PQSAT  ! saturation vapor 
                                                        ! specific humidity
                                                        ! with respect to
                                                        ! water (kg/kg)
!
!*       0.2   Declarations of local variables
!
REAL                           :: ZFOES  ! saturation vapor 
                                                        ! pressure
                                                        ! (Pascal) 
!
REAL                           :: ZWORK1
REAL                           :: ZWORK2
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_THERMOS:QSATW_0D',0,ZHOOK_HANDLE)
!
ZWORK2 = XRD/XRV
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
ZFOES = EXP( XALPW - XBETAW/PT - XGAMW*LOG(PT)  )
ZWORK1    = ZFOES/PP
!
!*       2.    COMPUTE SATURATION HUMIDITY
!              ---------------------------
!
PQSAT = ZWORK2*ZWORK1 / (1.+(ZWORK2-1.)*ZWORK1)
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_THERMOS:QSATW_0D',1,ZHOOK_HANDLE)
!
END FUNCTION QSATW_0D
!-------------------------------------------------------------------------------
!
!     ######################################
      FUNCTION QSATW_1D(PT,PP) RESULT(PQSAT)
!     ######################################
!
!!****  *QSATW * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor 
!     pressure from temperature 
!      
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor 
!!    pressure of the triple point es(Tt) (XESTT), i.e  
!!     
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!  
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt) 
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!  
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function  
!!      
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH 
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CSTS
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, DIMENSION(:), INTENT(IN)                :: PT     ! Temperature
                                                        ! (Kelvin)
REAL, DIMENSION(:), INTENT(IN)                :: PP     ! Pressure
                                                        ! (Pa)
REAL, DIMENSION(SIZE(PT))                   :: PQSAT  ! saturation vapor 
                                                        ! specific humidity
                                                        ! with respect to
                                                        ! water (kg/kg)
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT))                   :: ZFOES  ! saturation vapor 
                                                        ! pressure
                                                        ! (Pascal) 
!
REAL, DIMENSION(SIZE(PT))                   :: ZWORK1
REAL                                        :: ZWORK2 
INTEGER                         :: JJ !loop index
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMOS:QSATW_1D',0,ZHOOK_HANDLE)
ZWORK2 = XRD/XRV
!
!cdir nodep
DO JJ=1,SIZE(PT)
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
  ZFOES(JJ) = EXP( XALPW - XBETAW/PT(JJ) - XGAMW*LOG(PT(JJ))  )
  ZWORK1(JJ)    = ZFOES(JJ)/PP(JJ)
!
!*       2.    COMPUTE SATURATION HUMIDITY
!              ---------------------------
!
  PQSAT(JJ) = ZWORK2*ZWORK1(JJ) / (1.+(ZWORK2-1.)*ZWORK1(JJ))
!
ENDDO
IF (LHOOK) CALL DR_HOOK('MODE_THERMOS:QSATW_1D',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END FUNCTION QSATW_1D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!     ######################################
      FUNCTION QSATW_2D(PT,PP,KMASK,KL) RESULT(PQSAT)
!     ######################################
!
!!****  *QSATW * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor 
!     pressure from temperature 
!      
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor 
!!    pressure of the triple point es(Tt) (XESTT), i.e  
!!     
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!  
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt) 
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!  
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function  
!!      
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH 
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR,  ONLY : XUNDEF
USE MODD_CSTS
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, DIMENSION(:,:), INTENT(IN)              :: PT     ! Temperature
                                                        ! (Kelvin)
REAL, DIMENSION(:,:), INTENT(IN)              :: PP     ! Pressure
                                                        ! (Pa)
!                                                        
INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL   :: KMASK
!                                                KMASK = Number of soil moisture layers (DIF option)
INTEGER,               INTENT(IN), OPTIONAL   :: KL
!                                                KL = Max number of soil moisture layers (DIF option)
!
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2))        :: PQSAT  ! saturation vapor 
                                                        ! specific humidity
                                                        ! with respect to
                                                        ! water (kg/kg)
!
!*       0.2   Declarations of local variables
!
INTEGER, DIMENSION(SIZE(PT,1)) :: IMASK
!
REAL            :: ZFOES  ! saturation vapor pressure (Pascal) 
INTEGER         :: JJ, JI, INL, IWORK   ! loop indexes
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMOS:QSATW_2D',0,ZHOOK_HANDLE)
!
IF(PRESENT(KMASK).AND.PRESENT(KL))THEN
  IMASK(:)=KMASK(:)
  INL=KL
ELSE
  IMASK(:)=SIZE(PT,2)  
  INL=SIZE(PT,2)
ENDIF
!
PQSAT(:,:)=XUNDEF
!
DO JJ=1,INL
  DO JI=1,SIZE(PT,1)
!
     IWORK=IMASK(JI)
     IF(JJ<=IWORK)THEN
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
         ZFOES = EXP( XALPW - XBETAW/PT(JI,JJ) - XGAMW*LOG(PT(JI,JJ))  )
!
!*       2.    COMPUTE SATURATION HUMIDITY
!              ---------------------------
!
         PQSAT(JI,JJ) = XRD/XRV*ZFOES/PP(JI,JJ) / (1.+(XRD/XRV-1.)*ZFOES/PP(JI,JJ))  
!
     ENDIF
!
  ENDDO
ENDDO
IF (LHOOK) CALL DR_HOOK('MODE_THERMOS:QSATW_2D',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END FUNCTION QSATW_2D
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
!     ######################################
      FUNCTION QSATSEAW_1D(PT,PP) RESULT(PQSAT)
!     ######################################
!
!!****  *QSATW * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor 
!     pressure from temperature over saline seawater
!      
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor 
!!    pressure of the triple point es(Tt) (XESTT), i.e  
!!    The reduction due to salinity is compute with the factor 0.98 (reduction of 2%)
!!     
!!         es(T)= 0.98*EXP( alphaw - betaw /T - gammaw Log(T) )
!!  
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt) 
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!  
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function  
!!      
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH 
!!      Zeng, X., Zhao, M., and Dickinson, R. E., 1998 : Intercomparaison of bulk
!!      aerodynamic algorithm for the computation of sea surface fluxes using
!!      TOGA COARE and TAO data. Journal of Climate, vol 11, n°10, pp 2628--2644
!!
!!
!!    AUTHOR
!!    ------
!!	C. Lebeaupin    * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    6/04/2005 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CSTS
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, DIMENSION(:), INTENT(IN)                :: PT     ! Temperature
                                                        ! (Kelvin)
REAL, DIMENSION(:), INTENT(IN)                :: PP     ! Pressure
                                                        ! (Pa)
REAL, DIMENSION(SIZE(PT))                   :: PQSAT  ! saturation vapor 
                                                        ! specific humidity
                                                        ! with respect to
                                                        ! water (kg/kg)
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT))                   :: ZFOES  ! saturation vapor 
                                                        ! pressure
                                                        ! (Pascal) 
!
INTEGER                         :: JJ   ! loop index
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMOS:QSATSEAW_1D',0,ZHOOK_HANDLE)
DO JJ = 1, SIZE(PT)
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
  ZFOES(JJ) = 0.98*EXP( XALPW - XBETAW/PT(JJ) - XGAMW*LOG(PT(JJ))  )
! vapor pressure reduction of 2% over saline seawater could have a significant 
! impact on the computation of surface latent heat flux under strong wind 
! conditions (Zeng et al, 1998). 
!
!*       2.    COMPUTE SATURATION HUMIDITY
!              ---------------------------
!
  PQSAT(JJ) = XRD/XRV*ZFOES(JJ)/PP(JJ)   &
                     / (1.+(XRD/XRV-1.)*ZFOES(JJ)/PP(JJ))  
!
ENDDO
IF (LHOOK) CALL DR_HOOK('MODE_THERMOS:QSATSEAW_1D',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END FUNCTION QSATSEAW_1D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     ##############################################################
      FUNCTION DQSATW_O_DT_1D(PT,PP,PQSAT) RESULT(PDQSAT)
!     ##############################################################
!
!!****  *QSATW * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor 
!     pressure from temperature 
!      
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor 
!!    pressure of the triple point es(Tt) (XESTT), i.e  
!!     
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!  
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt) 
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!
!!      Finally, dqsat / dT  (T) is computed.
!!  
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function  
!!      
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH 
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CSTS       
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, DIMENSION(:),  INTENT(IN)             :: PT     ! Temperature
                                                          ! (Kelvin)
REAL, DIMENSION(:),  INTENT(IN)             :: PP     ! Pressure
                                                          ! (Pa)
REAL, DIMENSION(:),  INTENT(IN)             :: PQSAT  ! saturation vapor 
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg))
REAL, DIMENSION(SIZE(PT))                   :: PDQSAT ! derivative according
                                                          ! to temperature of
                                                          ! saturation vapor 
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg))
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT))  :: ZFOES  ! saturation vapor 
                                                          ! pressure
                                                          ! (Pascal) 
!
REAL                       :: ZWORK1
REAL, DIMENSION(SIZE(PT))  :: ZWORK2   ! loop index
INTEGER :: JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMOS:DQSATW_O_DT_1D',0,ZHOOK_HANDLE)
ZWORK1=XRD/XRV
!
DO JJ=1,SIZE(PT)
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
  ZFOES(JJ) = PP(JJ) / (1.+ZWORK1*(1./PQSAT(JJ)-1.))
  ZWORK2(JJ) = ZFOES(JJ) / PP(JJ)
!
!*       2.    DERIVATION ACCORDING TO TEMPERATURE
!              -----------------------------------
!
  PDQSAT(JJ) = PQSAT(JJ) / (1.+(ZWORK1-1.)*ZWORK2(JJ) ) &
                     * (XBETAW/PT(JJ)**2 - XGAMW/PT(JJ))  
!
ENDDO
IF (LHOOK) CALL DR_HOOK('MODE_THERMOS:DQSATW_O_DT_1D',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END FUNCTION DQSATW_O_DT_1D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     ##############################################################
      FUNCTION DQSATI_O_DT_1D(PT,PP,PQSAT) RESULT(PDQSAT)
!     ##############################################################
!
!!****  *QSATW * - function to compute saturation vapor humidity from
!!                 temperature (with respect to ice)
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor 
!     pressure from temperature 
!      
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor 
!!    pressure of the triple point es(Tt) (XESTT), i.e  
!!     
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!  
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt) 
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!
!!      Finally, dqsat / dT  (T) is computed.
!!  
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function  
!!      
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH 
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CSTS       
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL,    DIMENSION(:), INTENT(IN)               :: PT     ! Temperature
                                                          ! (Kelvin)
REAL,    DIMENSION(:), INTENT(IN)               :: PP     ! Pressure
                                                          ! (Pa)
REAL,    DIMENSION(:), INTENT(IN)               :: PQSAT  ! saturation vapor 
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg))
REAL,    DIMENSION(SIZE(PT))                    :: PDQSAT ! derivative according
                                                          ! to temperature of
                                                          ! saturation vapor 
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg))
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT))                       :: ZFOES  ! saturation vapor 
                                                          ! pressure
                                                          ! (Pascal) 
!
REAL   :: ZWORK1
INTEGER                 ::   JJ  ! loop index
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMOS:DQSATI_O_DT_1D',0,ZHOOK_HANDLE)
ZWORK1=XRD/XRV
DO JJ = 1,SIZE(PT)
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
  ZFOES(JJ) = PP(JJ) / (1.+ZWORK1*(1./PQSAT(JJ)-1.))
!
!*       3.    DERIVATION ACCORDING TO TEMPERATURE
!              -----------------------------------
!
  PDQSAT(JJ) = PQSAT(JJ) / (1.+(ZWORK1-1.)*ZFOES(JJ)/PP(JJ) ) &
                     * (XBETAI/PT(JJ)**2 - XGAMI/PT(JJ))  
!
ENDDO
IF (LHOOK) CALL DR_HOOK('MODE_THERMOS:DQSATI_O_DT_1D',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END FUNCTION DQSATI_O_DT_1D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!     ######################################
      FUNCTION QSATI_1D(PT,PP) RESULT(PQSAT)
!     ######################################
!
!!****  *QSATI * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor 
!     pressure from temperature 
!      
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor 
!!    pressure of the triple point es(Tt) (XESTT), i.e  
!!     
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!  
!!     with :
!!       alphaw (XALPI) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt) 
!!       betaw (XBETAI) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMI) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!  
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPI   : Constant for saturation vapor pressure function
!!        XBETAI  : Constant for saturation vapor pressure function
!!        XGAMI   : Constant for saturation vapor pressure function  
!!      
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH 
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CSTS       
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, DIMENSION(:), INTENT(IN)                :: PT     ! Temperature
                                                        ! (Kelvin)
REAL, DIMENSION(:), INTENT(IN)                :: PP     ! Pressure
                                                        ! (Pa)
REAL, DIMENSION(SIZE(PT))                   :: PQSAT  ! saturation vapor 
                                                        ! specific humidity
                                                        ! with respect to
                                                        ! water (kg/kg)
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT))                   :: ZFOES  ! saturation vapor 
                                                        ! pressure
                                                        ! (Pascal) 
!
INTEGER   :: JJ !loop index
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMOS:QSATI_1D',0,ZHOOK_HANDLE)
DO JJ = 1, SIZE(PT)
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
  ZFOES(JJ) = EXP( XALPI - XBETAI/PT(JJ) - XGAMI*LOG(PT(JJ))  )
!
!*       2.    COMPUTE SATURATION HUMIDITY
!              ---------------------------
!
  PQSAT(JJ) = XRD/XRV*ZFOES(JJ)/PP(JJ)   &
                     / (1.+(XRD/XRV-1.)*ZFOES(JJ)/PP(JJ))  
!
ENDDO
IF (LHOOK) CALL DR_HOOK('MODE_THERMOS:QSATI_1D',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END FUNCTION QSATI_1D
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!     ######################################
      FUNCTION QSATI_2D(PT,PP,KMASK,KL) RESULT(PQSAT)
!     ######################################
!
!!****  *QSATI * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor 
!     pressure from temperature 
!      
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor 
!!    pressure of the triple point es(Tt) (XESTT), i.e  
!!     
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!  
!!     with :
!!       alphaw (XALPI) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt) 
!!       betaw (XBETAI) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMI) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!  
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPI   : Constant for saturation vapor pressure function
!!        XBETAI  : Constant for saturation vapor pressure function
!!        XGAMI   : Constant for saturation vapor pressure function  
!!      
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH 
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR,  ONLY : XUNDEF
USE MODD_CSTS       
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, DIMENSION(:,:), INTENT(IN)            :: PT     ! Temperature
                                                      ! (Kelvin)
REAL, DIMENSION(:,:), INTENT(IN)            :: PP     ! Pressure
                                                      ! (Pa)
!                                                        
INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL   :: KMASK
!                                                KMASK = Number of soil moisture layers (DIF option)
INTEGER,               INTENT(IN), OPTIONAL   :: KL
!                                                KL = Max number of soil moisture layers (DIF option)
!                                                      
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2))      :: PQSAT  ! saturation vapor 
                                                      ! specific humidity
                                                      ! with respect to
                                                      ! water (kg/kg)
!
!*       0.2   Declarations of local variables
!
REAL            :: ZFOES  ! saturation vapor pressure (Pascal) 
!
INTEGER, DIMENSION(SIZE(PT,1)) :: IMASK
!
INTEGER         :: JJ, JI, INL, IWORK   ! loop indexes
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMOS:QSATI_2D',0,ZHOOK_HANDLE)
!
IF(PRESENT(KMASK))THEN
  IMASK(:)=KMASK(:)    
  INL=KL
ELSE
  IMASK(:)=SIZE(PT,2)  
  INL=SIZE(PT,2)
ENDIF
!
PQSAT(:,:)=XUNDEF
!
DO JJ=1,INL
  DO JI=1,SIZE(PT,1)
!
     IWORK=IMASK(JI)
     IF(JJ<=IWORK)THEN
!  
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
         ZFOES = EXP( XALPI - XBETAI/PT(JI,JJ) - XGAMI*LOG(PT(JI,JJ))  )
!
!*       2.    COMPUTE SATURATION HUMIDITY
!              ---------------------------
!
         PQSAT(JI,JJ) = XRD/XRV*ZFOES/PP(JI,JJ) / (1.+(XRD/XRV-1.)*ZFOES/PP(JI,JJ))  
!
     ENDIF
!
  ENDDO
ENDDO
IF (LHOOK) CALL DR_HOOK('MODE_THERMOS:QSATI_2D',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END FUNCTION QSATI_2D
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END MODULE MODE_THERMOS
