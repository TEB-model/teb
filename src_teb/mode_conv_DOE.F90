!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!###################
MODULE MODE_CONV_DOE
!###################
!
!!****  *MODE_CONV_DOE* -
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
!!    G. Pigeon
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/08/12 
!
!--------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!--------------------------------------------------------------------------------
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
INTERFACE CHTC_VERT_DOE
        MODULE PROCEDURE CHTC_VERT_DOE_0D
        MODULE PROCEDURE CHTC_VERT_DOE_1D
END INTERFACE
INTERFACE CHTC_UP_DOE
        MODULE PROCEDURE CHTC_UP_DOE_0D
        MODULE PROCEDURE CHTC_UP_DOE_1D
END INTERFACE
INTERFACE CHTC_SMOOTH_LEE_DOE
        MODULE PROCEDURE CHTC_SMOOTH_LEE_DOE_0D
        MODULE PROCEDURE CHTC_SMOOTH_LEE_DOE_1D
END INTERFACE
INTERFACE CHTC_SMOOTH_WIND_DOE
        MODULE PROCEDURE CHTC_SMOOTH_WIND_DOE_0D
        MODULE PROCEDURE CHTC_SMOOTH_WIND_DOE_1D
END INTERFACE
INTERFACE CHTC_ROUGH_DOE
        MODULE PROCEDURE CHTC_ROUGH_DOE_0D
        MODULE PROCEDURE CHTC_ROUGH_DOE_1D
END INTERFACE
!
CONTAINS
!----------------------------
!#############################################
FUNCTION CHTC_VERT_DOE_0D(PTS, PTA) RESULT(PCHTC)
!#############################################
!
!!****  *CHTC_VERT_DOE* - 
!!
!!    PURPOSE
!!    -------
!     function to compute convective surface coefficient for
!     a vertical surface from surface temperature and air
!     temperature
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    REFERENCE
!!    ---------
!!    EnergyPlus, Engineering Reference, DOE-2 model for convection on outside
!!    surfaces, p68
!!
!!    AUTHOR
!!    ------
!!      G. Pigeon       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/00/2012 
!
!-------------------------------------------------------------------------------
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, INTENT(IN)                :: PTS     ! Surface temperature (Kelvin)
REAL, INTENT(IN)                :: PTA     ! Air temperature (Kelvin)
REAL                            :: PCHTC   ! Convective heat transfer coefficient [W/(m2.K)]
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_VERT_DOE',0,ZHOOK_HANDLE)
!
!*       1.    COMPUTE THE CHTC
!              ----------------
!
PCHTC = 1.31 * (ABS(PTA - PTS))**(1./3.) 
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_VERT_DOE',1,ZHOOK_HANDLE)
!
END FUNCTION CHTC_VERT_DOE_0D
!#########################
!#############################################
FUNCTION CHTC_VERT_DOE_1D(PTS, PTA) RESULT(PCHTC)
!#############################################
!
!!****  *CHTC_VERT_DOE* - 
!!
!!    PURPOSE
!!    -------
!     function to compute convective surface coefficient for
!     a vertical surface from surface temperature and air
!     temperature
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    REFERENCE
!!    ---------
!!    EnergyPlus, Engineering Reference, DOE-2 model for convection on outside
!!    surfaces, p68
!!
!!    AUTHOR
!!    ------
!!      G. Pigeon       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/00/2012 
!
!-------------------------------------------------------------------------------
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, DIMENSION(:), INTENT(IN)                :: PTS     ! Surface temperature (Kelvin)
REAL, DIMENSION(:), INTENT(IN)                :: PTA     ! Air temperature (Kelvin)
REAL, DIMENSION(SIZE(PTS))                    :: PCHTC   ! Convective heat transfer coefficient [W/(m2.K)]
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_VERT_DOE',0,ZHOOK_HANDLE)
!
!*       1.    COMPUTE THE CHTC
!              ----------------
!
PCHTC = 1.31 * (ABS(PTA - PTS))**(1./3.) 
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_VERT_DOE',1,ZHOOK_HANDLE)
!
END FUNCTION CHTC_VERT_DOE_1D
!#########################
!
!#############################################
FUNCTION CHTC_UP_DOE_1D(PTS, PTA) RESULT(PCHTC)
!#############################################
!
!!****  *CHTC_UP_DOE* - 
!!
!!    PURPOSE
!!    -------
!     function to compute convective surface coefficient for
!     an upward surface from surface temperature and air
!     temperature
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    REFERENCE
!!    ---------
!!    EnergyPlus, Engineering Reference, DOE-2 model for convection on outside
!!    surfaces, p68
!!
!!    AUTHOR
!!    ------
!!      G. Pigeon       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/00/2012 
!
!-------------------------------------------------------------------------------
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
REAL, DIMENSION(:), INTENT(IN)                :: PTS     ! Surface temperature (Kelvin)
REAL, DIMENSION(:), INTENT(IN)                :: PTA     ! Air temperature (Kelvin)
REAL, DIMENSION(SIZE(PTS))                    :: PCHTC   ! Convective heat transfer coefficient [W/(m2.K)]
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_UP_DOE',0,ZHOOK_HANDLE)
!
!*       1.    COMPUTE THE CHTC
!              ----------------
!
WHERE(PTA < PTS)
   PCHTC  = 9.482 * (PTS - PTA)**(1./3.) / (7.283 - 1)
ELSE WHERE
   PCHTC  = 1.810 * (PTA - PTS)**(1./3.) / (1.382 + 1)  
END WHERE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_UP_DOE',1,ZHOOK_HANDLE)
!
END FUNCTION CHTC_UP_DOE_1D
!#######################
!#############################################
FUNCTION CHTC_UP_DOE_0D(PTS, PTA) RESULT(PCHTC)
!#############################################
!
!!****  *CHTC_UP_DOE* - 
!!
!!    PURPOSE
!!    -------
!     function to compute convective surface coefficient for
!     an upward surface from surface temperature and air
!     temperature
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    REFERENCE
!!    ---------
!!    EnergyPlus, Engineering Reference, DOE-2 model for convection on outside
!!    surfaces, p68
!!
!!    AUTHOR
!!    ------
!!      G. Pigeon       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/00/2012 
!
!-------------------------------------------------------------------------------
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
REAL, INTENT(IN)                :: PTS     ! Surface temperature (Kelvin)
REAL, INTENT(IN)                :: PTA     ! Air temperature (Kelvin)
REAL                            :: PCHTC   ! Convective heat transfer coefficient [W/(m2.K)]
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_UP_DOE',0,ZHOOK_HANDLE)
!
!*       1.    COMPUTE THE CHTC
!              ----------------
!
IF (PTA < PTS) THEN
   PCHTC  = 9.482 * (PTS - PTA)**(1./3.) / (7.283 - 1)
ELSE 
   PCHTC  = 1.810 * (PTA - PTS)**(1./3.) / (1.382 + 1)  
ENDIF
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_UP_DOE',1,ZHOOK_HANDLE)
!
END FUNCTION CHTC_UP_DOE_0D
!#######################
!
!#############################################
FUNCTION CHTC_DOWN_DOE(PTS, PTA) RESULT(PCHTC)
!#############################################
!
!!****  *CHTC_DOWN_DOE* - 
!!
!!    PURPOSE
!!    -------
!     function to compute convective surface coefficient for
!     an downward surface from surface temperature and air
!     temperature
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    REFERENCE
!!    ---------
!!    EnergyPlus, Engineering Reference, DOE-2 model for convection on outside
!!    surfaces, p68
!!
!!    AUTHOR
!!    ------
!!      G. Pigeon       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/00/2012 
!
!-------------------------------------------------------------------------------
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
REAL, DIMENSION(:), INTENT(IN)                :: PTS     ! Surface temperature (Kelvin)
REAL, DIMENSION(:), INTENT(IN)                :: PTA     ! Air temperature (Kelvin)
REAL, DIMENSION(SIZE(PTS))                    :: PCHTC   ! Convective heat transfer coefficient [W/(m2.K)]
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_DOWN_DOE',0,ZHOOK_HANDLE)
!
!*       1.    COMPUTE THE CHTC
!              ----------------
!
WHERE(PTA > PTS)
   PCHTC  = 9.482 * (PTA - PTS)**(1./3.) / (7.283 - 1)
ELSE WHERE
   PCHTC  = 1.810 * (PTS - PTA)**(1./3.) / (1.382 + 1)  
END WHERE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_DOWN_DOE',1,ZHOOK_HANDLE)
!
END FUNCTION CHTC_DOWN_DOE
!#########################
!
!#############################################
FUNCTION CHTC_SMOOTH_LEE_DOE_0D(PCHTCN, PVMOD) RESULT(PCHTC)
!#############################################
!
!!****  *CHTC_SMOOTH_LEE_DOE* - 
!!
!!    PURPOSE
!!    -------
!     function to compute convective surface coefficient for
!     a leeward smooth surface from the natural convection coef and the
!     wind speed
!
!!**  METHOD
!!    ------
!!
!!    from EnergyPlus Engineering Reference, average the leeward/windward coef 
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    REFERENCE
!!    ---------
!!    EnergyPlus, Engineering Reference, DOE-2 model for convection on outside
!!    surfaces, p68
!!
!!    AUTHOR
!!    ------
!!      G. Pigeon       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/00/2012 
!
!-------------------------------------------------------------------------------
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
REAL,INTENT(IN)                :: PCHTCN  ! Convective heat transfer coefficient for natural conv. [W/(m2.K)]
REAL,INTENT(IN)                :: PVMOD   ! wind speed (m/s)
REAL                           :: PCHTC   ! Convective heat transfer coefficient [W/(m2.K)]
!
!*       0.2   Declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
REAL :: ZALEE, ZBLEE ! coef for leeward facade
!----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_SMOOTH_LEE_DOE',0,ZHOOK_HANDLE)
!
!*       1.    INITIALIZE COEF
!              ----------------
!
ZALEE = 2.86
ZBLEE = 0.617

!*       2.    COMPUTE THE CHTC
!              ----------------
!
PCHTC = SQRT(PCHTCN**2+(ZALEE*PVMOD**ZBLEE)**2)
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_SMOOTH_LEE_DOE',1,ZHOOK_HANDLE)
!
END FUNCTION CHTC_SMOOTH_LEE_DOE_0D
!#########################

!#############################################
FUNCTION CHTC_SMOOTH_LEE_DOE_1D(PCHTCN, PVMOD) RESULT(PCHTC)
!#############################################
!
!!****  *CHTC_SMOOTH_LEE_DOE* - 
!!
!!    PURPOSE
!!    -------
!     function to compute convective surface coefficient for
!     a leeward smooth surface from the natural convection coef and the
!     wind speed
!
!!**  METHOD
!!    ------
!!
!!    from EnergyPlus Engineering Reference, average the leeward/windward coef 
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    REFERENCE
!!    ---------
!!    EnergyPlus, Engineering Reference, DOE-2 model for convection on outside
!!    surfaces, p68
!!
!!    AUTHOR
!!    ------
!!      G. Pigeon       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/00/2012 
!
!-------------------------------------------------------------------------------
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
REAL, DIMENSION(:), INTENT(IN)                :: PCHTCN  ! Convective heat transfer coefficient for natural conv. [W/(m2.K)]
REAL, DIMENSION(:), INTENT(IN)                :: PVMOD   ! wind speed (m/s)
REAL, DIMENSION(SIZE(PCHTCN))                 :: PCHTC   ! Convective heat transfer coefficient [W/(m2.K)]
!
!*       0.2   Declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
REAL :: ZALEE, ZBLEE ! coef for windward/leeward facade
!----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_SMOOTH_LEE_DOE',0,ZHOOK_HANDLE)
!
!*       1.    INITIALIZE COEF
!              ----------------
!
ZALEE = 2.86
ZBLEE = 0.617

!*       2.    COMPUTE THE CHTC
!              ----------------
!
PCHTC = SQRT(PCHTCN**2+(ZALEE*PVMOD**ZBLEE)**2)
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_SMOOTH_LEE_DOE',1,ZHOOK_HANDLE)
!
END FUNCTION CHTC_SMOOTH_LEE_DOE_1D
!#########################
!#############################################
FUNCTION CHTC_SMOOTH_WIND_DOE_0D(PCHTCN, PVMOD) RESULT(PCHTC)
!#############################################
!
!!****  *CHTC_SMOOTH_WIND_DOE* - 
!!
!!    PURPOSE
!!    -------
!     function to compute convective surface coefficient for
!     a smooth surface from the natural convection coef and the
!     wind speed
!
!!**  METHOD
!!    ------
!!
!!    from EnergyPlus Engineering Reference, average the leeward/windward coef 
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    REFERENCE
!!    ---------
!!    EnergyPlus, Engineering Reference, DOE-2 model for convection on outside
!!    surfaces, p68
!!
!!    AUTHOR
!!    ------
!!      G. Pigeon       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/00/2012 
!
!-------------------------------------------------------------------------------
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
REAL,INTENT(IN)                :: PCHTCN  ! Convective heat transfer coefficient for natural conv. [W/(m2.K)]
REAL,INTENT(IN)                :: PVMOD   ! wind speed (m/s)
REAL                           :: PCHTC   ! Convective heat transfer coefficient [W/(m2.K)]
!
!*       0.2   Declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
REAL :: ZAWIND, ZBWIND ! coef for windward/leeward facade
!----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_SMOOTH_WIND_DOE',0,ZHOOK_HANDLE)
!
!*       1.    INITIALIZE COEF
!              ----------------
!
ZAWIND = 2.38
ZBWIND = 0.89

!*       2.    COMPUTE THE CHTC
!              ----------------
!
PCHTC = SQRT(PCHTCN**2+(ZAWIND*PVMOD**ZBWIND)**2)
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_SMOOTH_WIND_DOE',1,ZHOOK_HANDLE)
!
END FUNCTION CHTC_SMOOTH_WIND_DOE_0D
!#########################

!#############################################
FUNCTION CHTC_SMOOTH_WIND_DOE_1D(PCHTCN, PVMOD) RESULT(PCHTC)
!#############################################
!
!!****  *CHTC_SMOOTH_WIND_DOE* - 
!!
!!    PURPOSE
!!    -------
!     function to compute convective surface coefficient for
!     a windward smooth surface from the natural convection coef and the
!     wind speed
!
!!**  METHOD
!!    ------
!!
!!    from EnergyPlus Engineering Reference, average the leeward/windward coef 
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    REFERENCE
!!    ---------
!!    EnergyPlus, Engineering Reference, DOE-2 model for convection on outside
!!    surfaces, p68
!!
!!    AUTHOR
!!    ------
!!      G. Pigeon       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/00/2012 
!
!-------------------------------------------------------------------------------
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
REAL, DIMENSION(:), INTENT(IN)                :: PCHTCN  ! Convective heat transfer coefficient for natural conv. [W/(m2.K)]
REAL, DIMENSION(:), INTENT(IN)                :: PVMOD   ! wind speed (m/s)
REAL, DIMENSION(SIZE(PCHTCN))                 :: PCHTC   ! Convective heat transfer coefficient [W/(m2.K)]
!
!*       0.2   Declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
REAL :: ZAWIND, ZBWIND ! coef for windward facade
!----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_SMOOTH_WIND_DOE',0,ZHOOK_HANDLE)
!
!*       1.    INITIALIZE COEF
!              ----------------
!
ZAWIND = 2.38
ZBWIND = 0.89

!*       2.    COMPUTE THE CHTC
!              ----------------
!
PCHTC = SQRT(PCHTCN**2+(ZAWIND*PVMOD**ZBWIND)**2) 
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_SMOOTH_WIND_DOE',1,ZHOOK_HANDLE)
!
END FUNCTION CHTC_SMOOTH_WIND_DOE_1D
!#########################
!#############################################
FUNCTION CHTC_ROUGH_DOE_0D(PCHTCN, PCHTCS, PROUGH) RESULT(PCHTC)
!#############################################
!
!!****  *CHTC_ROUGH_DOE* - 
!!
!!    PURPOSE
!!    -------
!     function to compute convective surface coefficient for
!     a rough surface from the natural convection coef, the
!     smooth surface convective coef and the roughness coef
!
!!**  METHOD
!!    ------
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    REFERENCE
!!    ---------
!!    EnergyPlus, Engineering Reference, DOE-2 model for convection on outside
!!    surfaces, p68
!!
!!    AUTHOR
!!    ------
!!      G. Pigeon       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/00/2012 
!
!-------------------------------------------------------------------------------
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
REAL, INTENT(IN)                :: PCHTCN  ! Convective heat transfer coefficient for natural conv. [W/(m2.K)]
REAL, INTENT(IN)                :: PCHTCS  ! Convective heat transfer coefficient for a smooth surface [W/(m2.K)]
REAL, INTENT(IN)                :: PROUGH  ! Convective heat transfer coefficient for a smooth surface [W/(m2.K)]
REAL                            :: PCHTC   ! Convective heat transfer coefficient [W/(m2.K)]
!
!*       0.2   Declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_ROUGH_DOE',0,ZHOOK_HANDLE)
!
!*       1.    COMPUTE THE CHTC
!              ----------------
!
PCHTC = PCHTCN + PROUGH * (PCHTCS - PCHTCN)
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_ROUGH_DOE',1,ZHOOK_HANDLE)
!
END FUNCTION CHTC_ROUGH_DOE_0D
!#########################
!#############################################
FUNCTION CHTC_ROUGH_DOE_1D(PCHTCN, PCHTCS, PROUGH) RESULT(PCHTC)
!#############################################
!
!!****  *CHTC_ROUGH_DOE* - 
!!
!!    PURPOSE
!!    -------
!     function to compute convective surface coefficient for
!     a rough surface from the natural convection coef, the
!     smooth surface convective coef and the roughness coef
!
!!**  METHOD
!!    ------
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    REFERENCE
!!    ---------
!!    EnergyPlus, Engineering Reference, DOE-2 model for convection on outside
!!    surfaces, p68
!!
!!    AUTHOR
!!    ------
!!      G. Pigeon       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/00/2012 
!
!-------------------------------------------------------------------------------
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
REAL, DIMENSION(:), INTENT(IN)                :: PCHTCN  ! Convective heat transfer coefficient for natural conv. [W/(m2.K)]
REAL, DIMENSION(:), INTENT(IN)                :: PCHTCS  ! Convective heat transfer coefficient for a smooth surface [W/(m2.K)]
REAL, DIMENSION(:), INTENT(IN)                :: PROUGH  ! Convective heat transfer coefficient for a smooth surface [W/(m2.K)]
REAL, DIMENSION(SIZE(PCHTCN))                 :: PCHTC   ! Convective heat transfer coefficient [W/(m2.K)]
!
!*       0.2   Declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_ROUGH_DOE',0,ZHOOK_HANDLE)
!
!*       1.    COMPUTE THE CHTC
!              ----------------
!
PCHTC = PCHTCN + PROUGH * (PCHTCS - PCHTCN)
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_ROUGH_DOE',1,ZHOOK_HANDLE)
!
END FUNCTION CHTC_ROUGH_DOE_1D

END MODULE MODE_CONV_DOE
