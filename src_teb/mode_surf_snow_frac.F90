!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ##########################
      MODULE MODE_SURF_SNOW_FRAC
!     ##########################
!
!!****  *MODE_SURF_SNOW_FRAC* -  module for routines to compute snow fraction
!!                               for surface schemes
!!
!!    PURPOSE
!!    -------
!    
!      The purpose of this routine is to store here all routines to compute
!     snow fractions for the TEB scheme. This allows to insure a coherent
!     way in retrieving snow fraction or snow contents.
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
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/03/99
!!     (B.Decharme) 12/03/08  Make sure PPSNV <= PPSNG
!--------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
!-------------------------------------------------------------------------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
CONTAINS
!-------------------------------------------------------------------------------
!
!     ###############################################
      FUNCTION SNOW_FRAC_GROUND(PWSNOW) RESULT(PPSNG)
!     ###############################################
!
USE MODD_SNOW_PAR, ONLY : XWCRN
IMPLICIT NONE
!
REAL, DIMENSION(:), INTENT(IN)  :: PWSNOW ! snow amount over natural areas (kg/m2)
REAL, DIMENSION(SIZE(PWSNOW))   :: PPSNG  ! snow fraction over bare ground
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_GROUND',0,ZHOOK_HANDLE)
PPSNG(:) = PWSNOW(:) / (PWSNOW(:)+XWCRN)       ! fraction of ground covered
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_GROUND',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW_FRAC_GROUND
!
!-------------------------------------------------------------------------------
!
!     ##########################################################
      FUNCTION WSNOW_FROM_SNOW_FRAC_GROUND(PPSNG) RESULT(PWSNOW)
!     ##########################################################
!
USE MODD_SNOW_PAR, ONLY : XWCRN
IMPLICIT NONE
!
REAL, DIMENSION(:), INTENT(IN)  :: PPSNG  ! snow fraction over bare ground
REAL, DIMENSION(SIZE(PPSNG))    :: PWSNOW ! snow amount over natural areas (kg/m2)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:WSNOW_FROM_SNOW_FRAC_GROUND',0,ZHOOK_HANDLE)
PWSNOW(:) = XWCRN * PPSNG(:) / (1. - PPSNG(:))
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:WSNOW_FROM_SNOW_FRAC_GROUND',1,ZHOOK_HANDLE)
!
END FUNCTION WSNOW_FROM_SNOW_FRAC_GROUND
!-------------------------------------------------------------------------------
!
!     #########################################################
      FUNCTION SNOW_FRAC_VEG(PPSNG,PWSNOW,PZ0VEG,PRHOS) RESULT(PPSNV)
!     #########################################################
!
USE MODD_SNOW_PAR, ONLY : XWSNV
IMPLICIT NONE
!
REAL, DIMENSION(:), INTENT(IN)  :: PPSNG  ! snow fraction over bare ground
REAL, DIMENSION(:), INTENT(IN)  :: PWSNOW ! snow amount over natural areas (kg/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PZ0VEG ! vegetation roughness length for momentum
REAL, DIMENSION(:), INTENT(IN)  :: PRHOS  ! snow density (kg/m3)
REAL, DIMENSION(SIZE(PWSNOW))   :: PPSNV  ! snow fraction over vegetation
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_VEG',0,ZHOOK_HANDLE)
PPSNV(:) = PWSNOW(:) / (PWSNOW(:)+PRHOS(:)*XWSNV*PZ0VEG(:))
! Make sure PPSNV <= PPSNG
PPSNV(:) = MIN(PPSNV(:),PPSNG(:))
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_VEG',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW_FRAC_VEG
!
!-------------------------------------------------------------------------------
!     **********************************************************
      FUNCTION SNOW_FRAC_VEG_A(P_PSNG,P_LAI,P_SNOWALB) RESULT(PPSNV)
!     **********************************************************
!
IMPLICIT NONE
!
REAL, DIMENSION(:), INTENT(IN)  :: P_LAI ! leaf area index
REAL, DIMENSION(:), INTENT(IN)  :: P_SNOWALB ! snow albedo
REAL, DIMENSION(:), INTENT(IN)  :: P_PSNG ! snow fraction over bare ground
REAL, DIMENSION(SIZE(P_LAI))    :: PPSNV  ! snow fraction over vegetation
!
!
!
! Definition of local variables
REAL, DIMENSION(SIZE(P_LAI))   :: FLAI  ! snow fraction over vegetation
REAL  RLAIMAX,RLAI,A1,A2
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_VEG_A',0,ZHOOK_HANDLE)
RLAIMAX=7.
RLAI=3.
A1=0.87
A2=0.84
FLAI(:)=1.
WHERE(P_LAI(:)>RLAI)
 FLAI(:)=1.-(P_LAI(:)/RLAIMAX)*(MAX(0.0,(A1-MAX(A2,P_SNOWALB(:))))/(A1-A2))
END WHERE
PPSNV(:)=P_PSNG(:)*FLAI(:)
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_VEG_A',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW_FRAC_VEG_A

!-------------------------------------------------------------------------------
!
!     ############################################################
      FUNCTION SNOW_FRAC_NAT(PWSNOW,PPSNG,PPSNV,PVEG) RESULT(PPSN)
!     ############################################################
!
IMPLICIT NONE
!
REAL, DIMENSION(:), INTENT(IN)  :: PWSNOW ! snow amount over natural areas (kg/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PPSNG  ! snow fraction over bare ground
REAL, DIMENSION(:), INTENT(IN)  :: PPSNV  ! snow fraction over vegetation
REAL, DIMENSION(:), INTENT(IN)  :: PVEG   ! vegetation fraction
REAL, DIMENSION(SIZE(PWSNOW))   :: PPSN   ! snow fraction over natural areas
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_NAT',0,ZHOOK_HANDLE)
PPSN(:) = (1-PVEG(:))*PPSNG(:) + PVEG(:)*PPSNV(:)
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_NAT',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW_FRAC_NAT
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
!     ##############################################################
      SUBROUTINE SNOW_FRAC_ROAD(PWSNOW_ROAD,OSNOW,PDN_ROAD,PDF_ROAD)
!     ##############################################################
!
USE MODD_SNOW_PAR, ONLY : XWCRN
!
REAL, DIMENSION(:), INTENT(IN)  :: PWSNOW_ROAD ! snow amount over roads (kg/m2) 
LOGICAL, DIMENSION(:), INTENT(IN)  :: OSNOW    ! T: snow-fall is occuring
REAL, DIMENSION(:), INTENT(OUT) :: PDN_ROAD    ! snow fraction over roads
REAL, DIMENSION(:), INTENT(OUT) :: PDF_ROAD    ! snow-free fraction over roads
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_ROAD',0,ZHOOK_HANDLE)
PDF_ROAD(:)     = 1.
PDN_ROAD(:)     = 0.
!
! due to the flatness of horizontal surfaces (compared to landscape and
! vegetation), the amount of snow necessary to cover the entire surface XWCRN
! is reduced (equal to 1kg/m2 instead of 10).
!
WHERE (PWSNOW_ROAD(:)>0. .OR. OSNOW)
  PDN_ROAD(:)     = MAX(MIN(PWSNOW_ROAD(:)/(PWSNOW_ROAD(:) + XWCRN*0.1) , 0.7), 0.01)
  PDF_ROAD(:)     = 1.-PDN_ROAD(:)
END WHERE
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_ROAD',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOW_FRAC_ROAD
!
!-------------------------------------------------------------------------------
!
!     ##############################################################
      SUBROUTINE SNOW_FRAC_ROOF(PWSNOW_ROOF,OSNOW,PDN_ROOF,PDF_ROOF)
!     ##############################################################
!
USE MODD_SNOW_PAR, ONLY : XWCRN
!
REAL, DIMENSION(:), INTENT(IN)  :: PWSNOW_ROOF ! snow amount over roofs (kg/m2) 
LOGICAL, DIMENSION(:), INTENT(IN)  :: OSNOW    ! T: snow-fall is occuring
REAL, DIMENSION(:), INTENT(OUT) :: PDN_ROOF    ! snow fraction over roofs
REAL, DIMENSION(:), INTENT(OUT) :: PDF_ROOF    ! snow-free fraction over roofs
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_ROOF',0,ZHOOK_HANDLE)
PDF_ROOF(:)     = 1.
PDN_ROOF(:)     = 0.
!
! due to the flatness of horizontal surfaces (compared to landscape and
! vegetation), the amount of snow necessary to cover the entire surface XWCRN
! is reduced (equal to 1kg/m2 instead of 10).
!
WHERE (PWSNOW_ROOF(:)>0. .OR. OSNOW)
  PDN_ROOF(:)     = MAX(PWSNOW_ROOF(:)/(PWSNOW_ROOF(:) + XWCRN*0.1),0.01)
  PDF_ROOF(:)     = 1.-PDN_ROOF(:)
END WHERE
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_ROOF',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOW_FRAC_ROOF
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! routines bidon pour tora
!
!     ########################################################
      FUNCTION SNOW_FRAC_NAT_1D(PWSNOW)RESULT(BIDON)
!     ########################################################
!
REAL, DIMENSION(:), INTENT(IN)  :: PWSNOW ! snow amount over natural areas (kg/m2)
REAL :: BIDON
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_NAT_1D',0,ZHOOK_HANDLE)
BIDON=PWSNOW(1)
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_NAT_1D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW_FRAC_NAT_1D
!
!-------------------------------------------------------------------------------
!
!     ########################################################
      FUNCTION SNOW_FRAC_NAT_2D(PWSNOW) RESULT(BIDON)
!     ########################################################

REAL :: BIDON
REAL, DIMENSION(:), INTENT(IN)  :: PWSNOW ! snow amount over natural areas (kg/m2)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_NAT_2D',0,ZHOOK_HANDLE)
BIDON=PWSNOW(1)
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_NAT_2D',1,ZHOOK_HANDLE)

END FUNCTION SNOW_FRAC_NAT_2D

!----------------------------------------------------------------------------------
!     ############################################################
      FUNCTION SNOW_FRAC_VEG_1D(PWSNOW) RESULT(BIDON)
!     ############################################################
REAL :: BIDON
REAL, DIMENSION(:), INTENT(IN)  :: PWSNOW ! snow amount over natural areas (kg/m2)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_VEG_1D',0,ZHOOK_HANDLE)
BIDON=PWSNOW(1)
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_VEG_1D',1,ZHOOK_HANDLE)
END FUNCTION SNOW_FRAC_VEG_1D
!
!-------------------------------------------------------------------------------
!
!     ############################################################
      FUNCTION SNOW_FRAC_VEG_2D(PWSNOW) RESULT(BIDON)
!     ############################################################

REAL :: BIDON
REAL, DIMENSION(:), INTENT(IN)  :: PWSNOW ! snow amount over natural areas (kg/m2)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_VEG_2D',0,ZHOOK_HANDLE)
BIDON=PWSNOW(1)
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_VEG_2D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW_FRAC_VEG_2D
!
!-------------------------------------------------------------------------------
!     ##################################################
      FUNCTION SNOW_FRAC_GROUND_1D(PWSNOW) RESULT(BIDON)
!     ##################################################
!
REAL :: BIDON      
REAL, DIMENSION(:), INTENT(IN)  :: PWSNOW ! snow amount over natural areas (kg/m2)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_GROUND_1D',0,ZHOOK_HANDLE)
BIDON=PWSNOW(1)
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_GROUND_1D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW_FRAC_GROUND_1D
!
!-------------------------------------------------------------------------------
!
!     ##################################################
      FUNCTION SNOW_FRAC_GROUND_2D(PWSNOW) RESULT(BIDON)
!     ##################################################
!
REAL :: BIDON
REAL, DIMENSION(:), INTENT(IN)  :: PWSNOW ! snow amount over natural areas (kg/m2)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_GROUND_2D',0,ZHOOK_HANDLE)
BIDON=PWSNOW(1)
IF (LHOOK) CALL DR_HOOK('MODE_SURF_SNOW_FRAC:SNOW_FRAC_GROUND_2D',1,ZHOOK_HANDLE)

!
END FUNCTION SNOW_FRAC_GROUND_2D
!

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
END MODULE MODE_SURF_SNOW_FRAC
