!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL licence version 2.1.
! See the following links for details:
! https://cecill.info/licences/Licence_CeCILL_V2.1-en.txt
! https://cecill.info/licences/Licence_CeCILL_V2.1-fr.txt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ######################
MODULE MODD_FORC_ATM
!     ######################
!
!!****  *MODD_FORC_ATM - declaration of atmospheric forcing variables
!!
!!    PURPOSE
!!    -------
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
!!	F. Habets   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       20/09/02
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!------------------------------------------------------------------------------
!
CHARACTER(LEN=6), DIMENSION(:),  ALLOCATABLE :: CSV       ! name of all scalar variables
REAL,             DIMENSION(:,:),ALLOCATABLE :: XDIR_ALB  ! direct albedo for each band
REAL,             DIMENSION(:,:),ALLOCATABLE :: XSCA_ALB  ! diffuse albedo for each band
REAL,             DIMENSION(:),  ALLOCATABLE :: XEMIS     ! emissivity
REAL,             DIMENSION(:),  ALLOCATABLE :: XTSRAD    ! radiative temperature

REAL, DIMENSION(:), ALLOCATABLE  :: XTSUN     ! solar time                    (s from midnight)
REAL, DIMENSION(:), ALLOCATABLE  :: XZREF     ! height of T,q forcing                 (m)
REAL, DIMENSION(:), ALLOCATABLE  :: XUREF     ! height of wind forcing                (m)
!
REAL, DIMENSION(:), ALLOCATABLE  :: XTA       ! air temperature forcing               (K)
REAL, DIMENSION(:), ALLOCATABLE  :: XQA       ! air specific humidity forcing         (kg/m3)
REAL, DIMENSION(:), ALLOCATABLE  :: XRHOA     ! air density forcing                   (kg/m3)
REAL, DIMENSION(:,:),ALLOCATABLE :: XSV       ! scalar variables
REAL, DIMENSION(:), ALLOCATABLE  :: XU        ! zonal wind                            (m/s)
REAL, DIMENSION(:), ALLOCATABLE  :: XV        ! meridian wind                         (m/s)
REAL, DIMENSION(:,:),ALLOCATABLE :: XDIR_SW   ! direct  solar radiation (on horizontal surf.)
!                                            !                                       (W/m2)
REAL, DIMENSION(:,:),ALLOCATABLE :: XSCA_SW   ! diffuse solar radiation (on horizontal surf.)
!                                            !                                       (W/m2)
REAL, DIMENSION(:), ALLOCATABLE  :: XSW_BANDS ! mean wavelength of each shortwave band (m)
REAL, DIMENSION(:), ALLOCATABLE  :: XZENITH   ! zenithal angle at t  (radian from the vertical)
REAL, DIMENSION(:), ALLOCATABLE  :: XZENITH2  ! zenithal angle at t+1(radian from the vertical)
REAL, DIMENSION(:), ALLOCATABLE  :: XAZIM     ! azimuthal angle      (radian from North, clockwise)
REAL, DIMENSION(:), ALLOCATABLE  :: XLW       ! longwave radiation (on horizontal surf.)
!                                            !                                       (W/m2)
REAL, DIMENSION(:), ALLOCATABLE  :: XPS       ! pressure at atmospheric model surface (Pa)
REAL, DIMENSION(:), ALLOCATABLE  :: XPA       ! pressure at forcing level             (Pa)
REAL, DIMENSION(:), ALLOCATABLE  :: XZS       ! atmospheric model orography           (m)
REAL, DIMENSION(:), ALLOCATABLE  :: XCO2      ! CO2 concentration in the air          (kg/kg)
REAL, DIMENSION(:), ALLOCATABLE  :: XSNOW     ! snow precipitation                    (kg/m2/s)
REAL, DIMENSION(:), ALLOCATABLE  :: XRAIN     ! liquid precipitation                  (kg/m2/s)
!
!
REAL, DIMENSION(:), ALLOCATABLE :: XSFTH     ! flux of heat                          (W/m2)
REAL, DIMENSION(:), ALLOCATABLE :: XSFTQ     ! flux of water vapor                   (kg/m2/s)
REAL, DIMENSION(:), ALLOCATABLE :: XSFU      ! zonal momentum flux                   (pa)
REAL, DIMENSION(:), ALLOCATABLE :: XSFV      ! meridian momentum flux                (pa)
REAL, DIMENSION(:), ALLOCATABLE :: XSFCO2    ! flux of CO2                           (kg/m2/s)
REAL, DIMENSION(:,:),ALLOCATABLE:: XSFTS     ! flux of scalar var.                   (kg/m2/s)
!
REAL, DIMENSION(:), ALLOCATABLE :: XPEW_A_COEF ! implicit coefficients
REAL, DIMENSION(:), ALLOCATABLE :: XPEW_B_COEF ! needed if HCOUPLING='I'
REAL, DIMENSION(:), ALLOCATABLE :: XPET_A_COEF
REAL, DIMENSION(:), ALLOCATABLE :: XPEQ_A_COEF
REAL, DIMENSION(:), ALLOCATABLE :: XPET_B_COEF
REAL, DIMENSION(:), ALLOCATABLE :: XPEQ_B_COEF

!------------------------------------------------------------------------------
!
END MODULE MODD_FORC_ATM

