!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL version 2.1
!SFX_LIC version 1. See LICENSE, Licence_CeCILL_V2.1-en.txt and Licence_CeCILL_V2.1-fr.txt  
!SFX_LIC for details. version 1.
!     ###########################
      MODULE MODD_CH_SNAP_n
!     ###########################
!
!!****  *MODD_CH_SNAP_n* - declaration of chemical emission data arrays
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     chemical emission data arrays.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!      D. Gazen   *L.A.*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/03/2001                      
!!      01/12/03    (D.Gazen) change emissions handling for surf. externalization
!!      M.Leriche 04/2014  change length of CHARACTER for emission 6->12
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_EFUTIL
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
TYPE CH_EMIS_SNAP_t
!
  INTEGER            :: NEMIS_NBR
!                          ! number of chemical pgd fields chosen by user
  CHARACTER(LEN=3)                         :: CCONVERSION
!                          ! Unit conversion code
  CHARACTER(LEN=5)                         :: CSNAP_TIME_REF
!                          ! Reference time for Snap temporal profiles
!                          !  'UTC  ' : UTC   time
!                          !  'SOLAR' : SOLAR time
!                          !  'LEGAL' : LEGAL time
!                          !

  CHARACTER(LEN=12), DIMENSION(:), POINTER :: CEMIS_NAME
!                          ! name of the chemical fields (emitted species)
  CHARACTER(LEN=40), DIMENSION(:), POINTER :: CEMIS_COMMENT
!                          ! comment on the chemical fields (emitted species)
!
  REAL,     DIMENSION(:,:,:), POINTER:: XEMIS_FIELDS_SNAP ! Emission factor for
!                                                         ! each chemical specie and
!                                                         ! each snap
  REAL,     DIMENSION(:,:),   POINTER:: XEMIS_FIELDS      ! Emission for each specie
!                                                         ! (at a given time taking into 
!                                                         ! account all snaps)
  REAL,     DIMENSION(:),     POINTER:: XDELTA_LEGAL_TIME ! Difference (in hours)) between
!                                                         ! Legal time and UTC time
  INTEGER            :: NEMIS_SNAP                        ! number of snaps
  INTEGER            :: NSNAP_M                           ! number of months
  INTEGER            :: NSNAP_D                           ! number of days
  INTEGER            :: NSNAP_H                           ! number of hours
  REAL,              DIMENSION(:,:,:), POINTER:: XSNAP_MONTHLY
  REAL,              DIMENSION(:,:,:), POINTER:: XSNAP_DAILY
  REAL,              DIMENSION(:,:,:), POINTER:: XSNAP_HOURLY
  REAL,              DIMENSION(:),     POINTER:: XCONVERSION ! conversion factor
!
  TYPE(PRONOSVAR_T),               POINTER     :: TSPRONOSLIST ! Head pointer on pronostic
!                                                              variables list
!-------------------------------------------------------------------------------
!
END TYPE CH_EMIS_SNAP_t



CONTAINS

!




SUBROUTINE CH_EMIS_SNAP_INIT(YCH_EMIS_SNAP)
TYPE(CH_EMIS_SNAP_t), INTENT(INOUT) :: YCH_EMIS_SNAP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_SNAP_n:CH_EMIS_FIELD_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YCH_EMIS_SNAP%CEMIS_COMMENT)
  NULLIFY(YCH_EMIS_SNAP%CEMIS_NAME)
  NULLIFY(YCH_EMIS_SNAP%XDELTA_LEGAL_TIME)
  NULLIFY(YCH_EMIS_SNAP%XEMIS_FIELDS)
  NULLIFY(YCH_EMIS_SNAP%XEMIS_FIELDS_SNAP)
  NULLIFY(YCH_EMIS_SNAP%XSNAP_DAILY)
  NULLIFY(YCH_EMIS_SNAP%XSNAP_HOURLY)
  NULLIFY(YCH_EMIS_SNAP%XSNAP_MONTHLY)
  NULLIFY(YCH_EMIS_SNAP%XCONVERSION)
YCH_EMIS_SNAP%CCONVERSION=' '
YCH_EMIS_SNAP%CSNAP_TIME_REF=' '
YCH_EMIS_SNAP%NEMIS_NBR=0
YCH_EMIS_SNAP%NEMIS_SNAP=0
YCH_EMIS_SNAP%NSNAP_M=0
YCH_EMIS_SNAP%NSNAP_D=0
YCH_EMIS_SNAP%NSNAP_H=0
IF (LHOOK) CALL DR_HOOK("MODD_CH_SNAP_n:CH_EMIS_FIELD_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE CH_EMIS_SNAP_INIT


END MODULE MODD_CH_SNAP_n

