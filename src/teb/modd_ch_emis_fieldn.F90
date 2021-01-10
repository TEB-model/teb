!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ###########################
      MODULE MODD_CH_EMIS_FIELD_n
!     ###########################
!
!!****  *MODD_CH_EMIS_FIELD_n* - declaration of chemical emission data arrays
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
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_EFUTIL
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
TYPE CH_EMIS_FIELD_t
!
  REAL               :: XTIME_SIMUL  = 0.
  INTEGER            :: NTIME_MAX
  INTEGER            :: NEMIS_NBR
!                          ! number of chemical pgd fields chosen by user
  CHARACTER(LEN=3) , DIMENSION(:), POINTER :: CEMIS_AREA
!                          ! areas where chemical pgd fields are defined
!                          ! 'ALL' : everywhere
!                          ! 'SEA' : where sea exists
!                          ! 'LAN' : where land exists
!                          ! 'WAT' : where inland water exists
!                          ! 'NAT' : where natural or agricultural areas exist
!                          ! 'TWN' : where town areas exist
!                          ! 'STR' : where streets are present
!                          ! 'BLD' : where buildings are present
!                          !
  CHARACTER(LEN=40), DIMENSION(:), POINTER :: CEMIS_COMMENT ! comment
  CHARACTER(LEN=40), DIMENSION(:), POINTER :: CEMIS_NAME
!                          ! name of the chemical pgd fields (emitted species)
!
  INTEGER,           DIMENSION(:), POINTER :: NEMIS_TIME   ! emission time
!
  REAL,              DIMENSION(:,:), POINTER:: XEMIS_FIELDS ! emission pgd fields values
!
  INTEGER                                          :: NEMISPEC_NBR ! Number of chemical species
!
  TYPE(EMISSVAR_T),  DIMENSION(:), POINTER :: TSEMISS      ! Offline emission struct array
!
  TYPE(PRONOSVAR_T),               POINTER     :: TSPRONOSLIST ! Head pointer on pronostic
!                                                              variables list
!-------------------------------------------------------------------------------
!
END TYPE CH_EMIS_FIELD_t



CONTAINS

!




SUBROUTINE CH_EMIS_FIELD_INIT(YCH_EMIS_FIELD)
TYPE(CH_EMIS_FIELD_t), INTENT(INOUT) :: YCH_EMIS_FIELD
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_EMIS_FIELD_N:CH_EMIS_FIELD_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YCH_EMIS_FIELD%CEMIS_AREA)
  NULLIFY(YCH_EMIS_FIELD%CEMIS_COMMENT)
  NULLIFY(YCH_EMIS_FIELD%CEMIS_NAME)
  NULLIFY(YCH_EMIS_FIELD%NEMIS_TIME)
  NULLIFY(YCH_EMIS_FIELD%XEMIS_FIELDS)
  NULLIFY(YCH_EMIS_FIELD%TSEMISS)
YCH_EMIS_FIELD%XTIME_SIMUL=0.
YCH_EMIS_FIELD%NEMIS_NBR=0
YCH_EMIS_FIELD%NTIME_MAX=-1
YCH_EMIS_FIELD%NEMISPEC_NBR=0
IF (LHOOK) CALL DR_HOOK("MODD_CH_EMIS_FIELD_N:CH_EMIS_FIELD_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE CH_EMIS_FIELD_INIT


END MODULE MODD_CH_EMIS_FIELD_n

