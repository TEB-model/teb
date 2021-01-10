!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ###########################
      MODULE MODD_EMIS_GR_FIELD_n
!     ###########################
!
!!****  *MODD_EMIS_GR_FIELD_n* - declaration of chemical emission data arrays
!!                               for model n
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     chemical emission data arrays for model n.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_EMIS_GR_FIELD)
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

TYPE EMIS_GR_FIELD_t
!
!
  INTEGER                                     :: NEMIS_GR_NBR
!                          ! number of chemical pgd fields chosen by user
  CHARACTER(LEN=3) , DIMENSION(:), POINTER:: CEMIS_GR_AREA
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
  CHARACTER(LEN=40), DIMENSION(:), POINTER:: CEMIS_GR_NAME
!                          ! name of the chemical pgd fields (emitted species)
!
  INTEGER,       DIMENSION(:),     POINTER:: NEMIS_GR_TIME   ! emission time
!
  REAL,          DIMENSION(:,:,:), POINTER:: XEMIS_GR_FIELDS ! emission pgd fields values
!
  INTEGER                                          :: NEMISPEC_NBR ! Number of chemical species
!
  TYPE(EMISSVAR_T),  DIMENSION(:), POINTER :: TSEMISS      ! Offline emission struct array
!
  TYPE(PRONOSVAR_T),               POINTER     :: TSPRONOSLIST ! Head pointer on pronostic
!                                                              variables list
!-------------------------------------------------------------------------------
!
END TYPE EMIS_GR_FIELD_t



CONTAINS

!




SUBROUTINE EMIS_GR_FIELD_INIT(YEMIS_GR_FIELD)
TYPE(EMIS_GR_FIELD_t), INTENT(INOUT) :: YEMIS_GR_FIELD
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_EMIS_GR_FIELD_N:EMIS_GR_FIELD_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YEMIS_GR_FIELD%CEMIS_GR_AREA)
  NULLIFY(YEMIS_GR_FIELD%CEMIS_GR_NAME)
  NULLIFY(YEMIS_GR_FIELD%NEMIS_GR_TIME)
  NULLIFY(YEMIS_GR_FIELD%XEMIS_GR_FIELDS)
YEMIS_GR_FIELD%NEMIS_GR_NBR=0
YEMIS_GR_FIELD%NEMISPEC_NBR=0
IF (LHOOK) CALL DR_HOOK("MODD_EMIS_GR_FIELD_N:EMIS_GR_FIELD_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE EMIS_GR_FIELD_INIT


END MODULE MODD_EMIS_GR_FIELD_n
