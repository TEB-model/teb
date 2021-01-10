!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #######################
      MODULE MODD_TYPE_EFUTIL
!     #######################
!
!!****  *MODD_TYPE_EFUTIL* - declaration of chemistry emissions types
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to define
!       the types of variables necessary for chemistry emissions. 
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE 
!!
!!    REFERENCE
!!    --------- 
!!      
!!       
!!    AUTHOR
!!    ------
!!      D. Gazen & C. Mari   *LA*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    31/10/00                   
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE 
!
!
INTEGER, PARAMETER :: JPNBCOEFFMAX=50
!
TYPE EMISSVAR_T
  CHARACTER(len=32)             :: CNAME ! Name of the EMISSION 
!
  INTEGER                       :: NTX   ! Current index on array 
!                                          NETIMES (needed in interpolation)
!
  INTEGER                       :: NPX   ! Index to use when t > tmax (periodic)
!                                         
  INTEGER,DIMENSION(:),POINTER  :: NETIMES   ! Emission times (ordered)
!
  INTEGER                       :: NDX       ! Current index on array XEMISDATA
!
  REAL,DIMENSION(:,:), POINTER  :: XEMISDATA ! Emission data array
!                                              with size (X*Y,NWS)  
!
  LOGICAL                       :: LREAD ! TRUE if emission data reading is needed
!
  INTEGER                       :: NWS   ! Read window size for I/O
!
  REAL, DIMENSION(:),POINTER    :: XFWORK ! Time interpolated emission data 
!                                           with size (X) 
END TYPE EMISSVAR_T
!
TYPE PRONOSVAR_T
  INTEGER                        :: NAMINDEX ! pronostic name is CNAMES(NAMINDEX)
!
  INTEGER                        :: NBCOEFF  ! number of coeff in agreg. equation
!                                      
  REAL,   DIMENSION(JPNBCOEFFMAX):: XCOEFF   ! agregation coeff array
!
  INTEGER,DIMENSION(JPNBCOEFFMAX):: NEFINDEX ! array of index on the EMISSVAR_T array
!
  TYPE(PRONOSVAR_T), POINTER     :: NEXT ! To build a PRONOSVAR_T link list
END TYPE PRONOSVAR_T
!
END MODULE MODD_TYPE_EFUTIL
!
