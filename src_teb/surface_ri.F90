!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     #########
    SUBROUTINE SURFACE_RI(PTG, PQS, PEXNS, PEXNA, PTA, PQA,   &
                               PZREF, PUREF, PDIRCOSZW, PVMOD, PRI )  
!   ######################################################################
!
!!****  *SURFACE_RI*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the richardson number near the ground
!         
!     
!!**  METHOD
!!    ------
!
!
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    MODD_CST
!!    MODD_GROUND_PAR
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!	V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    22/09/98 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,     ONLY : XRV, XRD, XG
USE MODD_SURF_ATM, ONLY : XRIMAX
USE MODI_WIND_THRESHOLD
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
REAL, DIMENSION(:), INTENT(IN)    :: PTG      ! surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PQS      ! surface specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS    ! surface exner function
REAL, DIMENSION(:), INTENT(IN)    :: PTA      ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA      ! specific humidity
                                              ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PEXNA    ! exner function
                                              ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD    ! module of the horizontal wind
!
REAL, DIMENSION(:), INTENT(IN)    :: PZREF    ! reference height of the first
                                              ! atmospheric level
REAL, DIMENSION(:), INTENT(IN)    :: PUREF    ! reference height of the wind
!                                             ! NOTE this is different from ZZREF
!                                             ! ONLY in stand-alone/forced mode,
!                                             ! NOT when coupled to a model (MesoNH)
REAL, DIMENSION(:), INTENT(IN)    :: PDIRCOSZW! Cosine of the angle between
!                                             ! the normal to the surface and
!                                             ! the vertical
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRI      ! Richardson number
!
!*      0.2    declarations of local variables
!
!
REAL, DIMENSION(SIZE(PTG))   :: ZTHVA, ZTHVS
REAL, DIMENSION(SIZE(PVMOD)) :: ZVMOD
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!       1.     Richardson number
!              -----------------
!                
!                                                 virtual potential        
!                                                 temperature at the 
!                                                 first atmospheric level and
!                                                 at the surface
!
IF (LHOOK) CALL DR_HOOK('SURFACE_RI',0,ZHOOK_HANDLE)
!
ZTHVA(:)=PTA(:)/PEXNA(:)*( 1.+(XRV/XRD-1.)*PQA(:) )   
ZTHVS(:)=PTG(:)/PEXNS(:)*( 1.+(XRV/XRD-1.)*PQS(:) )
!                                                 
ZVMOD(:) = WIND_THRESHOLD(PVMOD(:),PUREF(:))
!
                                                ! Richardson's number
PRI(:) = XG * PDIRCOSZW(:) * PUREF(:) * PUREF(:)              &
          * (ZTHVA(:)-ZTHVS(:)) / (0.5 * (ZTHVA(:)+ZTHVS(:)) )  &
          / (ZVMOD(:)*ZVMOD(:)) /PZREF(:)  
!
PRI(:) = MIN(PRI(:),XRIMAX)
!
IF (LHOOK) CALL DR_HOOK('SURFACE_RI',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE SURFACE_RI
