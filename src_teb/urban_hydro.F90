!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE URBAN_HYDRO(PWS_ROOF_MAX,PWS_ROAD_MAX, PWS_ROOF, PWS_ROAD,  &
                             PRR, PIRRIG_ROAD, PTSTEP, PBLD, PLE_ROOF,     &
                             PLE_ROAD, PRUNOFF_STRLROOF, PRUNOFF_ROAD      )
!   ##########################################################################
!
!!****  *URBAN_HYDRO*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the evolution of prognostic water reservoirs
!     of urbanized areas.
!         
!     
!!**  METHOD
!     ------
!
!
!   The roof reservoir runoff goes directly into the road reservoir.
!
!   Runoff occurs for road reservoir (too much water), as well as drainage
!   (evacuation system, typical time scale: 1 day)
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
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!      V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    23/01/98 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XLVTT
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
REAL, DIMENSION(:), INTENT(IN) :: PWS_ROOF_MAX    ! maximum deepness of roof water reservoir
REAL, DIMENSION(:), INTENT(IN) :: PWS_ROAD_MAX    ! maximum deepness of road water reservoir

REAL, DIMENSION(:), INTENT(INOUT) :: PWS_ROOF     ! roof water reservoir
REAL, DIMENSION(:), INTENT(INOUT) :: PWS_ROAD     ! road water reservoir
REAL, DIMENSION(:), INTENT(IN)    :: PRR          ! rain rate
REAL, DIMENSION(:), INTENT(IN)    :: PIRRIG_ROAD  ! watering rate for roads
REAL,               INTENT(IN)    :: PTSTEP       ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PBLD         ! fraction of buildings
REAL, DIMENSION(:), INTENT(IN)    :: PLE_ROOF     ! latent heat flux over roof
REAL, DIMENSION(:), INTENT(IN)    :: PLE_ROAD     ! latent heat flux over road
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRUNOFF_STRLROOF ! runoff (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT)   :: PRUNOFF_ROAD ! runoff (kg/m2/s)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*      0.2    declarations of local variables
!
!
!-------------------------------------------------------------------------------
!
!*      1.     Roof reservoir evolution
!              ------------------------
!
!
!                                           evolution of the water reservoir
!                                           (if we don't consider the runoff)
!                                           PRR in kg/m2/s therefore PWS in mm
!
IF (LHOOK) CALL DR_HOOK('URBAN_HYDRO',0,ZHOOK_HANDLE)
PWS_ROOF(:) =  PWS_ROOF(:)                                   &
               - PTSTEP * ( PLE_ROOF(:) / XLVTT - PRR(:) )  
!
!                                           Ws_town must be positive
!
PWS_ROOF(:) = MAX(0., PWS_ROOF(:))
!
!                                           if Ws_town > Ws_town_max,
!                                           there is runoff
!
PRUNOFF_STRLROOF(:) = MAX(0., (PWS_ROOF(:) - PWS_ROOF_MAX(:)) / PTSTEP )
!
PWS_ROOF(:) = MIN(PWS_ROOF(:), PWS_ROOF_MAX(:))
!
!-------------------------------------------------------------------------------
!
!*      2.     Road reservoir evolution
!              ------------------------
!
!
!                                           evolution of the water reservoir
!                                           (if we don't consider the runoff)
!                                           PRR in kg/m2/s therefore PWS in mm
!
PWS_ROAD(:) =  PWS_ROAD(:)                                 &
               - PTSTEP * ( PLE_ROAD(:) / XLVTT -  PRR(:) - PIRRIG_ROAD(:) )  
!
!                                           Ws_town must be positive
!
PWS_ROAD(:) = MAX(0., PWS_ROAD(:))
!
!                                           if Ws_town > Ws_town_max,
!                                           there is runoff
!
PRUNOFF_ROAD(:) = MAX(0., (PWS_ROAD(:) - PWS_ROAD_MAX(:)) / PTSTEP )
!
PWS_ROAD(:) = MIN(PWS_ROAD(:), PWS_ROAD_MAX(:))
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('URBAN_HYDRO',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE URBAN_HYDRO
