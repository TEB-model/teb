!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE BLD_E_BUDGET( OTI_EVOL, PTSTEP, PBLD, PWALL_O_HOR,      &
                             PRHOA, PT_ROOF, PT_WALL, PTI_BLD, PTS_FLOOR )  
!   ##########################################################################
!
!!****  *BLD_E_BUDGET*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the evoultion of the temperature of inside building air
        
!     
!!**  METHOD
!     ------
!
!     The resistance term between the surfaces and the room is given
!     by a standard value, which mimics both the convection
!     and the radiative interactions in the room.
!     This explains the very low resistance. It is used to compute
!     the evolution of the surfaces only.
!     This resistance value is 0.123 Km/W  (typical for inside surfaces).
!     (ENVIRONMENTAL SCIENCE IN BUILDING, 3rd Edition, Randall McMullan,
!      THE MACMILLAN PRESS Limited).
!
!
!
!     On the contrary, the evolution of the air temperature is mainly
!     governed by the convection (considering the low radiative absorption
!     of the air itself).
!     In order to have a simple formulation, a diurnal cycle is assumed,
!     with a force restore formulation.
!
!     The floor temperature is fixed
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
!!      Original    24/08/00 
!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XTT, XCPD, XDAY
USE MODD_SURF_PAR,ONLY : XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
LOGICAL,              INTENT(IN)   :: OTI_EVOL      ! true --> internal temp. of
!                                                   !      of buildings evolves
!                                                   ! false--> it is fixed
REAL,                 INTENT(IN)   :: PTSTEP        ! time step
REAL, DIMENSION(:),   INTENT(IN)   :: PBLD          ! building fraction
REAL, DIMENSION(:),   INTENT(IN)   :: PWALL_O_HOR   ! wall surf. / hor. surf.
REAL, DIMENSION(:),   INTENT(IN)   :: PRHOA         ! air density
                                                    ! at the lowest level
REAL, DIMENSION(:,:), INTENT(IN)   :: PT_ROOF       ! roof layers temperatures
REAL, DIMENSION(:,:), INTENT(IN)   :: PT_WALL       ! wall layers temperatures
REAL, DIMENSION(:),   INTENT(INOUT):: PTI_BLD       ! building air temperature
                                                    ! computed with its equation evolution
REAL, DIMENSION(:),   INTENT(IN)  :: PTS_FLOOR     ! floor surface temperature
!
!*      0.2    declarations of local variables
!
!
REAL                           :: ZTAU         ! temporal filter period
!
INTEGER                        :: IROOF        ! number of roof layers
INTEGER                        :: IWALL        ! number of wall layers
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('BLD_E_BUDGET',0,ZHOOK_HANDLE)
!
!*      1.   initializations
!            ---------------
!
IROOF = SIZE(PT_ROOF,2)
IWALL = SIZE(PT_WALL,2)
!
!!! 27/01/2012 passé dans TEB
!! PTS_FLOOR(:)= 19. + XTT
!!! 27/01/2012 passé dans TEB
!
!*      2.   no evolution of interior temperature if OTI_EVOL=.FALSE.
!            --------------------------------------------------------
!
IF (.NOT. OTI_EVOL .AND. LHOOK) CALL DR_HOOK('BLD_E_BUDGET',1,ZHOOK_HANDLE)
IF (.NOT. OTI_EVOL) RETURN
!
!*      3.   evolution of the internal temperature
!            -------------------------------------
!
ZTAU = XDAY
!
WHERE (PBLD(:) .GT. 0.)
  PTI_BLD(:) = PTI_BLD(:) * (ZTAU-PTSTEP)/ZTAU       &
            + ( PT_ROOF(:,IROOF) * PBLD       (:)    &
              + PT_WALL(:,IWALL) * PWALL_O_HOR(:)    &
              + PTS_FLOOR(:)      * PBLD       (:) )  &
             / (  2. * PBLD(:)  +  PWALL_O_HOR(:) ) * PTSTEP / ZTAU
ELSEWHERE
  PTI_BLD   (:) = PTS_FLOOR(:)
ENDWHERE
!
!
!*      5.   internal temperature set to a minimum value (heating)
!            -----------------------------------------------------
! 
PTI_BLD(:) = MAX( PTI_BLD(:) , PTS_FLOOR (:) )
!
IF (LHOOK) CALL DR_HOOK('BLD_E_BUDGET',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE BLD_E_BUDGET
