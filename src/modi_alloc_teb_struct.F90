!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL licence
!SFX_LIC version 2.1. See Licence_CeCILL_V2.1-en.txt and Licence_CeCILL_V2.1-fr.txt  
!SFX_LIC for details.
MODULE MODI_ALLOC_TEB_STRUCT
!
INTERFACE
!
!   ##########################################################################
    SUBROUTINE ALLOC_TEB_STRUCT(KROOF_LAYER,KROAD_LAYER,KWALL_LAYER,KFLOOR_LAYER,&
                               DTCO, G, TOP, T, BOP, B, TPN, TIR, DMT, GDM, GRM )
!   ##########################################################################
!
!!****  *ALLOC_TEB_STRUCT*
!!
!!    PURPOSE
!!    -------
!
!       To allocate all the structures used in Teb_garden and subsequent routines.
!
!!**  METHOD
!     ------
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
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
!!    Original    03/2018
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t, DATA_COVER_INIT
USE MODD_SFX_GRID_n, ONLY : GRID_t, GRID_INIT
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t, TEB_OPTIONS_INIT
USE MODD_TEB_n, ONLY : TEB_t, TEB_INIT
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t, BEM_OPTIONS_INIT
USE MODD_BEM_n, ONLY : BEM_t, BEM_INIT
USE MODD_TEB_PANEL_n, ONLY : TEB_PANEL_t, TEB_PANEL_INIT
USE MODD_TEB_IRRIG_n, ONLY : TEB_IRRIG_t, TEB_IRRIG_INIT
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_t, DIAG_MISC_TEB_INIT
USE MODD_SURFEX_n, ONLY : TEB_GARDEN_MODEL_t
USE MODD_SURFEX_n, ONLY : TEB_GREENROOF_MODEL_t
!
!
IMPLICIT NONE
!
!*      0.1    Declarations of arguments
!
!
!
!
INTEGER, INTENT(IN) :: KROOF_LAYER   ! number of roof layers
INTEGER, INTENT(IN) :: KROAD_LAYER   ! number of road layers
INTEGER, INTENT(IN) :: KWALL_LAYER   ! number of wall layers
INTEGER, INTENT(IN) :: KFLOOR_LAYER  ! number of floor layers
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(GRID_t), INTENT(INOUT) :: G
TYPE(TEB_OPTIONS_t),   INTENT(INOUT) :: TOP
TYPE(TEB_t),           INTENT(INOUT) :: T
TYPE(BEM_OPTIONS_t),   INTENT(INOUT) :: BOP
TYPE(BEM_t),           INTENT(INOUT) :: B
TYPE(TEB_PANEL_t),     INTENT(INOUT) :: TPN
TYPE(TEB_IRRIG_t),     INTENT(INOUT) :: TIR
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: DMT
TYPE(TEB_GARDEN_MODEL_t), INTENT(INOUT) :: GDM
TYPE(TEB_GREENROOF_MODEL_t), INTENT(INOUT) :: GRM
!
END SUBROUTINE ALLOC_TEB_STRUCT
!
END INTERFACE
!
END MODULE MODI_ALLOC_TEB_STRUCT
