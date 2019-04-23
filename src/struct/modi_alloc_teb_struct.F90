!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
MODULE MODI_ALLOC_TEB_STRUCT
!
INTERFACE
!
!   ##########################################################################
    SUBROUTINE ALLOC_TEB_STRUCT(KROOF_LAYER,KROAD_LAYER,KWALL_LAYER,KFLOOR_LAYER, &
                                TOP,T,BOP,B,TPN,TIR,DMT                           )
!   ##########################################################################
!
!
!*      0.1    Declarations of arguments
!
!
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_TEB_PANEL_n, ONLY : TEB_PANEL_t
USE MODD_TEB_IRRIG_n, ONLY : TEB_IRRIG_t
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_t
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KROOF_LAYER   ! number of roof layers
INTEGER, INTENT(IN) :: KROAD_LAYER   ! number of road layers
INTEGER, INTENT(IN) :: KWALL_LAYER   ! number of wall layers
INTEGER, INTENT(IN) :: KFLOOR_LAYER  ! number of floor layers
TYPE(TEB_OPTIONS_t),   INTENT(INOUT) :: TOP
TYPE(TEB_t),           INTENT(INOUT) :: T
TYPE(BEM_OPTIONS_t),   INTENT(INOUT) :: BOP
TYPE(BEM_t),           INTENT(INOUT) :: B
TYPE(TEB_PANEL_t),     INTENT(INOUT) :: TPN
TYPE(TEB_IRRIG_t),     INTENT(INOUT) :: TIR
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: DMT
!
END SUBROUTINE ALLOC_TEB_STRUCT
!
END INTERFACE
!
END MODULE MODI_ALLOC_TEB_STRUCT
