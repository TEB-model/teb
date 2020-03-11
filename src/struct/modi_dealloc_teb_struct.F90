!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL licence
!SFX_LIC version 2.1. See Licence_CeCILL_V2.1-en.txt and Licence_CeCILL_V2.1-fr.txt  
!SFX_LIC for details.
MODULE MODI_DEALLOC_TEB_STRUCT
!
INTERFACE
!
!   ##########################################################################
    SUBROUTINE DEALLOC_TEB_STRUCT(TOP,T,BOP,B,TPN,TIR,DMT)
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
TYPE(TEB_OPTIONS_t),   INTENT(INOUT) :: TOP
TYPE(TEB_t),           INTENT(INOUT) :: T
TYPE(BEM_OPTIONS_t),   INTENT(INOUT) :: BOP
TYPE(BEM_t),           INTENT(INOUT) :: B
TYPE(TEB_PANEL_t),     INTENT(INOUT) :: TPN
TYPE(TEB_IRRIG_t),     INTENT(INOUT) :: TIR
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: DMT
!
END SUBROUTINE DEALLOC_TEB_STRUCT
!
END INTERFACE
!
END MODULE MODI_DEALLOC_TEB_STRUCT
