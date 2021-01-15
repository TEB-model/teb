!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL version 2.1
!SFX_LIC version 1. See LICENSE, Licence_CeCILL_V2.1-en.txt and Licence_CeCILL_V2.1-fr.txt  
!SFX_LIC for details. version 1.
MODULE MODD_SLT_n

!Purpose: 
!Declare variables and constants necessary to do the sea salt calculations
!Here are only the variables which depend on the grid!
!
!Author: Alf Grini / Pierre Tulet
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
TYPE SLT_t
  REAL,DIMENSION(:), POINTER     :: XEMISRADIUS_SLT             ! Number median radius for each source mode
  REAL,DIMENSION(:), POINTER     :: XEMISSIG_SLT                ! sigma for each source mode
END TYPE SLT_t
!
CONTAINS
!
SUBROUTINE SLT_INIT(YSLT)
TYPE(SLT_t), INTENT(INOUT) :: YSLT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_SLT_N:SLT_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YSLT%XEMISRADIUS_SLT)
  NULLIFY(YSLT%XEMISSIG_SLT)
IF (LHOOK) CALL DR_HOOK("MODD_SLT_N:SLT_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE SLT_INIT
!
END MODULE MODD_SLT_n
