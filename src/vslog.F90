!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL version 2.1
!SFX_LIC version 1. See LICENSE, Licence_CeCILL_V2.1-en.txt and Licence_CeCILL_V2.1-fr.txt  
!SFX_LIC for details. version 1.
      SUBROUTINE VSLOG(PA,PLOG,N)
!
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE PARKIND1  ,ONLY : JPRB
!
!
!

!   COMPUTES THE LOGARITHM

      IMPLICIT NONE

      INTEGER :: N
      REAL :: PA(N), PLOG(N)

      INTEGER :: J
      REAL(KIND=JPRB) :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('VSLOG',0,ZHOOK_HANDLE)
      DO J=1,N
        PLOG(J) = LOG(PA(J))
      END DO
      IF (LHOOK) CALL DR_HOOK('VSLOG',1,ZHOOK_HANDLE)

      END
