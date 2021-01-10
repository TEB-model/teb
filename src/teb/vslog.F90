!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
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
