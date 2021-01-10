!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
MODULE MODD_DST_n

!Purpose: 
!Declare variables and constants necessary to do the dust calculations
!Here are only the variables which depend on the grid!

!Author: Alf Grini <alf.grini@cnrm.meteo.fr>

!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE DST_t
  !
  INTEGER, DIMENSION(:), POINTER     :: NVT_DST           !MASK: dust vegetation number to vegetation number
  INTEGER, DIMENSION(:), POINTER   :: NSIZE_PATCH_DST   !Number of points for a patch and a vegetation class
  INTEGER, DIMENSION(:,:), POINTER :: NR_PATCH_DST      !Mask from patch-points to dust-points
  REAL,DIMENSION(:), POINTER         :: Z0_EROD_DST       !Roughness length momentum over erodible dust emitter sfc
  CHARACTER(LEN=6), DIMENSION(:), POINTER  :: CSV_DST     !Name of scalar variables 
  REAL, DIMENSION(:,:),POINTER    :: XSFDST             !Dust variables to be send to output
  REAL, DIMENSION(:,:),POINTER    :: XSFDSTM            !Dust variables to be send to output
  REAL,DIMENSION(:), POINTER   :: XEMISRADIUS_DST         !Number median radius for each source mode
  REAL,DIMENSION(:), POINTER   :: XEMISSIG_DST            !sigma for each source mode
  REAL,DIMENSION(:), POINTER   :: XMSS_FRC_SRC            !Mass fraction of each source mode  
  !
END TYPE DST_t
!
TYPE DST_NP_t
!
TYPE(DST_t), DIMENSION(:), POINTER :: AL=>NULL()
!
END TYPE DST_NP_t
!
CONTAINS
!
SUBROUTINE DST_INIT(YDST)
TYPE(DST_t), INTENT(INOUT) :: YDST
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DST_N:DST_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YDST%NVT_DST)
  NULLIFY(YDST%NSIZE_PATCH_DST)
  NULLIFY(YDST%NR_PATCH_DST)
  NULLIFY(YDST%Z0_EROD_DST)
  NULLIFY(YDST%CSV_DST)
  NULLIFY(YDST%XSFDST)
  NULLIFY(YDST%XSFDSTM)
  NULLIFY(YDST%XEMISRADIUS_DST)
  NULLIFY(YDST%XEMISSIG_DST)
  NULLIFY(YDST%XMSS_FRC_SRC)
IF (LHOOK) CALL DR_HOOK("MODD_DST_N:DST_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DST_INIT
!
SUBROUTINE DST_NP_INIT(YDST_NP,KPATCH)
TYPE(DST_NP_t), INTENT(INOUT) :: YDST_NP
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK("MODD_DST_N:DST_NP_INIT",0,ZHOOK_HANDLE)
!
IF (ASSOCIATED(YDST_NP%AL)) THEN
  DO JP = 1,KPATCH
    CALL DST_INIT(YDST_NP%AL(JP))
  ENDDO
  DEALLOCATE(YDST_NP%AL)
ELSE
  ALLOCATE(YDST_NP%AL(KPATCH))
  DO JP = 1,KPATCH
    CALL DST_INIT(YDST_NP%AL(JP))
  ENDDO
ENDIF
!
IF (LHOOK) CALL DR_HOOK("MODD_DST_N:DST_NP_INIT",1,ZHOOK_HANDLE)
!
END SUBROUTINE DST_NP_INIT
!
END MODULE MODD_DST_n
