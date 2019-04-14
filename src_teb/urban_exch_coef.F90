!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE URBAN_EXCH_COEF(HZ0H, PZ0_O_Z0H, PTG, PQS, PEXNS, PEXNA, PTA, PQA,   &
                             PZREF, PUREF, PVMOD, PZ0,                            &
                             PRI, PCD, PCDN, PAC, PRA, PCH                        )  
!          #######################################################################
!
!!****  *URBAN_DRAG*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the surface drag coefficients over roof, road or town
!     according several scientific formulations
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
!!      Original    01/2009    from urban_drag.f90 (modified by S.  Leroyer at CMC)
!         01/2009 (S. Leroyer) option (HZ0H) for z0h applied on roof, road and town
!!      B. Decharme    06/2009 limitation of Ri
!!      B. Decharme    09/2009 limitation of Ri in surface_ri.F90
!
!-------------------------------------------------------------------------------

USE MODI_SURFACE_RI
USE MODI_SURFACE_CD
USE MODI_SURFACE_AERO_COND
USE MODI_WIND_THRESHOLD
!
USE MODD_CSTS, ONLY : XKARMAN
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_FLXSURF3BX
!
IMPLICIT NONE
!
 CHARACTER(LEN=6)                  :: HZ0H     ! TEB option for z0h roof & road
REAL,               INTENT(IN)    :: PZ0_O_Z0H! z0/z0h ratio used in Mascart (1995)
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
REAL, DIMENSION(:), INTENT(IN)    :: PZ0      ! roughness length for momentum
REAL, DIMENSION(:), INTENT(IN)    :: PZREF    ! reference height of the first
                                              ! atmospheric level
REAL, DIMENSION(:), INTENT(IN)    :: PUREF    ! reference height of the wind
!                                             ! NOTE this is different from ZZREF
!                                             ! ONLY in stand-alone/forced mode,
!                                             ! NOT when coupled to a model (MesoNH)
REAL, DIMENSION(:), INTENT(OUT)   :: PRI      ! Richardson number
!
REAL, DIMENSION(:), INTENT(OUT)   :: PCD      ! drag coefficient for momentum
REAL, DIMENSION(:), INTENT(OUT)   :: PCDN     ! neutral drag coefficient for momentum
REAL, DIMENSION(:), INTENT(OUT)   :: PAC      ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PRA      ! aerodynamical resistance
REAL, DIMENSION(:), INTENT(OUT)   :: PCH      ! drag coefficient for heat
!
!* local variables
!
REAL, DIMENSION(SIZE(PZ0))         :: ZZ0H     ! roughness length for heat
!* BRUT82 & KAND07 cases
REAL,DIMENSION(SIZE(PTA)) :: cmu, ctu, rib,ftemp,fvap,ilmo  ! temporary var for
REAL,DIMENSION(SIZE(PTA)) :: ue,fcor, hBL,lzz0,lzz0t,fm, fh ! flxsurf3
REAL,DIMENSION(SIZE(PTA)) :: z0h_roof,z0h_town,z0h_road     ! local thermal roughness
REAL,DIMENSION(SIZE(PTA)) :: zustar, zta, ztg
REAL,DIMENSION(SIZE(PTA)) :: ZVMOD                          ! wind
INTEGER N
!
!* MASC95 case
REAL,DIMENSION(SIZE(PTA)) :: ZDIRCOSZW     ! orography slope cosine
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*      1.     Mascart 1995 exchange coefficients
!              ----------------------------------
!
IF (LHOOK) CALL DR_HOOK('URBAN_EXCH_COEF',0,ZHOOK_HANDLE)
IF (HZ0H=='MASC95') THEN
! 
  ZZ0H = PZ0 / PZ0_O_Z0H ! fixed ratio for MASC95
  ZDIRCOSZW=1.           ! no orography slope effect taken into account in TEB
!
  CALL SURFACE_RI(PTG, PQS, PEXNS, PEXNA, PTA, PQA,   &
                    PZREF, PUREF, ZDIRCOSZW, PVMOD, PRI )  
!
  CALL SURFACE_CD(PRI, PZREF, PUREF, PZ0, ZZ0H, PCD, PCDN)
!
  CALL SURFACE_AERO_COND(PRI, PZREF, PUREF, PVMOD, PZ0, ZZ0H, PAC, PRA, PCH)
!
!
!*      2.     Brutsaert 1982  or Kanda 2007 exchange coefficients
!              ---------------------------------------------------
!
ELSEIF(HZ0H=='BRUT82' .OR. HZ0H=='KAND07')THEN
  ! initialisations
  fcor(:)=1.0372462E-04
!RJ: can be removed
  N=SIZE(PTA)
  !
  ! Set a minimum threshold to the wind
  ZVMOD(:) = WIND_THRESHOLD(PVMOD(:),PUREF(:))
  !
  ! First guess of u*
  ZUSTAR(:) = 0.4 * ZVMOD(:) / LOG( PUREF/PZ0(:) )
  !
  IF (HZ0H=='KAND07') THEN ! Kanda 2007
    ZZ0H(:)= PZ0(:) * 7.4 * EXP( - 1.29 *( PZ0(:)*zustar(:)/1.461e-5)**0.25)
  ELSEIF (HZ0H=='BRUT82') THEN ! Brutsaert 1982
    ZZ0H(:)= PZ0(:) * 7.4 * EXP( - 2.46 *( PZ0(:)*zustar(:)/1.461e-5)**0.25)
  ENDIF

  ZTA = PTA/PEXNA
  ZTG = PTG/PEXNS
  CALL FLXSURF3BX( cmu, ctu, PRI,ftemp,fvap,ilmo,    &
                   ue, fcor, ZTA, PQA,           &
                   PUREF, PZREF, ZVMOD, ZTG, PQS,&
                   hBL, PZ0,ZZ0H,                      &
                   lzz0, lzz0t, fm, fh,N               )  
!
  PCD(:) = (cmu(:)/ue(:))**2  
  PCDN(:) = (XKARMAN/LOG(PUREF(:)/PZ0(:)))**2
  PAC(:) = (cmu(:)*ctu(:)/ue(:)**2) * ZVMOD(:)
  PRA(:) = 1. / PAC(:)
  PCH(:) = 1. / (PRA(:) * ZVMOD(:))
!
END IF
IF (LHOOK) CALL DR_HOOK('URBAN_EXCH_COEF',1,ZHOOK_HANDLE)

END SUBROUTINE URBAN_EXCH_COEF
