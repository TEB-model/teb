!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-------------------------------------- LICENCE BEGIN ------------------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------------------
!copyright (C) 2001  MSC-RPN COMM  %%%RPNPHY%%%
!!!S/P  FLXSURF3
!
      SUBROUTINE FLXSURF3BX(CMU, CTU, RIB, FTEMP, FVAP, ILMO,           &
     &                   UE, FCOR, TA , QA , ZU, ZT, VA,                &
     &                   TG , QG , H , Z0 , Z0T,                        &
     &                   LZZ0, LZZ0T, FM, FH, N )
!
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE PARKIND1  ,ONLY : JPRB
!RJ: added modi, after freeform conversion
      USE MODI_VSLOG
!
      IMPLICIT NONE
      INTEGER :: N
      REAL :: CMU(N),CTU(N),RIB(N),FCOR(N),ILMO(N)
      REAL :: FTEMP(N),FVAP(N),TA(N),QA(N),ZU(N),ZT(N),VA(N)
      REAL :: TG(N),QG(N),H(N),Z0(N),UE(N)
      REAL :: Z0T(N),LZZ0(N),LZZ0T(N)
      REAL :: fm(N),fh(N)
!
!Author
!          Y.Delage (Jul 1990)
!Revision
! 001      G. Pellerin (Jun 94) New function for unstable case
! 002      G. Pellerin (Jui 94) New formulation for stable case
! 003      B. Bilodeau (Nov 95) Replace VK by KARMAN
! 004      M. Desgagne (Dec 95) Add safety code in function ff
!                               and ensures that RIB is non zero
! 005      R. Sarrazin (Jan 96) Correction for H
! 006      C. Girard (Nov 95) - Diffuse T instead of Tv
! 007      G. Pellerin (Feb 96) Revised calculation for H (stable)
! 008      G. Pellerin (Feb 96) Remove corrective terms to CTU
! 009      Y. Delage and B. Bilodeau (Jul 97) - Cleanup
! 010      Y. Delage (Feb 98) - Addition of HMIN
! 011      D. Talbot and Y. Delage (Jan 02) -
!             Correct bug of zero divide by dg in loop 35
! 012      Y. Delage (Oct 03) - Set top of surface layer at ZU +Z0
!                   - Output UE instead of UE**2 and rename subroutine
!                   - Change iteration scheme for stable case
!                   - Introduce log-linear profile for near-neutral stable cases
!                   - set VAMIN inside flxsurf and initialise ILMO and H
!                   - Put stability functions into local functions via stabfunc.h
! 013      Y. Delage (Sep 04) - Input of wind and temperature/humidity
!                                at different levels
! 014      R. McTaggart-Cowan and B. Bilodeau (May 2006) -
!             Clean up stabfunc.h
! 015      L. Spacek (Dec 07) - Correction of the log-linear profile
!                               Double precision for rib calculations
!
!Object
!          to calculate surface layer transfer coefficients and fluxes
!
!Arguments
!
!          - Output -
! CMU      transfer coefficient of momentum times UE
! CTU      transfer coefficient of temperature times UE
! RIB      bulk Richardson number
! FTEMP    temperature flux
! FVAP     vapor flux
! ILMO     (1/length of Monin-Obukov)
! UE       friction velocity
! H        height of the boundary layer
! FM       momentum stability function
! FH       heat stability function
! LZZ0     log ((zu+z0)/z0)
! LZZ0T    log ((zt+z0)/z0t)
!
!          - Input -
! FCOR     Coriolis factor
! ZU       height of wind input (measured from model base at topo height + Z0)
! ZT       height of temperature and humidity input
! TA       potential temperature at ZT
! QA       specific humidity at ZT
! VA       wind speed at ZU
! TG       surface temperature
! QG       specific humidity at the surface
! Z0       roughness length for momentum      flux calculations
! Z0T      roughness length for heat/moisture flux calculations
! N        horizontal dimension
!
!
!RJ #include "surfcon.h"
!RJ       LOGICAL :: INIT
!     PHYSICAL CONSTANTS
       REAL,PARAMETER :: CPD      =.100546e+4        ! J K-1 kg-1    ! specific heat of dry air
       REAL,PARAMETER :: CPV      =.186946e+4        ! J K-1 kg-1    ! specific heat of water vapour
       REAL,PARAMETER :: RGASD    =.28705e+3         ! J K-1 kg-1    ! gas constant for dry air
       REAL,PARAMETER :: RGASV    =.46151e+3         ! J K-1 kg-1    ! gas constant for water vapour
       REAL,PARAMETER :: TRPL     =.27316e+3         ! K             ! triple point of water
       REAL,PARAMETER :: TCDK     =.27315e+3         !               ! conversion from kelvin to celsius
       REAL,PARAMETER :: RAUW     =.1e+4             !               ! density of liquid H2O
       REAL,PARAMETER :: EPS1     =.62194800221014   !               ! RGASD/RGASV
       REAL,PARAMETER :: EPS2     =.3780199778986    !               ! 1 - EPS1
       REAL,PARAMETER :: DELTA    =.6077686814144    !               ! 1/EPS1 - 1
       REAL,PARAMETER :: CAPPA    =.28549121795      !               ! RGASD/CPD
       REAL,PARAMETER :: TGL      =.27316e+3         ! K             ! ice temperature in the atmosphere
       REAL,PARAMETER :: CONSOL   =.1367e+4          ! W m-2         ! solar constant
       REAL,PARAMETER :: GRAV     =.980616e+1        ! M s-2         ! gravitational acceleration
       REAL,PARAMETER :: RAYT     =.637122e+7        ! M             ! mean radius of the earth
       REAL,PARAMETER :: STEFAN   =.566948e-7        ! J m-2 s-1 K-4 ! Stefan-Boltzmann constant
       REAL,PARAMETER :: PI       =.314159265359e+1  !               ! PI constant = ACOS(-1)
       REAL,PARAMETER :: OMEGA    =.7292e-4          ! s-1           ! angular speed of rotation of the earth
       REAL,PARAMETER :: KNAMS    =.514791           !               ! conversion from knots to m/s
       REAL,PARAMETER :: STLO     =.6628486583943e-3 ! K s2 m-2      ! Schuman-Newell Lapse Rate
       REAL,PARAMETER :: KARMAN   =.35               !               ! Von Karman constant
       REAL,PARAMETER :: RIC      =.2                !               ! Critical Richardson number
       REAL,PARAMETER :: CHLC     =.2501e+7          ! J kg-1        ! latent heat of condensation
       REAL,PARAMETER :: CHLF     =.334e+6           ! J kg-1        ! latent heat of fusion
       REAL,PARAMETER :: T1S      =.27316e+3         ! K             ! constant used to calculate L/Cp in fcn HTVOCP
       REAL,PARAMETER :: T2S      =.25816e+3         ! K             ! constant used to calculate L/Cp in fcn HTVOCP
       REAL,PARAMETER :: AW       =.3135012829948e+4 !               ! constant used to calculate L/Cp in fcn HTVOCP
       REAL,PARAMETER :: BW       =.2367075766316e+1 !               ! constant used to calculate L/Cp in fcn HTVOCP
       REAL,PARAMETER :: AI       =.2864887713087e+4 !               ! constant used to calculate L/Cp in fcn HTVOCP
       REAL,PARAMETER :: BI       =.166093131502     !               ! constant used to calculate L/Cp in fcn HTVOCP
       REAL,PARAMETER :: SLP      =.6666666666667e-1 !               ! constant used to calculate L/Cp in fcn HTVOCP

!RJ #include "consphy.h"
!     INITIALIZES THE CONSTANTS FOR THE COMMONS OF THE FLXSURF3 ROUTINE FROM
!     CANADIAN METEOROLOGICAL CENTER
      REAL,PARAMETER :: AS    = 12.
      REAL,PARAMETER :: ASX   = 5.
      REAL,PARAMETER :: CI    = 40.
      REAL,PARAMETER :: BS    = 1.0
      REAL,PARAMETER :: BETA  = 1.0
      REAL,PARAMETER :: FACTN = 1.2
      REAL,PARAMETER :: HMIN  = 30.
      REAL,PARAMETER :: ANGMAX= 0.85
      REAL,PARAMETER :: RAC3  = SQRT(3.)
!
!*
!
      INTEGER,PARAMETER :: JDBL=8
!
      INTEGER :: J
      INTEGER :: IT
      INTEGER,PARAMETER :: ITMAX = 3
      REAL,PARAMETER :: HMAX = 1500.0
      REAL,PARAMETER :: CORMIN = 0.7E-4
      REAL,PARAMETER :: EPSLN = 1.0e-05
      REAL,PARAMETER :: VAMIN = 0.1
      REAL :: CM,CT,ZP
      REAL :: F,G,DG
      REAL :: hi,HE,HS,unsl
      REAL(KIND=JDBL) :: DTHV,TVA,TVS
      REAL :: HL,U
      REAL :: CS,XX,XX0,YY,YY0
      REAL :: ZB,DD,ILMOX
      REAL :: DF,ZZ,betsasx
      REAL :: aa,bb,cc
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
      DF(ZZ)=(1-ZZ*hi)*sqrt(1+4*AS*BETA*unsl*ZZ/(1-ZZ*hi))
      CS=AS*2.5
      betsasx=1./asx
!
      IF (LHOOK) CALL DR_HOOK('FLXSURF3BX',0,ZHOOK_HANDLE)
!
      DO J=1,N
        LZZ0 (J)=1+ZU(J)/Z0(J)
        LZZ0T(J)=(ZT(J)+Z0(J))/Z0T(J)
      ENDDO
!
      call vslog(LZZ0T,LZZ0T,N)
      call vslog(LZZ0 ,LZZ0 ,N)
!
      DO J=1,N
!
!  CALCULATE THE RICHARDSON NUMBER
        ZP=ZU(J)**2/(ZT(J)+Z0(J)-Z0T(J))
        u=max(vamin,va(j))
        tva=(1.0_JDBL+DELTA*QA(J))*TA(J)
        tvs=(1.0_JDBL+DELTA*QG(J))*TG(J)
        dthv=tva-tvs
        RIB(J)=GRAV/(tvs+0.5_JDBL*dthv)*ZP*dthv/(u*u)
        if (rib(j)>=0.0_JDBL) rib(j) = max(rib(j), EPSLN)
        if (rib(j)<0.0_JDBL) rib(j) = min(rib(j),-EPSLN)
!
!  FIRST APPROXIMATION TO ILMO
        IF(RIB(J)>0.0_JDBL)  THEN
           FM(J)=LZZ0(J)+CS*RIB(J)/max(2*z0(j),1.0_JDBL)
           FH(J)=BETA*(LZZ0T(J)+CS*RIB(J))/                             &
     &          max(sqrt(z0(j)*z0t(j)),1.0_JDBL)
           ILMO(J)=RIB(J)*FM(J)*FM(J)/(ZP*FH(J))
           F=MAX(ABS(FCOR(J)),CORMIN)
           H(J)=BS*sqrt(KARMAN*u/(ILMO(J)*F*fm(j)))
        ELSE
           FM(J)=LZZ0(J)-min(0.7_JDBL+log(1-rib(j)),LZZ0(J)-1)
           FH(J)=BETA*(LZZ0T(J)-min(0.7_JDBL+log(1-rib(j)),LZZ0T(J)-1))
           ILMO(J)=RIB(J)*FM(J)*FM(J)/(ZP*FH(J))
        ENDIF
      ENDDO
!
! - - - - - - - - -  BEGINNING OF ITERATION LOOP - - - - - - - - - - -
      DO 35 IT=1,ITMAX
      DO 35 J=1,N
        u=max(vamin,va(j))
        ZP=ZU(J)**2/(ZT(J)+Z0(J)-Z0T(J))
        IF(RIB(J)>0.0_JDBL)  THEN
!----------------------------------------------------------------------
!  STABLE CASE
        ILMO(J)=max(EPSLN,ILMO(J))
        hl=(ZU(J)+10*Z0(J))*FACTN
        F=MAX(ABS(FCOR(J)),CORMIN)
        hs=BS*sqrt(KARMAN*u/(ILMO(J)*F*fm(j)))
        H(J)=MAX(HMIN,hs,hl,factn/(4*AS*BETA*ilmo(j)))
        hi=1/h(j)
!CDIR IEXPAND
        fm(J)=LZZ0(J)+psi(ZU(J)+Z0(J),hi,ilmo(j))-psi(Z0(J),hi,ilmo(j))
!CDIR IEXPAND
        fh(J)=BETA*(LZZ0T(J)+psi(ZT(J)+Z0(J),hi,ilmo(j))-psi(Z0T(J),hi, &
     &       ilmo(j)))
        unsl=ILMO(J)
        DG=-ZP*FH(J)/(FM(J)*FM(J))*(1+beta*(DF(ZT(J)+Z0(J))-DF(Z0T(J)))/&
     &     (2*FH(J))-(DF(ZU(J)+Z0(J))-DF(Z0(J)))/FM(J))
!----------------------------------------------------------------------
!  UNSTABLE CASE
      ELSE
        ILMO(J)=MIN(0.,ILMO(J))
!CDIR IEXPAND
        FM(J)=fmi(zu(j)+z0(j),z0 (j),lzz0 (j),ilmo(j),xx,xx0)
!CDIR IEXPAND
        FH(J)=fhi(zt(j)+z0(j),z0t(j),lzz0t(j),ilmo(j),yy,yy0)
         DG=-ZP*FH(J)/(FM(J)*FM(J))*(1+beta/FH(J)*(1/YY-1/YY0)-2/FM(J)* &
     &               (1/XX-1/XX0))
      ENDIF
!----------------------------------------------------------------------
      IF(IT<ITMAX) THEN
             G=RIB(J)-FH(J)/(FM(J)*FM(J))*ZP*ILMO(J)
             ILMO(J)=ILMO(J)-G/DG
      ENDIF
   35 CONTINUE
! - - - - - -  - - - END OF ITERATION LOOP - - - - - - - - - - - - - -
!
      DO 80 J=1,N
      u=max(vamin,va(j))
      if(asx<as) then
!----------------------------------------------------------------------
!  CALCULATE ILMO AND STABILITY FUNCTIONS FROM LOG-LINEAR PROFILE
!     (SOLUTION OF A QUADRATIC EQATION)
!
        zb=zu(j)/(zt(j)+z0(j)-z0t(j))
!  DISCRIMINANT
        dd=(beta*lzz0t(j)*zb)**2-4*rib(j)*asx*lzz0(j)*                  &
     &       (beta*lzz0t(j)*zb-lzz0(j))
        if(rib(j)>0.0_JDBL.and.rib(j)<betsasx.and.dd>=0.) then
!  COEFFICIENTS
           aa=asx*asx*rib(j)-asx
           bb=-beta*lzz0t(j)*zb+2*rib(j)*asx*lzz0(j)
           cc=rib(j)*lzz0(j)**2
!  SOLUTION
           if(bb>=0)then
              ilmox=(-bb-sqrt(dd))/(2*zu(j)*aa)
           else
              ilmox=2*cc/(zu(j)*(-bb+sqrt(dd)))
           endif
           if(ilmox<ilmo(j)) then
              ilmo(j)=ilmox
              fm(j)=lzz0(j)+asx*zu(j)*ilmox
              fh(j)=beta*lzz0t(j)+asx*(zt(j)+z0(j)-z0t(j))*ilmox
           endif
        endif
      endif
!----------------------------------------------------------------------
        CM=KARMAN/FM(J)
        CT=KARMAN/FH(J)
        UE(J)=u*CM
        CMU(J)=CM*UE(J)
        CTU(J)=CT*UE(J)
        if (rib(j)>0.0_JDBL) then
!          stable case
              H(J)=MIN(H(J),hmax)
        else
!          unstable case
              F=MAX(ABS(FCOR(J)),CORMIN)
              he=max(HMIN,0.3_JDBL*UE(J)/F)
              H(J)=MIN(he,hmax)
        endif
        FTEMP(J)=-CTU(J)*(TA(J)-TG(J))
        FVAP(J)=-CTU(J)*(QA(J)-QG(J))
   80 CONTINUE
      IF (LHOOK) CALL DR_HOOK('FLXSURF3BX',1,ZHOOK_HANDLE)
      CONTAINS
!RJ: inlining directly
!RJ #include "stabfunc2.h"
!
!   Internal function FMI
!   Stability function for momentum in the unstable regime (ilmo<0)
!   Reference: Delage Y. and Girard C. BLM 58 (19-31) Eq. 19
!
      FUNCTION FMI(Z2,Z02,LZZ02,ILMO2,X,X0)
      IMPLICIT NONE
!
      REAL              :: FMI
      REAL, INTENT(IN ) :: Z2,Z02,LZZ02,ILMO2
      REAL, INTENT(OUT) :: X,X0
!
      X =(1-CI*Z2 *BETA*ILMO2)**(0.16666666)
      X0=(1-CI*Z02*BETA*ILMO2)**(0.16666666)
      FMI=LZZ02+LOG((X0+1)**2*SQRT(X0**2-X0+1)*(X0**2+X0+1)**1.5        &
     &               /((X+1)**2*SQRT(X**2-X+1)*(X**2+X+1)**1.5))        &
     &              +RAC3*ATAN(RAC3*((X**2-1)*X0-(X0**2-1)*X)/          &
     &              ((X0**2-1)*(X**2-1)+3*X*X0))
!
      RETURN
      END FUNCTION FMI
!
!   Internal function FHI
!   Stability function for heat and moisture in the unstable regime (ilmo<0)
!   Reference: Delage Y. and Girard C. BLM 58 (19-31) Eq. 17
!
      FUNCTION FHI(Z2,Z0T2,LZZ0T2,ILMO2,Y,Y0)
      IMPLICIT NONE
!
      REAL              :: FHI
      REAL, INTENT(IN ) :: Z2,Z0T2,LZZ0T2,ILMO2
      REAL, INTENT(OUT) :: Y,Y0
!
      Y =(1-CI*Z2  *BETA*ILMO2)**(0.33333333)
      Y0=(1-CI*Z0T2*BETA*ILMO2)**(0.33333333)
      FHI=BETA*(LZZ0T2+1.5*LOG((Y0**2+Y0+1)/(Y**2+Y+1))+RAC3*           &
     &        ATAN(RAC3*2*(Y-Y0)/((2*Y0+1)*(2*Y+1)+3)))
!
      RETURN
      END FUNCTION FHI
!
!   Internal function psi
!   Stability function for momentum in the stable regime (unsl>0)
!   Reference :  Y. Delage, BLM, 82 (p23-48) (Eqs.33-37)
!
      FUNCTION PSI(Z2,HI2,ILMO2)
      IMPLICIT NONE
!
      REAL :: PSI
      REAL :: a,b,c,d
      REAL, INTENT(IN ) :: ILMO2,Z2,HI2
!
      d = 4*AS*BETA*ILMO2
      c = d*hi2 - hi2**2
      b = d - 2*hi2
      a = sqrt(1 + b*z2 - c*z2**2)
      psi = 0.5 * (a-z2*hi2-log(1+b*z2*0.5+a)-                          &
     &            b/(2*sqrt(c))*asin((b-2*c*z2)/d))
!
      RETURN
      END FUNCTION PSI
!
      END SUBROUTINE FLXSURF3BX
