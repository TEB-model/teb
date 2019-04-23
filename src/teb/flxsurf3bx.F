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
***S/P  FLXSURF3
*
      SUBROUTINE FLXSURF3BX(CMU, CTU, RIB, FTEMP, FVAP, ILMO,
     X                    UE, FCOR, TA , QA , ZU, ZT, VA,
     Y                    TG , QG , H , Z0 , Z0T,
     %                    LZZ0, LZZ0T, FM, FH, N )

#include "impnone.cdk"
      INTEGER N
      REAL CMU(N),CTU(N),RIB(N),FCOR(N),ILMO(N)
      REAL FTEMP(N),FVAP(N),TA(N),QA(N),ZU(N),ZT(N),VA(N)
      REAL TG(N),QG(N),H(N),Z0(N),UE(N)
      REAL Z0T(N),LZZ0(N),LZZ0T(N)
      REAL fm(N),fh(N)
*
*Author
*          Y.Delage (Jul 1990)
*Revision
* 001      G. Pellerin (Jun 94) New function for unstable case
* 002      G. Pellerin (Jui 94) New formulation for stable case
* 003      B. Bilodeau (Nov 95) Replace VK by KARMAN
* 004      M. Desgagne (Dec 95) Add safety code in function ff
*                               and ensures that RIB is non zero
* 005      R. Sarrazin (Jan 96) Correction for H
* 006      C. Girard (Nov 95) - Diffuse T instead of Tv
* 007      G. Pellerin (Feb 96) Revised calculation for H (stable)
* 008      G. Pellerin (Feb 96) Remove corrective terms to CTU
* 009      Y. Delage and B. Bilodeau (Jul 97) - Cleanup
* 010      Y. Delage (Feb 98) - Addition of HMIN
* 011      D. Talbot and Y. Delage (Jan 02) - 
*             Correct bug of zero divide by dg in loop 35
* 012      Y. Delage (Oct 03) - Set top of surface layer at ZU +Z0
*                   - Output UE instead of UE**2 and rename subroutine
*                   - Change iteration scheme for stable case
*                   - Introduce log-linear profile for near-neutral stable cases
*                   - set VAMIN inside flxsurf and initialise ILMO and H
*                   - Put stability functions into local functions via stabfunc.cdk
* 013      Y. Delage (Sep 04) - Input of wind and temperature/humidity
*                                at different levels
* 014      R. McTaggart-Cowan and B. Bilodeau (May 2006) -
*             Clean up stabfunc.cdk
* 015      L. Spacek (Dec 07) - Correction of the log-linear profile
*                               Double precision for rib calculations
*
*Object
*          to calculate surface layer transfer coefficients and fluxes
*
*Arguments
*
*          - Output -
* CMU      transfer coefficient of momentum times UE
* CTU      transfer coefficient of temperature times UE
* RIB      bulk Richardson number
* FTEMP    temperature flux
* FVAP     vapor flux
* ILMO     (1/length of Monin-Obukov)
* UE       friction velocity 
* H        height of the boundary layer
* FM       momentum stability function
* FH       heat stability function
* LZZ0     log ((zu+z0)/z0)
* LZZ0T    log ((zt+z0)/z0t)
*
*          - Input -
* FCOR     Coriolis factor
* ZU       height of wind input (measured from model base at topo height + Z0)
* ZT       height of temperature and humidity input
* TA       potential temperature at ZT
* QA       specific humidity at ZT
* VA       wind speed at ZU
* TG       surface temperature
* QG       specific humidity at the surface
* Z0       roughness length for momentum      flux calculations
* Z0T      roughness length for heat/moisture flux calculations
* N        horizontal dimension
*
*
#include "surfcon.cdk"
#include "consphy.cdk"
*
**
*
      INTEGER J
      INTEGER IT,ITMAX
      REAL HMAX,CORMIN,EPSLN
      REAL CM,CT,ZP
      REAL F,G,DG
      REAL hi,HE,HS,unsl
      REAL*8 DTHV,TVA,TVS
      REAL HL,VAMIN,U
      REAL CS,XX,XX0,YY,YY0
      REAL ZB,DD,ILMOX
      REAL DF,ZZ,betsasx
      REAL aa,bb,cc
      SAVE HMAX,CORMIN,EPSLN
      SAVE ITMAX
      SAVE VAMIN
*
      DATA CORMIN, HMAX /0.7E-4 ,  1500.0/
      DATA ITMAX / 3  /
      DATA VAMIN /0.1/
      DATA EPSLN / 1.0e-05 /
*
      DF(ZZ)=(1-ZZ*hi)*sqrt(1+4*AS*BETA*unsl*ZZ/(1-ZZ*hi))
      CS=AS*2.5
      betsasx=1./asx
*
      DO J=1,N
        LZZ0 (J)=1+ZU(J)/Z0(J)
        LZZ0T(J)=(ZT(J)+Z0(J))/Z0T(J)
      ENDDO
*
      call vslog(LZZ0T,LZZ0T,N)
      call vslog(LZZ0 ,LZZ0 ,N)
*
      DO J=1,N
*
*  CALCULATE THE RICHARDSON NUMBER
        ZP=ZU(J)**2/(ZT(J)+Z0(J)-Z0T(J))
        u=max(vamin,va(j))
        tva=(1.d0+DELTA*QA(J))*TA(J)
        tvs=(1.d0+DELTA*QG(J))*TG(J)
        dthv=tva-tvs
        RIB(J)=GRAV/(tvs+0.5*dthv)*ZP*dthv/(u*u)
        if (rib(j).ge.0.0) rib(j) = max(rib(j), EPSLN)
        if (rib(j).lt.0.0) rib(j) = min(rib(j),-EPSLN)
*
*  FIRST APPROXIMATION TO ILMO
        IF(RIB(J).GT.0.)  THEN
           FM(J)=LZZ0(J)+CS*RIB(J)/max(2*z0(j),1.0)
           FH(J)=BETA*(LZZ0T(J)+CS*RIB(J))/
     1           max(sqrt(z0(j)*z0t(j)),1.0)
           ILMO(J)=RIB(J)*FM(J)*FM(J)/(ZP*FH(J))
           F=MAX(ABS(FCOR(J)),CORMIN)
           H(J)=BS*sqrt(KARMAN*u/(ILMO(J)*F*fm(j)))
        ELSE
           FM(J)=LZZ0(J)-min(0.7+log(1-rib(j)),LZZ0(J)-1)
           FH(J)=BETA*(LZZ0T(J)-min(0.7+log(1-rib(j)),LZZ0T(J)-1))
           ILMO(J)=RIB(J)*FM(J)*FM(J)/(ZP*FH(J))
        ENDIF
      ENDDO
*
* - - - - - - - - -  BEGINNING OF ITERATION LOOP - - - - - - - - - - -
      DO 35 IT=1,ITMAX
      DO 35 J=1,N
        u=max(vamin,va(j))
        ZP=ZU(J)**2/(ZT(J)+Z0(J)-Z0T(J))
        IF(RIB(J).GT.0.)  THEN
*----------------------------------------------------------------------
*  STABLE CASE
        ILMO(J)=max(EPSLN,ILMO(J))
        hl=(ZU(J)+10*Z0(J))*FACTN
        F=MAX(ABS(FCOR(J)),CORMIN)
        hs=BS*sqrt(KARMAN*u/(ILMO(J)*F*fm(j)))
        H(J)=MAX(HMIN,hs,hl,factn/(4*AS*BETA*ilmo(j)))
        hi=1/h(j)
*CDIR IEXPAND
        fm(J)=LZZ0(J)+psi(ZU(J)+Z0(J),hi,ilmo(j))-psi(Z0(J),hi,ilmo(j))
*CDIR IEXPAND
        fh(J)=BETA*(LZZ0T(J)+psi(ZT(J)+Z0(J),hi,ilmo(j))-psi(Z0T(J),hi,
     1        ilmo(j)))
        unsl=ILMO(J)
        DG=-ZP*FH(J)/(FM(J)*FM(J))*(1+beta*(DF(ZT(J)+Z0(J))-DF(Z0T(J)))/
     1      (2*FH(J))-(DF(ZU(J)+Z0(J))-DF(Z0(J)))/FM(J))
*----------------------------------------------------------------------
*  UNSTABLE CASE
      ELSE
        ILMO(J)=MIN(0.,ILMO(J))
*CDIR IEXPAND
        FM(J)=fmi(zu(j)+z0(j),z0 (j),lzz0 (j),ilmo(j),xx,xx0)
*CDIR IEXPAND
        FH(J)=fhi(zt(j)+z0(j),z0t(j),lzz0t(j),ilmo(j),yy,yy0)
         DG=-ZP*FH(J)/(FM(J)*FM(J))*(1+beta/FH(J)*(1/YY-1/YY0)-2/FM(J)*
     %                (1/XX-1/XX0))
      ENDIF
*----------------------------------------------------------------------
      IF(IT.LT.ITMAX) THEN
             G=RIB(J)-FH(J)/(FM(J)*FM(J))*ZP*ILMO(J)
             ILMO(J)=ILMO(J)-G/DG
      ENDIF
   35 CONTINUE
* - - - - - -  - - - END OF ITERATION LOOP - - - - - - - - - - - - - -
*
      DO 80 J=1,N
      u=max(vamin,va(j))
      if(asx.lt.as) then
*----------------------------------------------------------------------
*  CALCULATE ILMO AND STABILITY FUNCTIONS FROM LOG-LINEAR PROFILE
*     (SOLUTION OF A QUADRATIC EQATION)
*
        zb=zu(j)/(zt(j)+z0(j)-z0t(j))
*  DISCRIMINANT
        dd=(beta*lzz0t(j)*zb)**2-4*rib(j)*asx*lzz0(j)*
     1        (beta*lzz0t(j)*zb-lzz0(j))
        if(rib(j).gt.0..and.rib(j).lt.betsasx.and.dd.ge.0.) then           
*  COEFFICIENTS
           aa=asx*asx*rib(j)-asx
           bb=-beta*lzz0t(j)*zb+2*rib(j)*asx*lzz0(j)
           cc=rib(j)*lzz0(j)**2
*  SOLUTION
           if(bb>=0)then
              ilmox=(-bb-sqrt(dd))
     1          /(2*zu(j)*aa)
           else
              ilmox=2*cc/(zu(j)*(-bb+sqrt(dd)))
           endif
           if(ilmox.lt.ilmo(j)) then
              ilmo(j)=ilmox
              fm(j)=lzz0(j)+asx*zu(j)*ilmox
              fh(j)=beta*lzz0t(j)+asx*(zt(j)+z0(j)-z0t(j))*ilmox
           endif
        endif
      endif
*----------------------------------------------------------------------
        CM=KARMAN/FM(J)
        CT=KARMAN/FH(J)
        UE(J)=u*CM
        CMU(J)=CM*UE(J)
        CTU(J)=CT*UE(J)
        if (rib(j).gt.0.0) then
*          stable case
              H(J)=MIN(H(J),hmax)
        else
*          unstable case
              F=MAX(ABS(FCOR(J)),CORMIN)
              he=max(HMIN,0.3*UE(J)/F)
              H(J)=MIN(he,hmax)
        endif
        FTEMP(J)=-CTU(J)*(TA(J)-TG(J))
        FVAP(J)=-CTU(J)*(QA(J)-QG(J))
   80 CONTINUE
      RETURN
      CONTAINS
#include "stabfunc2.cdk"
      END
