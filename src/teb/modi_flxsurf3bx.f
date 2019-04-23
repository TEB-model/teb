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
!     ######spl
      MODULE MODI_FLXSURF3BX
      INTERFACE
      SUBROUTINE FLXSURF3BX(CMU, CTU, RIB, FTEMP, FVAP, ILMO,
     1                   UE, FCOR, TA , QA , ZU, ZT, VA,
     1                   TG , QG , H , Z0 , Z0T,
     1                   LZZ0, LZZ0T, FM, FH, N )
!
      INTEGER N
      REAL CMU(N),CTU(N),RIB(N),FCOR(N),ILMO(N)
      REAL FTEMP(N),FVAP(N),TA(N),QA(N),ZU(N),ZT(N),VA(N)
      REAL TG(N),QG(N),H(N),Z0(N),UE(N)
      REAL Z0T(N),LZZ0(N),LZZ0T(N)
      REAL fm(N),fh(N)
      END
      END INTERFACE
      END MODULE MODI_FLXSURF3BX
