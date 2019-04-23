!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #############################################################
SUBROUTINE WINDOW_DATA(KI, B)
!     #############################################################
!
!
!!****  *WINDOW_DATA*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the radiative properties of the window
!         
!     
!!**  METHOD
!!     ------
!! uses U_WIN (window conductivity) and SHGC (window solar heat gain coef)
!! to compute the normal incident transmittance/reflectance
!! a factor to normal incident tran/ref is applied to take into account all
!! possible incidence between 18deg and 72deg
!!
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    MODD_CST
!!
!!      
!!    REFERENCE
!!    ---------
!!
!! Simple Window Model from EnergyPlus Engineering Reference V7, 2012, p. 217.
!!      
!!    AUTHOR
!!    ------
!!
!!      B. Bueno           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    11/10
!!      G. Pigeon & J. Le Bras  09/12 apply new factor to normal incident tran/ref
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_BEM_n, ONLY : BEM_t
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INTEGER,             INTENT(IN)  :: KI       ! number of points
TYPE(BEM_t), INTENT(INOUT) :: B
!
REAL, DIMENSION(KI) :: ZR_IW      ! interior film resistance under winter conditions
REAL, DIMENSION(KI) :: ZR_OW      ! exterior fild resistance uder winter conditions
REAL, DIMENSION(KI) :: ZR_IS      ! interior film resistance under summer conditions
REAL, DIMENSION(KI) :: ZR_OS      ! exterior fild resistance uder summer conditions
REAL, DIMENSION(KI) :: ZFRAC      ! inward flowing fraction
REAL, DIMENSION(KI) :: ZTRAN_WIN_NORM   ! Window transmittance for normal incident beam
REAL, DIMENSION(KI) :: ZALB_WIN_NORM   ! Window reflectance for normal incident beam
REAL,DIMENSION(KI)  :: ZFAC_TRAN ! transmittance factor of the window to the normal incidence
REAL,DIMENSION(KI)  :: ZFAC_ALB  ! reflectance (albedo) factor of the window to normal incidence
REAL(KIND=JPRB) :: ZHOOK_HANDLE

REAL:: TA, RA, TBCD, RBCD, TD, RD, TE, RE, TFGHI, RFGHI, TFH, RFH, TJ, RJ
!

!* initialization and determine transmittance and reflectance factor of the window

!Curve AT temp Et        RA = R temp Et
TA = 0.9237767331    ; RA = 1.6911265959
TBCD = 0.8802058937  ; RBCD = 1.5239470129
TD = 0.893430955     ; RD = 1.0888507586
TE = 0.8743096978    ; RE = 1.5494187256
TFH = 0.8424427652   ; RFH = 1.0983010317
TFGHI = 0.8327695761 ; RFGHI = 1.3316554732
TJ = 0.7889001341    ; RJ = 1.0837577691

WHERE (B%XU_WIN <= 1.42) !from 3 building BLD, TOUR and MI-2005
   !ZONE 1   
   WHERE (B%XSHGC>=0.45)   
      ZFAC_TRAN = TE 
         ZFAC_ALB = RE 
   END WHERE
   !ZONE 2
   WHERE (B%XSHGC<0.45 .AND. B%XSHGC>0.35)
      ZFAC_TRAN = TJ+(B%XSHGC-0.35)*(TE-TJ)/(0.45-0.35) 
         ZFAC_ALB = RJ+(B%XSHGC-0.35)*(RE-RJ)/(0.45-0.35)    
   END WHERE
   !ZONE 3
   WHERE(B%XSHGC<=0.35)
      ZFAC_TRAN = TJ  
         ZFAC_ALB = RJ
   END WHERE   
END WHERE


WHERE (B%XU_WIN>1.42 .AND. B%XU_WIN<1.70)
   !ZONE 4   
   WHERE (B%XSHGC>=0.55)
      ZFAC_TRAN = TE 
         ZFAC_ALB = RE
   END WHERE
   !ZONE 5
   WHERE (B%XSHGC>0.5 .AND. B%XSHGC<0.55)
      ZFAC_TRAN = TE+(B%XU_WIN-1.42)*((TFGHI+(B%XSHGC-0.5)*(TE-TFGHI)/(0.55-0.50))-TE)/(1.70-1.42)  
         ZFAC_ALB = RE+(B%XU_WIN-1.42)*((RFGHI+(B%XSHGC-0.5)*(RE-RFGHI)/(0.55-0.50))-RE)/(1.70-1.42)
   END WHERE
   !ZONE 6
   WHERE (B%XSHGC>=0.45 .AND. B%XSHGC<=0.5)
      ZFAC_TRAN =  TE+(B%XU_WIN-1.42)*(TFGHI-TE)/(1.70-1.42)
         ZFAC_ALB = RE+(B%XU_WIN-1.42)*(RFGHI-RE)/(1.70-1.42)
   END WHERE
   !ZONE 7
   WHERE (B%XSHGC>0.35 .AND. B%XSHGC<0.45)
      ZFAC_TRAN = (TJ+(B%XSHGC-0.35)*(TE-TJ)/(0.45-0.35))+ &
              (B%XU_WIN-1.42)*(TFGHI-(TJ+(B%XSHGC-0.35)*(TE-TJ)/(0.45-0.35)))/(1.70-1.42)  
         
      ZFAC_ALB = (RJ+(B%XSHGC-0.35)*(RE-RJ)/(0.45-0.35))+&
              (B%XU_WIN-1.42)*(RFGHI-(RJ+(B%XSHGC-0.35)*(RE-RJ)/(0.45-0.35)))/(1.70-1.42) 
   END WHERE
   !ZONE 8
   WHERE (B%XSHGC>=0.30 .AND. B%XSHGC<=0.35)
      ZFAC_TRAN =  TJ+(B%XU_WIN-1.42)*(TFGHI-TJ)/(1.70-1.42)
         ZFAC_ALB = RJ+(B%XU_WIN-1.42)*(RFGHI-RJ)/(1.70-1.42)
   END WHERE
   !ZONE 9
   WHERE (B%XSHGC>0.25 .AND. B%XSHGC<0.35)
      ZFAC_TRAN = TJ+(B%XU_WIN-1.42)*((TFH+(B%XSHGC-0.25)*(TJ-TFH)/(0.30-0.25))-TJ)/(1.70-1.42)  
         ZFAC_ALB = RJ+(B%XU_WIN-1.42)*((RFH+(B%XSHGC-0.25)*(RJ-RFH)/(0.30-0.25))-RJ)/(1.70-1.42)
   END WHERE
   !ZONE 10
   WHERE (B%XSHGC<=0.25)
      ZFAC_TRAN =  TJ+(B%XU_WIN-1.42)*(TFH-TJ)/(1.70-1.42)
         ZFAC_ALB = RJ+(B%XU_WIN-1.42)*(RFH-RJ)/(1.70-1.42)
   END WHERE


END WHERE


WHERE (B%XU_WIN>=1.70 .AND. B%XU_WIN<=3.41)
   !ZONE 11   
   WHERE (B%XSHGC>=0.55)
      ZFAC_TRAN = TE  
         ZFAC_ALB = RE
   END WHERE
   !ZONE 12
   WHERE (B%XSHGC > 0.5 .AND. B%XSHGC<0.55)
      ZFAC_TRAN = TFGHI+(B%XSHGC-0.50)*(TE-TFGHI)/(0.55-0.50) 
         ZFAC_ALB = RFGHI+(B%XSHGC-0.50)*(RE-RFGHI)/(0.55-0.50) 
   END WHERE    
   !ZONE 13
   WHERE (B%XSHGC>=0.30 .AND. B%XSHGC<=0.50)
      ZFAC_TRAN = TFGHI  
         ZFAC_ALB = RFGHI
   END WHERE

   !ZONE 14
   WHERE (B%XSHGC > 0.25 .AND. B%XSHGC<0.30)
      ZFAC_TRAN = TFH+(B%XSHGC-0.25)*(TFGHI-TFH)/(0.30-0.25) 
         ZFAC_ALB = RFGHI+(B%XSHGC-0.25)*(RFGHI-RFH)/(0.30-0.25) 
   END WHERE 
   !ZONE 15
   WHERE (B%XSHGC<=0.25)
      ZFAC_TRAN = TFH  
         ZFAC_ALB = RFH
   END WHERE

END WHERE

WHERE (B%XU_WIN>3.41 .AND. B%XU_WIN<4.54)
   !ZONE 16
   WHERE (B%XSHGC>=0.65)
      ZFAC_TRAN =  TE+(B%XU_WIN-3.41)*(TA-TE)/(4.54-3.41)
         ZFAC_ALB = RE+(B%XU_WIN-3.41)*(RA-RE)/(4.54-3.41)
   END WHERE
   !ZONE 17
   WHERE (B%XSHGC>0.60 .AND. B%XSHGC<0.65)
      ZFAC_TRAN = TE+(B%XU_WIN-3.41)*((TBCD+(B%XSHGC-0.60)*(TA-TBCD)/(0.65-0.60))-TE)/(4.54-3.41)  
         ZFAC_ALB = RE+(B%XU_WIN-3.41)*((RBCD+(B%XSHGC-0.60)*(RA-RBCD)/(0.65-0.60))-RE)/(4.54-3.41)
   END WHERE
   !ZONE 18
   WHERE (B%XSHGC>=0.55 .AND. B%XSHGC<=0.60)
      ZFAC_TRAN =  TE+(B%XU_WIN-3.41)*(TBCD-TE)/(4.54-3.41)
         ZFAC_ALB = RE+(B%XU_WIN-3.41)*(RBCD-RE)/(4.54-3.41)
   END WHERE
   !ZONE 19
   WHERE (B%XSHGC>0.50 .AND. B%XSHGC<0.55)
      ZFAC_TRAN = (TFGHI+(B%XSHGC-0.50)*(TE-TFGHI)/(0.55-0.50)) + &
                            (B%XU_WIN-3.41)*(TBCD-(TFGHI+(B%XSHGC-0.50)*(TE-TFGHI)/(0.55-0.50)))/(4.54-3.41)  
         
      ZFAC_ALB = (RFGHI+(B%XSHGC-0.50)*(RE-RFGHI)/(0.55-0.50)) + &
                           (B%XU_WIN-3.41)*(RBCD-(RFGHI+(B%XSHGC-0.50)*(RE-RFGHI)/(0.55-0.50)))/(4.54-3.41) 
   END WHERE
   !ZONE 20
   WHERE (B%XSHGC>=0.45 .AND. B%XSHGC<=0.50)
      ZFAC_TRAN =  TFGHI+(B%XU_WIN-3.41)*(TBCD-TFGHI)/(4.54-3.41)
         ZFAC_ALB = RFGHI+(B%XU_WIN-3.41)*(RBCD-RFGHI)/(4.54-3.41)
   END WHERE
   !ZONE 21
   WHERE (B%XSHGC>0.30 .AND. B%XSHGC<0.45)
      ZFAC_TRAN = TFGHI+(B%XU_WIN-3.41)*((TA+(B%XSHGC-0.30)*(TBCD-TA)/(0.45-0.30))-TFGHI)/(4.54-3.41)  
         ZFAC_ALB = RFGHI+(B%XU_WIN-3.41)*((RA+(B%XSHGC-0.30)*(RBCD-RA)/(0.45-0.30))-RFGHI)/(4.54-3.41)
   END WHERE
   !ZONE 22
   WHERE (B%XSHGC>0.25 .AND. B%XSHGC<0.30)
      ZFAC_TRAN = (TFH+(B%XSHGC-0.25)*(TFGHI-TFH)/(0.30-0.25)) + &
                  (B%XU_WIN-3.41)*(TD-(TFH+(B%XSHGC-0.30)*(TFGHI-TFH)/(0.30-0.25)))/(4.54-3.41)  
         ZFAC_ALB = (RFH+(B%XSHGC-0.25)*(RFGHI-RFH)/(0.30-0.25)) + &
                    (B%XU_WIN-3.41)*(RD-(RFH+(B%XSHGC-0.30)*(RFGHI-RFH)/(0.30-0.25)))/(4.54-3.41)
   END WHERE
   !ZONE 23
   WHERE (B%XSHGC<=0.25)
      ZFAC_TRAN =  TFH+(B%XU_WIN-3.41)*(TD-TFH)/(4.54-3.41)
         ZFAC_ALB = RFH+(B%XU_WIN-3.41)*(RD-RFH)/(4.54-3.41)
   END WHERE

END WHERE


WHERE (B%XU_WIN>=4.54)
   !ZONE 24
   WHERE (B%XSHGC>=0.65)
      ZFAC_TRAN = TA  
         ZFAC_ALB = RA
   END WHERE
   !ZONE 25
   WHERE (B%XSHGC >= 0.60 .AND. B%XSHGC<=0.65)
      ZFAC_TRAN = TBCD+(B%XSHGC-0.60)*(TA-TBCD)/(0.65-0.60) 
         ZFAC_ALB = RBCD+(B%XSHGC-0.60)*(RA-RBCD)/(0.65-0.60) 
   END WHERE 
   !ZONE 26
   WHERE (B%XSHGC>=0.45 .AND. B%XSHGC<=0.60)
      ZFAC_TRAN = TBCD  
         ZFAC_ALB = RBCD
   END WHERE
   !ZONE 27
   WHERE (B%XSHGC >= 0.30 .AND. B%XSHGC<=0.45)
      ZFAC_TRAN = TD+(B%XSHGC-0.30)*(TBCD-TD)/(0.45-0.30) 
         ZFAC_ALB = RD+(B%XSHGC-0.30)*(RBCD-RD)/(0.45-0.30) 
   END WHERE
   !ZONE 28
   WHERE (B%XSHGC<=0.30)
      ZFAC_TRAN = TD  
         ZFAC_ALB = RD
   END WHERE
   

END WHERE

!*      A.     Determine glass-to-glass U-factor
!              ---------------------------------
!
IF (LHOOK) CALL DR_HOOK('WINDOW_DATA',0,ZHOOK_HANDLE)
!
WHERE (B%XU_WIN(:) < 5.85)
  ZR_IW(:) = 1. / (0.359073*log(B%XU_WIN(:)) + 6.949915)
ELSEWHERE
  ZR_IW(:) = 1. / (1.788041*    B%XU_WIN(:)  - 2.886625)
END WHERE
!
ZR_OW(:) = 1. / (0.025342*B%XU_WIN(:) + 29.163853)
!
B%XUGG_WIN(:) = 1./(1./B%XU_WIN(:) - ZR_IW(:) - ZR_OW(:))
!
!*      B.     Determine layer solar transmittance
!              -----------------------------------
!
WHERE (B%XU_WIN(:) > 4.5) 
  WHERE (B%XSHGC(:) < 0.7206)
    ZTRAN_WIN_NORM(:) = 0.939998 * B%XSHGC(:)**2 + 0.20332 * B%XSHGC(:)
  ELSE WHERE
    ZTRAN_WIN_NORM(:) = 1.30415 * B%XSHGC(:) - 0.30515
  END WHERE
END WHERE 

WHERE (B%XU_WIN(:) < 3.4)
  WHERE (B%XSHGC(:) <= 0.15)
    ZTRAN_WIN_NORM(:) = 0.41040 * B%XSHGC(:)
  ELSE WHERE
    ZTRAN_WIN_NORM(:) = 0.085775 * B%XSHGC(:)**2 + 0.963954 * B%XSHGC(:) - 0.084958
  END WHERE
END WHERE

WHERE (B%XU_WIN(:) >= 3.4 .AND. B%XU_WIN(:) <=4.5)
   WHERE(B%XSHGC(:) <= 0.15)
      ZTRAN_WIN_NORM(:) = 0.5* (0.939998 * B%XSHGC(:)**2 + (0.20332+0.41040) * B%XSHGC(:))
   ELSE WHERE
      WHERE(B%XSHGC(:) > 0.15 .AND. B%XSHGC(:) < 0.7206)
         ZTRAN_WIN_NORM(:) = 0.5 * &
                 ((0.939998+0.085775) * B%XSHGC(:)**2 + (0.20332 + 0.963954) * B%XSHGC(:) - 0.084958)
      ELSE WHERE
         ZTRAN_WIN_NORM(:) = 0.5 * &
                 (0.085775 * B%XSHGC(:)**2 + (0.963954+1.30415) * B%XSHGC(:) - (0.084958+0.30515))
      END WHERE
   END WHERE
END WHERE
!
!
!*      C.     Determine layer solar reflectance
!              ---------------------------------
!
WHERE (B%XU_WIN(:) > 4.5)
  ZR_IS(:) = 1. / ( 29.436546*(B%XSHGC(:)-ZTRAN_WIN_NORM(:))**3         &
                  - 21.943415*(B%XSHGC(:)-ZTRAN_WIN_NORM(:))**2         &
                  + 9.945872 *(B%XSHGC(:)-ZTRAN_WIN_NORM(:)) + 7.426151 )
!
  ZR_OS (:) = 1./ ( 2.225824*(B%XSHGC(:)-ZTRAN_WIN_NORM(:)) + 20.57708 )
ELSE WHERE
   WHERE(B%XU_WIN(:) >= 3.4)
      ZR_IS(:) = 0.5 * (                                                &
                1. / ( 29.436546*(B%XSHGC(:)-ZTRAN_WIN_NORM(:))**3        &
                  - 21.943415*(B%XSHGC(:)-ZTRAN_WIN_NORM(:))**2           &
                  + 9.945872 *(B%XSHGC(:)-ZTRAN_WIN_NORM(:)) + 7.426151 ) &
              + 1./ (199.8208128*(B%XSHGC(:)-ZTRAN_WIN_NORM(:))**3        &
                  - 90.639733*(B%XSHGC(:)-ZTRAN_WIN_NORM(:))**2           &
                  + 19.737055*(B%XSHGC(:)-ZTRAN_WIN_NORM(:)) + 6.766575) ) 
!
      ZR_OS (:) = 0.5 * ( &
                  1./ ( 2.225824*(B%XSHGC(:)-ZTRAN_WIN_NORM(:)) + 20.57708 ) &
                + 1./ ( 5.763355*(B%XSHGC(:)-ZTRAN_WIN_NORM(:)) + 20.541528 ) )
!
   ELSE WHERE
      ZR_IS (:) = 1./ (199.8208128*(B%XSHGC(:)-ZTRAN_WIN_NORM(:))**3       &
                  - 90.639733*(B%XSHGC(:)-ZTRAN_WIN_NORM(:))**2        &
                  + 19.737055*(B%XSHGC(:)-ZTRAN_WIN_NORM(:)) + 6.766575)
!
      ZR_OS(:) = 1./ ( 5.763355*(B%XSHGC(:)-ZTRAN_WIN_NORM(:)) + 20.541528 )
   END WHERE
!
END WHERE
!
ZFRAC(:) = (ZR_OS(:) + 0.5*(1./B%XUGG_WIN(:))) /               &
           (ZR_OS(:) + 1./B%XUGG_WIN(:) + ZR_IS(:))
!
ZALB_WIN_NORM(:) = 1. - ZTRAN_WIN_NORM(:) - (B%XSHGC(:)-ZTRAN_WIN_NORM(:)) / ZFRAC(:)
!
!*      D.     Determine layer solar absortance
!              --------------------------------
!

B%XTRAN_WIN(:) = ZTRAN_WIN_NORM(:) * ZFAC_TRAN 
B%XALB_WIN(:) = ZALB_WIN_NORM(:) * ZFAC_ALB
!
B%XABS_WIN(:) = 1. - B%XALB_WIN(:) - B%XTRAN_WIN(:)
!
IF (LHOOK) CALL DR_HOOK('WINDOW_DATA',1,ZHOOK_HANDLE)
!
END SUBROUTINE WINDOW_DATA
