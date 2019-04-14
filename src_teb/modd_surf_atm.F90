!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ####################
      MODULE MODD_SURF_ATM
!     ####################
!
!!****  *MODD_SURF_ATM - declaration of surface ATM
!!
!!    PURPOSE
!!    -------
!     Declaration of surface parameters
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	P. Le Moigne *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       10/2007
!
!*       0.   DECLARATIONS
!             ------------
!
!!      B.Decharme     04/2009 Add flag used to Read/Write precipitation forcing from/into the restart file for ARPEGE/ALADIN run
!!      B.Decharme     08/2009 Add flag used to know if you use SURFEX in the Earth System Model
!!      B. Decharme    09/2012 New wind implicitation key option
!
IMPLICIT NONE
!
!-----------------------------------------------------------------------------------------------------
REAL       :: XCISMIN  ! minimum wind shear
REAL       :: XVMODMIN ! minimum wind speed
LOGICAL    :: LALDTHRES! activate aladin threshold for wind
!
LOGICAL    :: LDRAG_COEF_ARP ! activate aladin formulation for Cd and Ch
LOGICAL    :: LALDZ0H
!
LOGICAL    :: LNOSOF   ! No parameterization of Subgrid Orography effects on atmospheric Forcing  
LOGICAL    :: LVERTSHIFT  ! vertical shift from atmospheric orography to surface orography
!
LOGICAL    :: LRW_PRECIP       ! flag used to Read/Write precipitation forcing from/into the restart file for ARPEGE/ALADIN run
LOGICAL    :: LSAVE_PRECIP     ! used with LRW_PRECIP
!
LOGICAL    :: LVZIUSTAR0_ARP   ! activate aladin formulation for zoh over sea
LOGICAL    :: LRRGUST_ARP      ! activate aladin formulation for CD CH, CDN, correction due to moist gustiness
LOGICAL    :: LCPL_ARP         ! activate aladin formulation for Cp and L
LOGICAL    :: LQVNPLUS         ! An option for the resolution of the surface temperature equation
!
LOGICAL    :: LCPL_ESM         ! Earth System Model Coupling key
!
 CHARACTER(LEN=3) :: CIMPLICIT_WIND ! wind implicitation option
!                                  ! 'OLD' = direct
!                                  ! 'NEW' = Taylor serie, order 1
!
REAL       :: XEDB
REAL       :: XEDC
REAL       :: XEDD
REAL       :: XEDK
REAL       :: XUSURIC
REAL       :: XUSURID
REAL       :: XUSURICL
REAL       :: XVCHRNK
REAL       :: XVZ0CM
REAL       :: XRIMAX
REAL       :: XDELTA_MAX ! Maximum fraction of the foliage covered by intercepted water for high vegetation
!
REAL       :: XWINDMIN ! minimum wind speed (canopy)
!
REAL       :: XRZHZ0M
REAL       :: XVZIUSTAR0
REAL       :: XRRSCALE
REAL       :: XRRGAMMA
REAL       :: XUTILGUST
!-----------------------------------------------------------------------------------------------------
!
END MODULE MODD_SURF_ATM
