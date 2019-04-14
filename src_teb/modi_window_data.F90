!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_WINDOW_DATA
INTERFACE
SUBROUTINE WINDOW_DATA(KI,XSHGC, XU_WIN, XALB_WIN, XABS_WIN, XUGG_WIN, XTRAN_WIN)
INTEGER,             INTENT(IN)  :: KI        ! number of points
REAL, DIMENSION(KI), INTENT(IN)  :: XSHGC    ! solar heat gain coef. of windows
REAL, DIMENSION(KI), INTENT(IN)  :: XU_WIN   ! window U-factor [K m W-2]
REAL, DIMENSION(KI), INTENT(OUT) :: XALB_WIN ! window albedo
REAL, DIMENSION(KI), INTENT(OUT) :: XABS_WIN ! window absortance
REAL, DIMENSION(KI), INTENT(OUT) :: XUGG_WIN ! window glass-to-glass U-factor [W m-2 K-1]
REAL, DIMENSION(KI), INTENT(OUT) :: XTRAN_WIN! window transmittance (-)
END SUBROUTINE WINDOW_DATA
END INTERFACE
END MODULE MODI_WINDOW_DATA
