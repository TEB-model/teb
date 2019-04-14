!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 1998-2013 Meteo-France
! This is part of the TEB software governed by the CeCILL-C licence version 1.
! See LICENCE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt for details.
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt
! http://www.cecill.info/licences/Licence_CeCILL-C_V1-fr.txt
! The CeCILL-C licence is compatible with L-GPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MODI_ROOF_IMPL_COEF
INTERFACE
SUBROUTINE ROOF_IMPL_COEF(PTSTEP, KROOF_LAYER, PD_ROOF, PTC_ROOF, PHC_ROOF, PT_ROOF, PTDEEP_A,PTDEEP_B)
REAL              ,   INTENT(IN)  :: PTSTEP      ! time step
INTEGER           ,   INTENT(IN)  :: KROOF_LAYER ! number of roof layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_ROOF   ! thickness            of each layer
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_ROOF  ! thermal conductivity of each layer
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_ROOF  ! heat capacity        of each layer
REAL, DIMENSION(:,:), INTENT(IN)  :: PT_ROOF   ! temperature          of each layer
REAL, DIMENSION(:),   INTENT(OUT) :: PTDEEP_A, PTDEEP_B     
END SUBROUTINE ROOF_IMPL_COEF
END INTERFACE
END MODULE MODI_ROOF_IMPL_COEF
