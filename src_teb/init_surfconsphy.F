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
       SUBROUTINE INIT_SURFCONSPHY

*     INITIALIZES THE CONSTANTS FOR THE COMMONS OF THE FLXSURF3 ROUTINE FROM
*     CANADIAN METEOROLOGICAL CENTER
*
#include "surfcon.cdk"
#include "consphy.cdk"

*     CONSTANTS RELATED TO VERTICAL DIFFUSION 

      AS    = 12.
      ASX   = 5.
      CI    = 40.
      BS    = 1.0
      BETA  = 1.0
      FACTN = 1.2
      HMIN  = 30.
      ANGMAX= 0.85
      RAC3  = SQRT(3.)

*     PHYSICAL CONSTANTS

       CPD      =.100546e+4        ! J K-1 kg-1    ! specific heat of dry air
       CPV      =.186946e+4        ! J K-1 kg-1    ! specific heat of water vapour
       RGASD    =.28705e+3         ! J K-1 kg-1    ! gas constant for dry air
       RGASV    =.46151e+3         ! J K-1 kg-1    ! gas constant for water vapour
       TRPL     =.27316e+3         ! K             ! triple point of water
       TCDK     =.27315e+3         !               ! conversion from kelvin to celsius
       RAUW     =.1e+4             !               ! density of liquid H2O
       EPS1     =.62194800221014   !               ! RGASD/RGASV
       EPS2     =.3780199778986    !               ! 1 - EPS1
       DELTA    =.6077686814144    !               ! 1/EPS1 - 1
       CAPPA    =.28549121795      !               ! RGASD/CPD
       TGL      =.27316e+3         ! K             ! ice temperature in the atmosphere
       CONSOL   =.1367e+4          ! W m-2         ! solar constant
       GRAV     =.980616e+1        ! M s-2         ! gravitational acceleration
       RAYT     =.637122e+7        ! M             ! mean radius of the earth
       STEFAN   =.566948e-7        ! J m-2 s-1 K-4 ! Stefan-Boltzmann constant
       PI       =.314159265359e+1  !               ! PI constant = ACOS(-1)
       OMEGA    =.7292e-4          ! s-1           ! angular speed of rotation of the earth
       KNAMS    =.514791           !               ! conversion from knots to m/s
       STLO     =.6628486583943e-3 ! K s2 m-2      ! Schuman-Newell Lapse Rate
       KARMAN   =.35               !               ! Von Karman constant
       RIC      =.2                !               ! Critical Richardson number
       CHLC     =.2501e+7          ! J kg-1        ! latent heat of condensation
       CHLF     =.334e+6           ! J kg-1        ! latent heat of fusion
       T1S      =.27316e+3         ! K             ! constant used to calculate L/Cp in fcn HTVOCP
       T2S      =.25816e+3         ! K             ! constant used to calculate L/Cp in fcn HTVOCP
       AW       =.3135012829948e+4 !               ! constant used to calculate L/Cp in fcn HTVOCP
       BW       =.2367075766316e+1 !               ! constant used to calculate L/Cp in fcn HTVOCP
       AI       =.2864887713087e+4 !               ! constant used to calculate L/Cp in fcn HTVOCP
       BI       =.166093131502     !               ! constant used to calculate L/Cp in fcn HTVOCP
       SLP      =.6666666666667e-1 !               ! constant used to calculate L/Cp in fcn HTVOCP

       END SUBROUTINE INIT_SURFCONSPHY
