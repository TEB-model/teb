! MinimalDX version 0.1.1 (https://www.github.com/dmey/minimal-dx).
! Copyright 2018 D. Meyer and R. Raustad. Licensed under MIT.
!
! Description
! This module currently contains one function (SimpleFanOnOff) for simulating a simple on/off fan.
! Modified by D. Meyer and R. Raustad (2018).
!
! Copyright notice
! EnergyPlus, Copyright (c) 1996-2019, The Board of Trustees of the University of Illinois,
! The Regents of the University of California, through Lawrence Berkeley National Laboratory
! (subject to receipt of any required approvals from the U.S. Dept. of Energy), Oak Ridge
! National Laboratory, managed by UT-Battelle, Alliance for Sustainable Energy, LLC, and other
! contributors. All rights reserved.
!
! NOTICE: This Software was developed under funding from the U.S. Department of Energy and the
! U.S. Government consequently retains certain rights. As such, the U.S. Government has been
! granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable,
! worldwide license in the Software to reproduce, distribute copies to the public, prepare
! derivative works, and perform publicly and display publicly, and to permit others to do so.
!
! Redistribution and use in source and binary forms, with or without modification, are permitted
! provided that the following conditions are met:
!
! (1) Redistributions of source code must retain the above copyright notice, this list of
!     conditions and the following disclaimer.
!
! (2) Redistributions in binary form must reproduce the above copyright notice, this list of
!     conditions and the following disclaimer in the documentation and/or other materials
!     provided with the distribution.
!
! (3) Neither the name of the University of California, Lawrence Berkeley National Laboratory,
!     the University of Illinois, U.S. Dept. of Energy nor the names of its contributors may be
!     used to endorse or promote products derived from this software without specific prior
!     written permission.
!
! (4) Use of EnergyPlus(TM) Name. If Licensee (i) distributes the software in stand-alone form
!     without changes from the version obtained under this License, or (ii) Licensee makes a
!     reference solely to the software portion of its product, Licensee must refer to the
!     software as "EnergyPlus version X" software, where "X" is the version number Licensee
!     obtained under this License and may not use a different name for the software. Except as
!     specifically required in this Section (4), Licensee shall not use in a company name, a
!     product name, in advertising, publicity, or other promotional activities any name, trade
!     name, trademark, logo, or other designation of "EnergyPlus", "E+", "e+" or confusingly
!     similar designation, without the U.S. Department of Energy's prior written consent.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
! AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
! OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.

module MinimalDXFan
  !+ This module currently contains one function (SimpleFanOnOff) for simulating
  !+ a simple on/off fan. Modified by D. Meyer and R. Raustad (2018).

  implicit none

  contains
  pure function GetOnOffFan(Mode, MotEff, FanPower, MotInAirFrac, InletAirEnthalpy, AirMassFlowRate) result(OutletAirEnthalpy)
    !+ Simplified version of SimOnOffFan subroutine in EnergyPlus
    !+ Given the mode of operation (on or off), fan motor efficiency, power of the fan, fraction of motor heat entering air stream
    !+ moist air enthaply of the air entering the fan, and mass flow rate, it returns the moist air enthaply of the air
    !+ after it has passed through the fan. It assumes sensible heating process only - i.e. it does not change the moisture
    !+ in the air - and the mass flow rate across the component remains the same.
    !+ Original Credits: NREL Energy Plus, Shirey, R. Raustad - FSEC, Brent Griffith,  Chandan Sharma, Rongpeng Zhang
    !+ Reference: ASHRAE HVAC 2 Toolkit, page 2-3 (FANSIM)
    !+ https://github.com/NREL/EnergyPlus/blob/d37252156cb0eef0cb9b1af5ce7dcd7423011649/src/EnergyPlus/Fans.cc#L1967-L1989
    !+ This fan does not change the moisture or Mass Flow across the component


    ! Function arguments:
    integer,intent(in)  ::  Mode                ! Mode of operation: 1 for on, 0 for off                [1]
    real,   intent(in)  ::  MotEff              ! Fan motor efficiency                                  [1]
    real,   intent(in)  ::  FanPower            ! Power of the fan to be simulated                      [W]
    real,   intent(in)  ::  MotInAirFrac        ! Fraction of motor heat entering air stream            [1]
    real,   intent(in)  ::  InletAirEnthalpy    ! Moist air enthaply of the air before entering the fan [J/kg]
    real,   intent(in)  ::  AirMassFlowRate     ! Mass flow rate through the Fan being simulated        [kg/Sec]
    real                ::  OutletAirEnthalpy   ! Moist air enthaply of the air after entering the fan  [J/kg]

    ! Local variables
    real    ::  FanShaftPower
    real    ::  PowerLossToAir

    if (Mode == 1 ) then
      ! Power delivered to shaft
      FanShaftPower = MotEff * FanPower
      PowerLossToAir = FanShaftPower + (FanPower - FanShaftPower) * MotInAirFrac
      OutletAirEnthalpy = InletAirEnthalpy + PowerLossToAir / AirMassFlowRate

    else if (Mode == 0 ) then
      ! Fan is off and not operating no power consumed and mass flow rate.
      OutletAirEnthalpy = InletAirEnthalpy

    else
      error stop 'Fan Mode must either be 1 for on or 0 for off'
    end if
  end function GetOnOffFan
end module MinimalDXFan
