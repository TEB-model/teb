! MinimalDX version 0.1.1 (https://www.github.com/dmey/minimal-dx).
! Copyright 2018 D. Meyer and R. Raustad. Licensed under MIT.

module PsychroWrapper
  !+ Wraps PsychroLib functions.

use Psychrometrics, only: PsyRhoAirFnPbTdbW, & ! Becomes GetMoistAirDensity
                          GetTWetBulbFromHumRatio => PsyTwbFnTdbWPb,        &
                          GetMoistAirEnthalpy => PsyHFnTdbW,                &
                          PsyWFnTdbH, & ! Becomes GetHumRatioFromEnthalpyAndTDryBulb
                          GetTDryBulbFromEnthalpyAndHumRatio => PsyTdbFnHW, &
                          GetRelHumFromHumRatio => PsyRhFnTdbWPb,           &
                          GetHumRatioFromTDewPoint => PsyWFnTdpPb,          &
                          PsyTsatFnHPb

implicit none

private

public :: InitPsychrometrics
public :: GetMoistAirDensity
public :: GetTWetBulbFromHumRatio
public :: GetMoistAirEnthalpy
public :: GetHumRatioFromEnthalpyAndTDryBulb
public :: GetTDryBulbFromEnthalpyAndHumRatio
public :: GetRelHumFromHumRatio
public :: GetHumRatioFromTDewPoint
public :: PsyTsatFnHPb

contains
  subroutine InitPsychrometrics()
  end subroutine InitPsychrometrics

  function GetMoistAirDensity(InletTDryBulb, InletHumRatio, OutdoorPressure) result(MoistAirDensity)
    real :: InletTDryBulb, InletHumRatio, OutdoorPressure, MoistAirDensity
    MoistAirDensity = PsyRhoAirFnPbTdbW(OutdoorPressure, InletTDryBulb, InletHumRatio)
  end function GetMoistAirDensity

  function GetHumRatioFromEnthalpyAndTDryBulb(hADP, tADP) result(HumRatio)
    real :: hADP, tADP, HumRatio
    HumRatio = PsyWFnTdbH(tADP, hADP)
  end function GetHumRatioFromEnthalpyAndTDryBulb
end module PsychroWrapper
