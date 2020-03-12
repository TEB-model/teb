! MinimalDX version 0.1.1 (https://www.github.com/dmey/minimal-dx).
! Copyright 2018 D. Meyer and R. Raustad. Licensed under MIT.

module PsychroWrapper
  !+ Wraps PsychroLib functions.

use psychrolib, only: SetUnitSystem, SI, GetMoistAirDensity, GetTWetBulbFromHumRatio, &
                      GetMoistAirEnthalpy, GetHumRatioFromEnthalpyAndTDryBulb,        &
                      GetTDryBulbFromEnthalpyAndHumRatio, GetRelHumFromHumRatio,      &
                      GetHumRatioFromTDewPoint
use EPlusPsychro, only: PsyTsatFnHPb

implicit none

private

public :: GetMoistAirDensity
public :: GetTWetBulbFromHumRatio
public :: GetMoistAirEnthalpy
public :: GetHumRatioFromEnthalpyAndTDryBulb
public :: GetTDryBulbFromEnthalpyAndHumRatio
public :: GetRelHumFromHumRatio
public :: GetHumRatioFromTDewPoint
public :: PsyTsatFnHPb
public :: InitPsychrometrics

contains
    subroutine InitPsychrometrics()
      call SetUnitSystem(SI)
    end subroutine InitPsychrometrics
end module PsychroWrapper
