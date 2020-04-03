! MinimalDX version 0.1.2 (https://www.github.com/dmey/minimal-dx).
! Copyright 2018-2020 D. Meyer and R. Raustad. Licensed under MIT.

program  test_cooling_point
  !+ Program to compare MinimalDX against original EnergyPlus implementation.

  ! Add InitializePsychRoutines here to fix the allocatable array is already allocated error
  use Psychrometrics, only: InitializePsychRoutines, PsyRhoAirFnPbTdbW, PsyWFnTdbRhPb
  ! EnergyPlus DX Coil model (original implementation)
  use EPlusWrapperCooling, only : SimCalcDoe2DXCoil
  ! DX Coil model (simplified) using EnergyPlus psychrometric routines
  use MinimalDXCoolingDriver, only: SimMinimalDXCooling

  implicit none

  ! The number one ('1') next to the variable description denotes that the variable is dimensionless - e.g.   COP [1]
  real :: OutdoorTDryBulb
      !+ Outdoor dry bulb air temperature `[°C]`
  real :: OutdoorHumRatio
      !+ Outdoor air humidity ratio `[kgH₂O kgAIR⁻¹]`
  real :: OutdoorPressure
      !+ Outdoor barometric pressure `[Pa]`
  real :: InletTDryBulb
      !+ Indoor (inlet) dry bulb air temperature `[°C]`
  real :: InletHumRatio
      !+ Indoor (inlet) air humidity ratio `[kgH₂O kgAIR⁻¹]`
  real :: RatedCOP
      !+ Rated Coefficient Of Performance (COP) `[1]`
  real :: RatedTotCap
      !+ Rated (total) system capacity `[W]`
  real :: SensibleCoolingLoad
  !+ Building sensible load to be met `[W]`
  real :: RatedAirMassFlowRate
      !+ Rated air mass flow rate `[kg s⁻¹]`
  real :: COP
      !+ Actual (calculated) Coefficient Of Performance (COP) `[1]`
  real :: TotalCoolingCapacity
      !+ Actual (calculated) total system capacity `[W]`
  real :: OutletTemperature
      !+ Actual (calculated) outlet air dry bulb temperature existing the cooling coil `[°C]`
  real :: OutletHumRatio
      !+ Actual (calculated) outlet air humidity ratio existing the cooling coil `[kgH₂O kgAIR⁻¹]`
  real :: ElecCoolingPower
      !+ Calculated electrical power consumed by the DX unit `[W]`
  real :: LatCoolingEnergyRate
      ! Total latent cooling energy rate extracted by the coil from the indoor environment `[J kg⁻¹]`
  real :: TotalCoolingEnergyRate
      !+ Total cooling power of the DX unit (energy rate extracted by DX unit from the indoor environment) `[W]`
  real :: TotalSensibleHeatOut
      !+ Total power rejected by the evaporator into the outdoor environment
      !+ i.e. TotalCoolingEnergyRate + ElecCoolingPower `[W]`

  ! call InitializePsychRoutines here to fix the allocatable array is already allocated error
  call InitializePsychRoutines

  OutdoorTDryBulb = 20.
  OutdoorPressure = 90000.
  OutdoorHumRatio = PsyWFnTdbRhPb(OutdoorTDryBulb, 0.5, OutdoorPressure)
  InletTDryBulb = 23.
  InletHumRatio = PsyWFnTdbRhPb(InletTDryBulb, 0.5, OutdoorPressure)
  RatedCOP = 3.
  RatedTotCap = 4000.
  SensibleCoolingLoad = 100.
  RatedAirMassFlowRate = 0.00005 * RatedTotCap * PsyRhoAirFnPbTdbW(OutdoorPressure, InletTDryBulb, InletHumRatio)

  call SimCalcDoe2DXCoil(OutdoorTDryBulb, OutdoorHumRatio, OutdoorPressure,                 &
                          InletTDryBulb, InletHumRatio,                                     &
                          RatedCOP, RatedTotCap, SensibleCoolingLoad, RatedAirMassFlowRate, &
                          COP, TotalCoolingCapacity,                                        &
                          OutletTemperature, OutletHumRatio,                                &
                          ElecCoolingPower, LatCoolingEnergyRate,    &
                          TotalCoolingEnergyRate, TotalSensibleHeatOut)

  call SimMinimalDXCooling(OutdoorTDryBulb, OutdoorHumRatio, OutdoorPressure,                 &
                            InletTDryBulb, InletHumRatio,                                     &
                            RatedCOP, RatedTotCap, SensibleCoolingLoad, RatedAirMassFlowRate, &
                            COP, TotalCoolingCapacity,                                        &
                            OutletTemperature, OutletHumRatio,                                &
                            ElecCoolingPower, LatCoolingEnergyRate,    &
                            TotalCoolingEnergyRate, TotalSensibleHeatOut)

end program test_cooling_point
