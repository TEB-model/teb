! MinimalDX version 0.1.2 (https://www.github.com/dmey/minimal-dx).
! Copyright 2018-2020 D. Meyer and R. Raustad. Licensed under MIT.

module MinimalDXCoolingDriver
  !+ Contains a simplified EnergyPlus subroutine for simulating the performance of a DX cooling coil.

  implicit none

  private
  public SimMinimalDXCooling

  contains
  subroutine SimMinimalDXCooling(OutdoorTDryBulb, OutdoorHumRatio, OutdoorPressure,                & ! I
                                 InletTDryBulb, InletHumRatio,                                     & ! I
                                 RatedCOP, RatedTotCap, SensibleCoolingLoad, RatedAirMassFlowRate, & ! I
                                 COP, TotalCoolingCapacity,                                        & ! O
                                 OutletTemperature, OutletHumRatio,                                & ! O
                                 ElecCoolingPower, LatCoolingEnergyRate,                           & ! O
                                 TotalCoolingEnergyRate, TotalSensibleHeatOut)                       ! O

    !+ Simplified EnergyPlus subroutine for calculating the performance of a DX cooling coil.
    !+ This is the subroutine to call from your program.

    use MinimalDXCooling, only: CalcMinimalDXCooling


    ! Subroutine arguments
    real, intent(in)    :: OutdoorTDryBulb
        !+ Outdoor dry bulb air temperature `[°C]`
    real, intent(in)    :: OutdoorHumRatio
        !+ Outdoor air humidity ratio `[kgH₂O kgAIR⁻¹]`
    real, intent(in)    :: OutdoorPressure
        !+ Outdoor barometric pressure `[Pa]`
    real, intent(in)    :: InletTDryBulb
        !+ Indoor (inlet) dry bulb air temperature `[°C]`
    real, intent(in)    :: InletHumRatio
        !+ Indoor (inlet) air humidity ratio `[kgH₂O kgAIR⁻¹]`
    real, intent(in)    :: RatedCOP
        !+ Rated Coefficient Of Performance (COP) `[1]`
    real, intent(in)    :: RatedTotCap
        !+ Rated (total) system capacity `[W]`
    real, intent(in)    :: SensibleCoolingLoad
        !+ Building sensible load to be met `[W]`
    real, intent(in)    :: RatedAirMassFlowRate
        !+ rated air mass flow rate `[kg s⁻¹]`
    real, intent(out)   :: COP
    !+ Actual (calculated) Coefficient Of Performance (COP) `[1]`
    real, intent(out)   :: TotalCoolingCapacity
        !+ Actual (calculated) total system capacity `[W]`
    real, intent(out)   :: OutletTemperature
        !+ Actual (calculated) outlet air dry bulb temperature existing the cooling coil `[°C]`
    real, intent(out)   :: OutletHumRatio
        !+ Actual (calculated) outlet air humidity ratio existing the cooling coil `[kgH₂O kgAIR⁻¹]`
    real, intent(out)   :: ElecCoolingPower
        !+ Calculated electrical power consumed by the DX unit `[W]`
    real, intent(out)   :: LatCoolingEnergyRate
        ! Total latent cooling energy rate extracted by the coil from the indoor environment `[J kg⁻¹]`
    real, intent(out)   :: TotalCoolingEnergyRate
        !+ Total cooling power of the DX unit (energy rate extracted by DX unit from the indoor environment) `[W]`
    real, intent(out)   :: TotalSensibleHeatOut
        !+ Total power rejected by the evaporator into the outdoor environment
        !+ i.e. TotalCoolingEnergyRate + ElecCoolingPower `[W]`

    ! Local variables
    real :: PartLoadRatio
    real :: SensCoolingEnergyRate
        !+ Sensible cooling power used to calculate the PLR. This is the maximum amount of sensible heat rate that the coil
        !+ is capable of extracting from the indoor environment for the specified conditions. `[W]`

    ! Get SensibleCoolingLoad and SensCoolingEnergyRates to calculate actual PartLoadRatio
    PartLoadRatio = 1.

    call CalcMinimalDXCooling(OutdoorTDryBulb, OutdoorHumRatio, OutdoorPressure,                & ! I
                              InletTDryBulb, InletHumRatio,                                     & ! I
                              RatedCOP, RatedTotCap, PartLoadRatio, RatedAirMassFlowRate,       & ! I
                              OutletTemperature, OutletHumRatio,                                & ! O
                              ElecCoolingPower, SensCoolingEnergyRate, LatCoolingEnergyRate,    & ! O
                              TotalCoolingEnergyRate, TotalSensibleHeatOut)

    ! Calculate PartLoadRatio from previous call
    PartLoadRatio = SensibleCoolingLoad / SensCoolingEnergyRate

    ! Keep PartLoadRatio  bounded between 0 and 1
    if ( PartLoadRatio < 0. ) PartLoadRatio = 0.
    if ( PartLoadRatio > 1. ) PartLoadRatio = 1.

    ! Call CalcDoe2DXCoil again with correct part load ratio
    call CalcMinimalDXCooling(OutdoorTDryBulb, OutdoorHumRatio, OutdoorPressure,                & ! I
                              InletTDryBulb, InletHumRatio,                                     & ! I
                              RatedCOP, RatedTotCap, PartLoadRatio, RatedAirMassFlowRate,       & ! I
                              OutletTemperature, OutletHumRatio,                                & ! O
                              ElecCoolingPower, SensCoolingEnergyRate, LatCoolingEnergyRate,    & ! O
                              TotalCoolingEnergyRate, TotalSensibleHeatOut)

    COP = TotalCoolingEnergyRate / ElecCoolingPower
    TotalCoolingCapacity = TotalCoolingEnergyRate / PartLoadRatio


  end subroutine SimMinimalDXCooling
end module MinimalDXCoolingDriver