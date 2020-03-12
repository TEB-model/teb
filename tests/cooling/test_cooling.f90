! MinimalDX version 0.1.1 (https://www.github.com/dmey/minimal-dx).
! Copyright 2018 D. Meyer and R. Raustad. Licensed under MIT.

program test_cooling
  !+ Program to compare MinimalDX against original EnergyPlus implementation.

  ! Add InitializePsychRoutines here to fix the allocatable array is already allocated error
  use Psychrometrics, only: InitializePsychRoutines
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

  integer         :: num_args, ioStatus = 0
  character(255)  :: InputFilePath, OutputFilePath, ModelName

  ! call InitializePsychRoutines here to fix the allocatable array is already allocated error
  call InitializePsychRoutines

  ! Check that the command line arguments are given and that are exactly 2
  num_args = command_argument_count()
  if ( num_args == 3 ) then
    call get_command_argument(1, InputFilePath)
    call get_command_argument(2, OutputFilePath)
    call get_command_argument(3, ModelName)
  else
    error stop "You must give the input and output file path as a command line argument"
  end if

  ! Open input file and skip first line -- it assumes that the file contains a header
  open(unit=10, file=trim(adjustl(InputFilePath)), status='old', action='read')
  read(10,*)

  ! Open output file and write header
  open (unit=20,file=trim(adjustl(OutputFilePath)), action="write", status="replace")
  write (20, '(1(A160))') "COP|1 TotalCoolingCapacity|W OutletTemperature|°C &
              &OutletHumRatio|kg/kg ElecCoolingPower|W TotalCoolingEnergyRate|W &
              &TotalSensibleHeatOut|W"

  do while (ioStatus == 0)
    read(10,*, iostat=ioStatus) OutdoorTDryBulb, OutdoorHumRatio, OutdoorPressure, &
         InletTDryBulb, InletHumRatio,  RatedCOP,                                  &
         RatedTotCap, SensibleCoolingLoad, RatedAirMassFlowRate

    if (ioStatus /= 0) then
      exit
    end if

    if (trim(adjustl(ModelName)) == 'EnergyPlus') then
      call SimCalcDoe2DXCoil(OutdoorTDryBulb, OutdoorHumRatio, OutdoorPressure,                &
                             InletTDryBulb, InletHumRatio,                                     &
                             RatedCOP, RatedTotCap, SensibleCoolingLoad, RatedAirMassFlowRate, &
                             COP, TotalCoolingCapacity,                                        &
                             OutletTemperature, OutletHumRatio,                                &
                             ElecCoolingPower, LatCoolingEnergyRate,                           &
                             TotalCoolingEnergyRate, TotalSensibleHeatOut)

    else if (trim(adjustl(ModelName)) == 'MinimalDX') then
      call SimMinimalDXCooling(OutdoorTDryBulb, OutdoorHumRatio, OutdoorPressure,                &
                               InletTDryBulb, InletHumRatio,                                     &
                               RatedCOP, RatedTotCap, SensibleCoolingLoad, RatedAirMassFlowRate, &
                               COP, TotalCoolingCapacity,                                        &
                               OutletTemperature, OutletHumRatio,                                &
                               ElecCoolingPower, LatCoolingEnergyRate,                           &
                               TotalCoolingEnergyRate, TotalSensibleHeatOut)


    else
      error stop 'The name is not recognized'
    end if

    ! Convert OutletTemperature from degC to K to avoid disc around 0 when comp with rel diff.
    write (20, '(7(ES24.15E2))') COP, TotalCoolingCapacity, OutletTemperature + 273.15, &
                OutletHumRatio, ElecCoolingPower, TotalCoolingEnergyRate, TotalSensibleHeatOut
    end do
  close(20)
  close(10)
end program test_cooling
