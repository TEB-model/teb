! MinimalDX version 0.1.2 (https://www.github.com/dmey/minimal-dx).
! Copyright 2018-2020 D. Meyer and R. Raustad. Licensed under MIT.

program make_test_data_cooling
  !+ Generates test data for testing the cooling coils.

  ! Using Energy Plus psychrometric routines
  use Psychrometrics, only: PsyWFnTdbRhPb, PsyRhoAirFnPbTdbW

  implicit none

  !+ @note
  !+ The number one ('1') next to the variable description denotes that the variable is dimensionless - e.g.   COP `[1]`

  real :: OutdoorTDryBulb
      !+ Outdoor dry bulb air temperature `[°C]`
  real :: OutdoorRelHum
      !+ Outdoor relative humidity `[1]`
  real :: OutdoorHumRatio
      !+ Outdoor air humidity ratio `[kgH₂O kgAIR⁻¹]`
  real :: OutdoorPressure
      !+ Outdoor barometric pressure `[Pa]`
  real :: InletTDryBulb
      !+ Indoor (inlet) dry bulb air temperature `[°C]`
  real :: InletRelHum
      !+ Indoor (inlet) relative humidity `[1]`
  real :: InletHumRatio
      !+ Indoor (inlet) air humidity ratio `[kgH₂O kgAIR⁻¹]`
  real :: RatedCOP
      !+ Rated Coefficient Of Performance (COP) `[1]`
  real :: RatedTotCap
      !+ Rated (total) system capacity `[W]`
  real :: SensibleCoolingLoad
      !+ Building sensible load to be met `[W]`
  real :: RatedAirMassFlowRate
      !+ Rated air mass flow rate `[kg s⁻¹]
  integer ::  IdxOutdoorTDryBulb
      !+ Index for OutdoorTDryBulb `[1]`
  integer ::  IdxOutdoorRelHum
      !+ Index for OutdoorRelHum `[1]`
  integer ::  IdxOutdoorPressure
      !+ Index for OutdoorPressure `[1]`
  integer ::  IdxInletTDryBulb
      !+ Index for InletTDryBulb `[1]`
  integer ::  IdxInletRelHum
      !+ Index for InletRelHum `[1]`
  integer ::  IdxRatedCOP
      !+ Index for RatedCOP `[1]`
  integer ::  IdxRatedTotCap
      !+ Index for RatedTotCap `[1]`
  integer ::  IdxSensibleCoolingLoad
      !+ Index for SensibleCoolingLoad `[1]`
  integer ::  IterationCounter
      !+ Iteration counter `[1]`
  integer ::  num_args
      !+ Number of command line arguments `[1]`
  character(255) ::  DataSetPath
      !+ Name of dataset including path `[1]`

  IterationCounter = 0

  ! Check that the command line argument is given and that it is exactly 1
  num_args = command_argument_count()
  if ( num_args == 1 ) then
    call get_command_argument(1, DataSetPath)
  else
    error stop "You must give the filename path of the data set you want to save as a command line argument"
  end if


  open (unit=10, file=trim(adjustl(DataSetPath)), action="write", status="replace")
  write (10, '(A220)') "OutdoorTDryBulb|°C OutdoorHumRatio|kg/kg OutdoorPressure|Pa &
                        &InletTDryBulb|°C InletHumRatio|kg/kg RatedCOP|1 &
                        &RatedTotCap|W SensibleCoolingLoad|W RatedAirMassFlowRate|m3/(kg*s) &
                        &OutdoorRelHum|[0-1] InletRelHum|[0-1]"

  do IdxOutdoorTDryBulb = 20, 40, 2
      do IdxOutdoorRelHum = 0, 100, 5
          do IdxOutdoorPressure = 60000, 120000, 10000
              ! The internal temperature is always greater than the outdoor temperature
              do IdxInletTDryBulb = 20, 40, 2
                  do IdxInletRelHum = 0, 100, 5
                      do IdxRatedCOP = 4, 4, 4
                          do IdxRatedTotCap = 4000, 4000, 4000
                              ! SensibleCoolingLoad scaled between 0% and 100% of RatedTotCap with steps of 2000 W
                              do IdxSensibleCoolingLoad   = IdxRatedTotCap * 0, IdxRatedTotCap * 1, IdxRatedTotCap * 0.5

                                  OutdoorTDryBulb     = DBLE(IdxOutdoorTDryBulb)
                                  OutdoorRelHum       = DBLE(IdxOutdoorRelHum) / 100
                                  OutdoorPressure     = DBLE(IdxOutdoorPressure)
                                  InletTDryBulb       = DBLE(IdxInletTDryBulb)
                                  InletRelHum         = DBLE(IdxInletRelHum) / 100
                                  RatedCOP            = DBLE(IdxRatedCOP)
                                  RatedTotCap         = DBLE(IdxRatedTotCap)
                                  SensibleCoolingLoad = DBLE(IdxSensibleCoolingLoad)

                                  OutdoorHumRatio = PsyWFnTdbRhPb(OutdoorTDryBulb, OutdoorRelHum, OutdoorPressure)
                                  InletHumRatio = PsyWFnTdbRhPb(InletTDryBulb, InletRelHum, OutdoorPressure)
                                  RatedAirMassFlowRate = 0.00005 * RatedTotCap * &
                                                          PsyRhoAirFnPbTdbW(OutdoorPressure, InletTDryBulb, InletHumRatio)

                                  write (10, '(11(ES19.10E2))') OutdoorTDryBulb, OutdoorHumRatio, OutdoorPressure,        &
                                                                  InletTDryBulb, InletHumRatio, RatedCOP,                   &
                                                                  RatedTotCap, SensibleCoolingLoad, RatedAirMassFlowRate,   &
                                                                  OutdoorRelHum, InletRelHum
                                  IterationCounter = IterationCounter + 1
                              end do
                          end do
                      end do
                  end do
              end do
          end do
      end do
  end do

  print *, 'Total number of iterations: ', IterationCounter

end program make_test_data_cooling