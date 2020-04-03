! MinimalDX version 0.1.2 (https://www.github.com/dmey/minimal-dx).
! Copyright 2018-2020 D. Meyer and R. Raustad. Licensed under MIT.

module EPlusWrapperCooling
  !+ Initialise CalcDoe2DXCoil from EnergyPlusFortran and CalcDXCoolingCoil

  implicit none

  private
  public :: SimCalcDoe2DXCoil

  contains

  subroutine  SimCalcDoe2DXCoil(OutdoorTDryBulb, OutdoorHumRatio, OutdoorPressure,                & ! I
                                InletTDryBulb, InletHumRatio,                                     & ! I
                                RatedCOP, RatedTotCap, SensibleCoolingLoad, RatedAirMassFlowRate, & ! I
                                COP, TotalCoolingCapacity,                                        & ! O
                                OutletTemperature, OutletHumRatio,                                & ! O
                                ElecCoolingPower, LatCoolingEnergyRate,    & ! O
                                TotalCoolingEnergyRate, TotalSensibleHeatOut)

    use DXCoils,        only:  CalcDoe2DXCoil, DXCoil, DXCoilFullLoadOutAirTemp,                 &
                                DXCoilFullLoadOutAirHumRat, DXCoilOutletTemp, DXCoilOutletHumRat, &
                                DXCoilPartLoadRatio, DXCoilFanOpMode
    use DataHeatBalance,only:  HeatReclaimDXCoil
    use DataEnvironment,only:  OutBaroPress, OutDryBulbTemp, OutHumRat, OutWetBulbTemp
    use Psychrometrics, only:  PsyTwbFnTdbWPb,  PsyHFnTdbW
    use CurveManager,   only:  NumCurves, PerfCurve
    use ScheduleManager,only:  Schedule

    integer, parameter :: r64=kind(1.0d0)

    ! Subroutine arguments
    ! The number one ('1') next to the variable description denotes that the variable is dimensionless - e.g.   COP [1]
    real(r64), intent(in)    :: OutdoorTDryBulb
        !+ Outdoor dry bulb air temperature `[°C]`
    real(r64), intent(in)    :: OutdoorHumRatio
        !+ Outdoor air humidity ratio `[kgH₂O kgAIR⁻¹]`
    real(r64), intent(in)    :: OutdoorPressure
        !+ Outdoor barometric pressure `[Pa]`
    real(r64), intent(in)    :: InletTDryBulb
        !+ Indoor (inlet) dry bulb air temperature `[°C]`
    real(r64), intent(in)    :: InletHumRatio
        !+ Indoor (inlet) air humidity ratio `[kgH₂O kgAIR⁻¹]`
    real(r64), intent(in)    :: RatedCOP
        !+ Rated Coefficient Of Performance (COP) `[1]`
    real(r64), intent(in)    :: RatedTotCap
        !+ Rated (total) system capacity `[W]`
    real(r64), intent(in)    :: SensibleCoolingLoad
    !+ Building sensible load to be met `[W]`
    real(r64), intent(in)    :: RatedAirMassFlowRate
        !+ Rated air mass flow rate `[kg s⁻¹]`
    real(r64), intent(out)   :: COP
        !+ Actual (calculated) Coefficient Of Performance (COP) `[1]`
    real(r64), intent(out)   :: TotalCoolingCapacity
        !+ Actual (calculated) total system capacity `[W]`
    real(r64), intent(out)   :: OutletTemperature
        !+ Actual (calculated) outlet air dry bulb temperature existing the cooling coil `[°C]`
    real(r64), intent(out)   :: OutletHumRatio
        !+ Actual (calculated) outlet air humidity ratio existing the cooling coil `[kgH₂O kgAIR⁻¹]`
    real(r64), intent(out)   :: ElecCoolingPower
        !+ Calculated electrical power consumed by the DX unit `[W]`
    real(r64), intent(out)   :: LatCoolingEnergyRate
    ! Total latent cooling energy rate extracted by the coil from the indoor environment `[J kg⁻¹]`
    real(r64), intent(out)   :: TotalCoolingEnergyRate
        !+ Total cooling power of the DX unit (energy rate extracted by DX unit from the indoor environment) `[W]`
    real(r64), intent(out)   :: TotalSensibleHeatOut
        !+ Total power rejected by the evaporator into the outdoor environment
        !+ i.e. TotalCoolingEnergyRate + ElecCoolingPower `[W]`

    ! Local variables
    real(r64)              :: PartLoadRatio                     ! Part load ratio                                   [1]
    integer,   parameter   :: DXCoilNum           = 1           ! The number of the DX coil to be simulated         [1]
    integer,   parameter   :: CompOp              = 1           ! Compressor operation; 1=on, 0=off                 [1]
    logical,   parameter   :: FirstHVACIteration  = .TRUE.      ! True if this is the first iteration of HVAC       [1]
    integer,   parameter   :: FanOpMode           = 2           ! Allows parent obj to control fan operation        [1]

    ! Allocate other arrays used within CalcDoe2DXCoil
    ! from DataHeatBalance
    allocate(HeatReclaimDXCoil(1))
    allocate(Schedule(-1:1))
    Schedule(1)%CurrentValue = 1.0

    ! Initialize EnergyPlus global variables
    OutDryBulbTemp  = OutdoorTDryBulb
    OutHumRat       = OutdoorHumRatio
    OutBaroPress    = OutdoorPressure
    OutWetBulbTemp  = PsyTwbFnTdbWPb(OutDryBulbTemp, OutHumRat, OutBaroPress)


    ! Set up the DX coil and initialize data
    allocate(DXCoil(1))
    allocate(DXCoilFullLoadOutAirTemp(1))
    allocate(DXCoilFullLoadOutAirHumRat(1))
    allocate(DXCoilOutletTemp(1))
    allocate(DXCoilOutletHumRat(1))
    allocate(DXCoilPartLoadRatio(1))
    allocate(DXCoilFanOpMode(1))

    DXCoil(1)%SchedPtr                      = 1
    DXCoil(1)%RatedTotCap(1)                = RatedTotCap
    DXCoil(1)%RatedSHR(1)                   = 0.8d0                 ! FIXME: find correct value for RatedSHR
    DXCoil(1)%RatedCBF(1)                   = 0.1d0                 ! FIXME: find correct value for RatedCBF
    DXCoil(1)%RatedEIR(1)                   = 1.0d0 / RatedCOP
    DXCoil(1)%BypassedFlowFrac              = 0.0d0
    DXCoil(1)%RatedAirMassFlowRate          = RatedAirMassFlowRate
    DXCoil(1)%InletAirMassFlowRate          = RatedAirMassFlowRate
    DXCoil(1)%InletAirTemp                  = InletTDryBulb
    DXCoil(1)%InletAirHumRat                = InletHumRatio

    DXCoil(1)%InletAirEnthalpy = PsyHFnTdbW(InletTDryBulb, InletHumRatio)

    ! Performance Curves
    NumCurves = 5
    allocate(PerfCurve(NumCurves))
    ! Capacity as a function of temperature
    DXCoil(1)%CCapFTemp(1)                  = 1
    DXCoil(1)%TotCapTempModFacCurveType(1)  = 4                 ! BiQuadratic curve
    PerfCurve(1)%InterpolationType          = 3                 ! EvaluateCurveToLimits
    PerfCurve(1)%CurveType                  = 4                 ! BiQuadratic
    PerfCurve(1)%Var1Min                    = -10.0d0           ! min inlet air WB temp limit
    PerfCurve(1)%Var2Min                    = -10.0d0           ! min outdoor air DB temp limit
    PerfCurve(1)%Var1Max                    = 50.0d0            ! max inlet air WB temp limit
    PerfCurve(1)%Var2Max                    = 50.0d0            ! max outdoor air DB temp limit
    PerfCurve(1)%CurveMin                   = 0.63d0             ! min value of curve output
    PerfCurve(1)%CurveMax                   = 1.57d0             ! max value of curve output
    PerfCurve(1)%Coeff1                     = 1.5509d0          !- Coefficient1 Constant
    PerfCurve(1)%Coeff2                     = -0.07505d0        !- Coefficient2 x
    PerfCurve(1)%Coeff3                     = 0.0031d0          !- Coefficient3 x2
    PerfCurve(1)%Coeff4                     = 0.0024d0          !- Coefficient4 y
    PerfCurve(1)%Coeff5                     = -0.00005d0        !- Coefficient5 y2
    PerfCurve(1)%Coeff6                     = -0.00043d0        !- Coefficient6 x*y
    PerfCurve(1)%CurveMaxPresent = .TRUE.
    PerfCurve(1)%CurveMinPresent = .TRUE.

    ! Capacity as a function of flow fraction
    DXCoil(1)%CCapFFlow(1)                  = 2
    PerfCurve(2)%InterpolationType          = 3                 ! EvaluateCurveToLimits
    PerfCurve(2)%CurveType                  = 3                 ! Quadratic
    PerfCurve(2)%Var1Min                    = 0.0d0             ! min flow fraction limit
    PerfCurve(2)%Var1Max                    = 1.5d0             ! max flow fraction limit
    PerfCurve(2)%CurveMin                   = 0.0d0             ! min value of curve output
    PerfCurve(2)%CurveMax                   = 2.0d0             ! max value of curve output
    PerfCurve(2)%Coeff1                     = 0.718605468d0     !- Coefficient1 Constant
    PerfCurve(2)%Coeff2                     = 0.410099989d0     !- Coefficient2 x
    PerfCurve(2)%Coeff3                     = -0.128705457d0    !- Coefficient3 x2

    ! EIR as a function of temperature
    DXCoil(1)%EIRFTemp(1)                   = 3
    DXCoil(1)%EIRTempModFacCurveType(1)     = 4                 ! BiQuadratic curve
    PerfCurve(3)%InterpolationType          = 3                 ! EvaluateCurveToLimits
    PerfCurve(3)%CurveType                  = 4                 ! BiQuadratic
    PerfCurve(3)%Var1Min                    = -10.0d0           ! min inlet air WB temp limit
    PerfCurve(3)%Var2Min                    = -10.0d0           ! min outdoor air DB temp limit
    PerfCurve(3)%Var1Max                    = 50.0d0            ! max inlet air WB temp limit
    PerfCurve(3)%Var2Max                    = 50.0d0            ! max outdoor air DB temp limit
    PerfCurve(3)%CurveMin                   = 0.83d0             ! min value of curve output
    PerfCurve(3)%CurveMax                   = 1.21d0             ! max value of curve output
    PerfCurve(3)%Coeff1                     = -0.30428d0        !- Coefficient1 Constant
    PerfCurve(3)%Coeff2                     = 0.11805d0         !- Coefficient2 x
    PerfCurve(3)%Coeff3                     = -0.00342d0        !- Coefficient3 x2
    PerfCurve(3)%Coeff4                     = -0.00626d0        !- Coefficient4 y
    PerfCurve(3)%Coeff5                     = 0.0007d0          !- Coefficient5 y2
    PerfCurve(3)%Coeff6                     = -0.00047d0        !- Coefficient6 x*y
    PerfCurve(3)%CurveMaxPresent = .TRUE.
    PerfCurve(3)%CurveMinPresent = .TRUE.

    ! EIR as a function of flow fraction
    DXCoil(1)%EIRFFlow(1)                   = 4
    PerfCurve(4)%InterpolationType          = 3                 ! EvaluateCurveToLimits
    PerfCurve(4)%CurveType                  = 3                 ! Quadratic
    PerfCurve(4)%Var1Min                    = 0.0d0             ! min flow fraction limit
    PerfCurve(4)%Var1Max                    = 1.5d0             ! max flow fraction limit
    PerfCurve(4)%CurveMin                   = 0.0d0             ! min value of curve output
    PerfCurve(4)%CurveMax                   = 2.0d0             ! max value of curve output
    PerfCurve(4)%Coeff1                     = 1.32299905d0      !- Coefficient1 Constant
    PerfCurve(4)%Coeff2                     = -0.477711207d0    !- Coefficient2 x
    PerfCurve(4)%Coeff3                     = 0.154712157d0     !- Coefficient3 x2

    ! Part load fraction as a function of part load ratio
    DXCoil(1)%PLFFPLR(1)                    = 5
    PerfCurve(5)%InterpolationType          = 3                 ! EvaluateCurveToLimits
    PerfCurve(5)%CurveType                  = 1                 ! Linear
    PerfCurve(5)%Var1Min                    = 0.0d0             ! min inlet air WB temp limit
    PerfCurve(5)%Var1Max                    = 1.0d0             ! max inlet air WB temp limit
    PerfCurve(5)%Coeff1                     = 0.90d0            !- Coefficient1 Constant
    PerfCurve(5)%Coeff2                     = 0.10d0            !- Coefficient2 x
    PerfCurve(5)%CurveMin                   = 0.0d0             ! min value of curve output
    PerfCurve(5)%CurveMax                   = 1.0d0             ! max value of curve output

    ! Run initial case with PartLoadRatio = 1 to get SensibleCoolingLoad and
    ! SensCoolingEnergyRates to calculate actual PartLoadRatio
    PartLoadRatio = 1.0d0
    call CalcDoe2DXCoil(DXCoilNum, CompOp, FirstHVACIteration, PartLoadRatio, FanOpMode)

    ! Calculate PartLoadRatio from previous call
    PartLoadRatio = SensibleCoolingLoad / DXCoil(1)%SensCoolingEnergyRate

    ! Keep PartLoadRatio  bounded between 0 and 1
    if ( PartLoadRatio < 0.0d0 ) PartLoadRatio = 0.0d0
    if ( PartLoadRatio > 1.0d0 ) PartLoadRatio = 1.0d0

    ! Call CalcDoe2DXCoil again with correct part load ratio
    call CalcDoe2DXCoil(DXCoilNum, CompOp, FirstHVACIteration, PartLoadRatio ,FanOpMode)

    ! Assign return variables of interest
    COP = DXCoil(1)%TotalCoolingEnergyRate / DXCoil(1)%ElecCoolingPower
    TotalCoolingCapacity = DXCoil(1)%TotalCoolingEnergyRate / DXCoil(1)%PartLoadRatio
    OutletTemperature = DXCoil(1)%OutletAirTemp
    OutletHumRatio = DXCoil(1)%OutletAirHumRat
    ElecCoolingPower = DXCoil(1)%ElecCoolingPower
    LatCoolingEnergyRate = DXCoil(1)%LatCoolingEnergyRate
    TotalCoolingEnergyRate = DXCoil(1)%TotalCoolingEnergyRate
    TotalSensibleHeatOut = DXCoil(1)%TotalCoolingEnergyRate + DXCoil(1)%ElecCoolingPower

    deallocate(HeatReclaimDXCoil)
    deallocate(Schedule)
    deallocate(DXCoil)
    deallocate(DXCoilFullLoadOutAirTemp)
    deallocate(DXCoilFullLoadOutAirHumRat)
    deallocate(DXCoilOutletTemp)
    deallocate(DXCoilOutletHumRat)
    deallocate(DXCoilPartLoadRatio)
    deallocate(DXCoilFanOpMode)
    deallocate(PerfCurve)

  end subroutine SimCalcDoe2DXCoil
end module EPlusWrapperCooling