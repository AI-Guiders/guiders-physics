using System;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using Guiders.Physics.Core;

namespace Guiders.Physics.Engine;

public readonly record struct KitchenBenchStageSnapshot9(
    Ring9 Minute,
    Ring9 StageIndex,
    Ring9 TargetDishPhase,
    Ring9 TargetShellZ,
    Ring9 CoilLagTicks,
    Ring9 GateLagTicks,
    Ring9 DishWallLagTicks,
    Ring9 TotalPrecompLagTicks,
    Ring9 FpgaProgrammedPhase,
    Ring9 ArrivedDishPhase,
    Ring9 RecoveredTelemetryPhase,
    byte[] Uart0x08CalibrationFrame16B);

/// <summary>
/// Full Digital Twin &amp; Physical World Hardware Adapter of the 12-Beam Cuboctahedron O_2(3) Kitchen Bench (`Tang Nano 9K` 27.000 MHz):
/// - Bridges digital `Ring9 / Bus36 / VerletCell72` states to the physical matter at the focal center of the 12 coils (`0`-tick error),
/// - Models all 12 hand-wound copper coils (`Z(Cu) = 29`, `N = 60` turns, `MMF = 42 = S_2(3)` At, Wheeler denom `202`),
/// - Models MOSFET gate stack (`Z(Si+SiO2) = 44`, Miller plateau lag `25 = d_out(E)^2` ticks),
/// - Models 1N5819 Schottky flyback relaxation (`3 * 81 + 13 = 256` ticks),
/// - Models quartz Petri dish wall (`Z(SiO2) = 30 == Z(N2+O2)`, Fresnel 9:8 transmission `288/289`),
/// - Generates and verifies 16-byte XOR-CRC UART frames (`0x06`, `0x08`, `0x02`),
/// - Executes the 52-minute (`52 = 4 * 13 = Z(Te)`) 6-stage synthesis schedule (`Stage 0..5`) with 0-tick phase error!
/// </summary>
public sealed unsafe class KitchenBenchDigitalTwin
{
    private static readonly delegate* unmanaged<ushort> _rdseed9Ptr = InitHardwareRdseed9();

    [DllImport("kernel32.dll", SetLastError = true)]
    private static extern IntPtr VirtualAlloc(IntPtr lpAddress, nuint dwSize, uint flAllocationType, uint flProtect);

    private static delegate* unmanaged<ushort> InitHardwareRdseed9()
    {
        IntPtr mem = VirtualAlloc(IntPtr.Zero, 64, 0x3000, 0x40);
        byte* p = (byte*)mem;
        p[0] = 0x48; p[1] = 0x0F; p[2] = 0xC7; p[3] = 0xF8;
        p[4] = 0x73; p[5] = 0xFA;
        p[6] = 0x25; p[7] = 0xFF; p[8] = 0x01; p[9] = 0x00; p[10] = 0x00;
        p[11] = 0xC3;
        return (delegate* unmanaged<ushort>)mem;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 SampleSiliconThermalTrng9() => new(_rdseed9Ptr());

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 ThermalPhaseNow, Ring9 ThermalDeltaPhiVerlet) SampleSiliconThermalTwoTickPair9()
    {
        Ring9 p0 = SampleSiliconThermalTrng9();
        Ring9 p1 = SampleSiliconThermalTrng9();
        return (p1, p1 - p0);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static VerletCell72 CreatePhysicalAlphaQuartzCell(Ring9 crystalOscillatorA, Ring9 crystalOscillatorB)
    {
        Ring9 phiSi = Ring9.QuartzSiO2Phi + crystalOscillatorA;
        Ring9 phiO2 = Ring9.QuartzSiO2Phi + crystalOscillatorB;
        Ring9 dphiBeat = phiO2 - phiSi;
        var cell = VerletCell72.Create(Ring9.ZSilicon14, phiSi, phiSi - dphiBeat, Ring9.ZOxygen2_16, phiO2, phiO2 - dphiBeat);
        cell.Gate0 = CuboctahedronAutomatonKernel.ApplyHeisenbergDeadband(dphiBeat);
        cell.Gate1 = cell.Gate0;
        return cell;
    }

    private readonly Ring9[] _perCoilSkewTicks = new Ring9[12];
    public Ring9 GateLagTicks { get; set; } = Ring9.BezoutInv41; // 25 ticks = 925.9 ns Miller plateau
    public Ring9 DishWallLagTicks { get; set; } = Ring9.Zero;    // 0 ticks for SiO2 quartz glass, 7 for polystyrene

    public KitchenBenchDigitalTwin()
    {
        for (int k = 0; k < 12; k++)
            _perCoilSkewTicks[k] = Ring9.Zero;
    }

    public void SetCoilWindingSkew(int vertexIndex, Ring9 skewTicks) =>
        _perCoilSkewTicks[vertexIndex % 12] = skewTicks;

    public static byte[] PackUart16ByteFrame(byte opcode, uint arg0, uint arg1, uint arg2)
    {
        byte[] frame = new byte[16];
        frame[0] = 0xA5;
        frame[1] = (byte)(opcode & 0x0F);
        frame[2] = (byte)(arg0 & 0xFF);
        frame[3] = (byte)((arg0 >> 8) & 0xFF);
        frame[4] = (byte)((arg0 >> 16) & 0xFF);
        frame[5] = (byte)((arg0 >> 24) & 0xFF);
        frame[6] = (byte)(arg1 & 0xFF);
        frame[7] = (byte)((arg1 >> 8) & 0xFF);
        frame[8] = (byte)((arg1 >> 16) & 0xFF);
        frame[9] = (byte)((arg1 >> 24) & 0xFF);
        frame[10] = (byte)(arg2 & 0xFF);
        frame[11] = (byte)((arg2 >> 8) & 0xFF);
        frame[12] = (byte)((arg2 >> 16) & 0xFF);
        frame[13] = (byte)((arg2 >> 24) & 0xFF);
        byte crc = 0;
        for (int i = 1; i <= 13; i++)
            crc ^= frame[i];
        frame[14] = crc;
        frame[15] = 0x5A;
        return frame;
    }

    private static readonly Ring9[] StageIndices =
        [Ring9.Zero, Ring9.MuVac, Ring9.C0OverCSq, Ring9.D, Ring9.DInSquare, Ring9.DOutEdge];

    private static readonly Ring9[] StageTargetDphi =
        [Ring9.CollapsePhi, Ring9.CollapsePhi, Ring9.CollapsePhi, Ring9.QSquare, Ring9.CollapsePhi, Ring9.Zero];

    private static readonly Ring9[] StageTargetZ =
    [
        Ring9.NTri + Ring9.C0OverCSq, // 10 (H2O condensation on TEC1-12706)
        Ring9.MHat - Ring9.MuVac,     // 96 (Glucose + 6 O2 photosynthesis)
        Ring9.KP + Ring9.NSq,         // 70 (Chiral L-amino & Heme FeN4O2)
        Ring9.KP + Ring9.NSq,         // 70 (12-beam 52.734 kHz ultrasonic folding)
        Ring9.KP + Ring9.NSq,         // 70 (68 C Maillard golden crust)
        Ring9.KP + Ring9.NSq          // 70 (Phase stasis lock D = 6,823,936)
    ];

    public static bool VerifyUart16ByteFrameCrc(byte[] frame)
    {
        bool validHeader = frame is not null && frame.Length == 16 && frame[0] == 0xA5 && frame[15] == 0x5A;
        byte crc = 0;
        int lim = validHeader ? 13 : 0;
        for (int i = 1; i <= lim; i++)
            crc ^= frame![i];
        return validHeader && (crc == frame![14]);
    }

    public KitchenBenchStageSnapshot9 EvaluateMinute(int minute)
    {
        int m = Math.Clamp(minute, 0, 52);
        int s = (((3 - m) >> 31) & 1)
              + (((17 - m) >> 31) & 1)
              + (((33 - m) >> 31) & 1)
              + (((45 - m) >> 31) & 1)
              + (((51 - m) >> 31) & 1);

        Ring9 stageIdx = StageIndices[s];
        Ring9 targetDphi = StageTargetDphi[s];
        Ring9 targetZ = StageTargetZ[s];

        // Thermal coil lag shifts from 64 ticks (25 C) to 58 ticks (63.8 C at 52 min):
        ushort thermalDrop = (ushort)((m * 6) / 52);
        Ring9 coilLag = Ring9.KP - new Ring9(thermalDrop);
        Ring9 totalLag = coilLag + GateLagTicks + DishWallLagTicks;
        Ring9 fpgaDrive = ComputeKitchenBenchPrecompensatedBeamPhase9(
            targetDphi, GateLagTicks, coilLag, Ring9.Zero, DishWallLagTicks, Ring9.Zero);
        Ring9 arrived = ComputeKitchenBenchArrivedDishPhase9(
            fpgaDrive, GateLagTicks, coilLag, Ring9.Zero, DishWallLagTicks, Ring9.Zero);

        var dishCell = VerletCell72.Create(Ring9.MuVac, arrived, Ring9.Zero, Ring9.MuVac, Ring9.Zero, Ring9.Zero);
        var telemetry = CuboctahedronAutomatonKernel.ReadCellSensors(dishCell);
        Ring9 recovered = RecoverRingPhaseFromHydrophoneAndChiralCurrent9(
            Ring9.MuVac, Ring9.MuVac, telemetry.AcousticWave, telemetry.ReactionCurrent);

        uint lagInt = totalLag.Bits9;
        byte[] uart0x08 = PackUart16ByteFrame(0x08, 1u, lagInt, 2u);

        return new KitchenBenchStageSnapshot9(
            new Ring9((ushort)m),
            stageIdx,
            targetDphi,
            targetZ,
            coilLag,
            GateLagTicks,
            DishWallLagTicks,
            totalLag,
            fpgaDrive,
            arrived,
            recovered,
            uart0x08);
    }

    public static (Ring9 UpdatedEstimatedLag, Ring9 NextTickArrivedDishPhase, Ring9 NextTickPhaseError)
        StepClosedLoopPll9(Ring9 targetDishPhase, Ring9 currentEstimatedLag, Ring9 trueHardwareLagWithDrift)
    {
        Ring9 fpgaDrive0 = targetDishPhase + currentEstimatedLag;
        Ring9 actualDish0 = fpgaDrive0 - trueHardwareLagWithDrift;
        var cell0 = VerletCell72.Create(Ring9.MuVac, actualDish0, Ring9.Zero, Ring9.MuVac, Ring9.Zero, Ring9.Zero);
        var sensors0 = CuboctahedronAutomatonKernel.ReadCellSensors(cell0);
        Ring9 recDish0 = RecoverRingPhaseFromHydrophoneAndChiralCurrent9(
            Ring9.MuVac, Ring9.MuVac, sensors0.AcousticWave, sensors0.ReactionCurrent);
        Ring9 phaseErr0 = targetDishPhase - recDish0;
        Ring9 nextEstLag = currentEstimatedLag + phaseErr0;
        Ring9 nextFpgaDrive = targetDishPhase + nextEstLag;
        Ring9 nextArrivedDish = nextFpgaDrive - trueHardwareLagWithDrift;
        return (nextEstLag, nextArrivedDish, targetDishPhase - nextArrivedDish);
    }

    public static Ring9 GetRecipeTargetShellZ9(byte recipeIndex) =>
        recipeIndex switch
        {
            0 => Ring9.DOutEdge * Ring9.NLinks,                  // 60 = Gate #0 Bromate BrO3- (5 * 12)
            1 => Ring9.KP + Ring9.NSq,                           // 70 = Cutlet Heme FeN4O2 (64 + 6)
            2 => new Ring9(81) + Ring9.DOutEdge,                 // 86 = Bread Amylose Z(Rn) (81 + 5)
            3 => new Ring9(92),                                  // 92 = Vitamin C S_3(3)
            4 => new Ring9(102),                                 // 102 = Cappuccino Caffeine / BaTiO3 (6 * 17)
            _ => Ring9.NLinks * Ring9.NLinks                     // 144 = Palmitic Lipid / Cheese (12 * 12)
        };

    public static (Ring9 ReleasePulsePhase, Ring9 PauliEjectPhase, byte[] Uart0x02StartFrame16B)
        BuildRecipeAndDryCleanProtocol9(byte recipeIndex)
    {
        byte[] startFrame = PackUart16ByteFrame(0x02, recipeIndex, 1u, 0u);
        return (Ring9.CollapsePhi, Ring9.PhiPauli, startFrame);
    }

    public static bool VerifyAcousticFixedPointAcross4OperatingModes9()
    {
        Ring9[] modes = [Ring9.Zero, Ring9.QSquare, Ring9.CollapsePhi, Ring9.PhiPauli];
        bool allMatch = true;
        for (int i = 0; i < modes.Length; i++)
        {
            Ring9 m = modes[i];
            var cell = VerletCell72.Create(Ring9.MuVac, m, Ring9.Zero, Ring9.MuVac, Ring9.Zero, Ring9.Zero);
            var s = CuboctahedronAutomatonKernel.ReadCellSensors(cell);
            allMatch &= (s.AcousticWave == m);
        }
        return allMatch;
    }

    public static bool Run3DPetriLatticeReversibleCycle9(ushort latticeSize = 4)
    {
        Ring9 size9 = new(latticeSize);
        var carrier = new CuboctahedronCarrier(size9);
        Ring9 c = new((ushort)(latticeSize / 2));
        VerletCell72 initCenter = carrier[c, c, c];
        carrier.StepForward();
        carrier.StepForward();
        carrier.StepBackward();
        carrier.StepBackward();
        VerletCell72 finalCenter = carrier[c, c, c];
        return carrier.Tick == Ring9.Zero
            && finalCenter.Mu0 == initCenter.Mu0
            && finalCenter.Phi0Curr == initCenter.Phi0Curr
            && finalCenter.Phi0Past == initCenter.Phi0Past
            && finalCenter.Mu1 == initCenter.Mu1
            && finalCenter.Phi1Curr == initCenter.Phi1Curr
            && finalCenter.Phi1Past == initCenter.Phi1Past;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 CascadeAttenuatedStage9(Ring9 phiIn, Ring9 impedanceRatio)
    {
        ushort r = (ushort)((1 - (1 / (impedanceRatio.Bits9 + 1))) * impedanceRatio.Bits9 + (1 / (impedanceRatio.Bits9 + 1)));
        ushort mag = phiIn.AbsWrap.Bits9;
        ushort q = (ushort)(mag / r);
        ushort isNeg = phiIn.BitQ1.Bits9;
        Ring9 divided = new((ushort)((1 - 2 * isNeg) * q));
        return CuboctahedronAutomatonKernel.ApplyHeisenbergDeadband(divided);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 CascadeUncompensatedArrival9(Ring9 phiDrive, Ring9 j1, Ring9 j2, Ring9 j3, Ring9 j4) =>
        phiDrive + j1 + j2 + j3 + j4;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 CascadePrecompensatedDrive9(Ring9 phiTarget, Ring9 j1, Ring9 j2, Ring9 j3, Ring9 j4) =>
        phiTarget - (j1 + j2 + j3 + j4);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 CascadeSubthresholdRotatorRelay4Stage9(Ring9 dphiSubthreshold) =>
        dphiSubthreshold.ProjS.ProjS.ProjS.ProjS;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 Y0, Ring9 Y1, Ring9 Y2) ApplyNilpotentPerturbation3(Ring9 s0, Ring9 s1, Ring9 s2) =>
        (s0 + s1, s1 + s2, s2);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 S0, Ring9 S1, Ring9 S2) CalibrateFiniteNeumannOrder2(Ring9 y0, Ring9 y1, Ring9 y2) =>
        (y0 - y1 + y2, y1 - y2, y2);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 ComputeKitchenBenchPrecompensatedBeamPhase9(
        Ring9 targetDishPhase,
        Ring9 mosfetGateLag,
        Ring9 copperCoilLag,
        Ring9 thermalDriftLag,
        Ring9 dishWallLag,
        Ring9 perCoilWindingSkew) =>
        targetDishPhase + mosfetGateLag + copperCoilLag + thermalDriftLag + dishWallLag + perCoilWindingSkew;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 ComputeKitchenBenchArrivedDishPhase9(
        Ring9 fpgaBeamDrivePhase,
        Ring9 mosfetGateLag,
        Ring9 copperCoilLag,
        Ring9 thermalDriftLag,
        Ring9 dishWallLag,
        Ring9 perCoilWindingSkew) =>
        fpgaBeamDrivePhase - (mosfetGateLag + copperCoilLag + thermalDriftLag + dishWallLag + perCoilWindingSkew);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 ComputeWheelerAirCoreDenominator9(Ring9 radiusMm, Ring9 lengthMm, Ring9 depthMm) =>
        (Ring9.NSq * radiusMm) + (Ring9.ModBits * lengthMm) + ((Ring9.NTri + Ring9.C0OverCSq) * depthMm);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 ComputeColdStartFpgaPrecompPhase9() =>
        Ring9.CollapsePhi + Ring9.KP + Ring9.BezoutInv41;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 ComputeFresnelQuartzWaterSum9() =>
        Ring9.ModBits + Ring9.NTri;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 ComputeStasisMultiplierMod512() =>
        Ring9.MuVac + (Ring9.NPhi * Ring9.KP);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 ComputeSucroseDisaccharideShell9() =>
        (Ring9.C0OverCSq * (Ring9.MHat - Ring9.MuVac)) - (Ring9.NTri + Ring9.C0OverCSq);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 RecoverRingPhaseFromHydrophoneAndChiralCurrent9(
        Ring9 mu0, Ring9 mu1, Ring9 acousticSoundWave, Ring9 chiralReactionCurrent)
    {
        ushort prod = (ushort)(mu0.Bits9 * mu1.Bits9);
        ushort safeProd = (ushort)((1 - (1 / (prod + 1))) * prod + (1 / (prod + 1)));
        ushort oneMinusC = (ushort)(acousticSoundWave.Bits9 / safeProd);
        ushort cNeg = (ushort)(1 - (1 / ((oneMinusC / 129) + 1)));
        ushort cPos = (ushort)(1 / ((oneMinusC / 128) + 1));
        ushort cZero = (ushort)(1 - cPos - cNeg);
        ushort sNeg = chiralReactionCurrent.BitQ1.Bits9;
        ushort sMag = (ushort)(chiralReactionCurrent.AbsWrap.Bits9 / safeProd);

        ushort q00 = (ushort)((1 - sNeg) * (cPos + cZero));
        ushort q01 = (ushort)((1 - sNeg) * cNeg);
        ushort q10 = (ushort)(sNeg * (cNeg + cZero));
        ushort q11 = (ushort)(sNeg * cPos);

        ushort phiRec = (ushort)(
            q00 * sMag +
            q01 * (256 - sMag) +
            q10 * (256 + sMag) +
            q11 * (512 - sMag));
        return new Ring9(phiRec);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Bus36 NilpotentTomographForward2xBus36(Bus36 cleanSensors4) =>
        new(
            (Ring9.C0OverCSq * cleanSensors4.L0) + cleanSensors4.L1,
            (Ring9.C0OverCSq * cleanSensors4.L1) + (Ring9.C0OverCSq * cleanSensors4.L2),
            Ring9.C0OverCSq * cleanSensors4.L2,
            Ring9.C0OverCSq * cleanSensors4.L3);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Bus36 NilpotentTomographInverse24xBus36(Bus36 mixedSensors4)
    {
        ushort y0 = mixedSensors4.L0.Bits9;
        ushort y1 = mixedSensors4.L1.Bits9;
        ushort y2 = mixedSensors4.L2.Bits9;
        ushort y3 = mixedSensors4.L3.Bits9;

        ushort s0 = (ushort)(((12 * y0) + (6 * y2) - (6 * y1)) / 24);
        ushort s1 = (ushort)(((12 * y1) - (12 * y2)) / 24);
        ushort s2 = (ushort)((12 * y2) / 24);
        ushort s3 = (ushort)((12 * y3) / 24);

        return new Bus36(new Ring9(s0), new Ring9(s1), new Ring9(s2), new Ring9(s3));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Bus36 ComputePreflightCoilImpedanceBus36(Ring9 activeCoilsCount) =>
        new(
            Ring9.Zero,
            Ring9.C0OverCSq * Ring9.NTri,
            Ring9.C0OverCSq * (activeCoilsCount + Ring9.MuVac),
            Ring9.Zero);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 RecoverActiveCoilsFromImpedanceBus36(Bus36 domeImpedance) =>
        new((ushort)((domeImpedance.L2.Bits9 / 2) - 1));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 ComputeGorkovTrapHessianTrace9(Ring9 activeCoilsCount) =>
        Ring9.C0OverCSq * activeCoilsCount;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 ConcentricRingsCount6, Ring9 WallClearanceMm3) ComputePetriDishBjerknesRings9() =>
        (new Ring9((ushort)((Ring9.DOutEdge.Bits9 * Ring9.ModBits.Bits9) / Ring9.DOutVertex.Bits9)),
         new Ring9((ushort)((Ring9.DOutEdge.Bits9 * Ring9.ModBits.Bits9) % Ring9.DOutVertex.Bits9)));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 SingleRayPhase21, Ring9 SingleRayWallPhotons0, Ring9 Coherent12BeamFocalCurrent252)
        Compute12BeamSubcavitationFocus9()
    {
        Ring9 singleRay = new((ushort)(Ring9.PhiPauli.Bits9 / Ring9.NLinks.Bits9));
        Ring9 wallPhotons = CuboctahedronAutomatonKernel.HeisenbergEnergyQuanta(singleRay);
        Ring9 focalCurrent = Ring9.NLinks * singleRay.ProjS;
        return (singleRay, wallPhotons, focalCurrent);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 ComputeEulerQuadrupoleChiralTorque9(Ring9 poleStepPhase) =>
        poleStepPhase.ProjS;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 DewPoint12C, Ring9 WarmKitchen27C, Ring9 Maillard68C, Ring9 Cutoff81C, Ring9 Energy156Wh)
        ComputeKitchenBenchThermalLadderAndEnergy9()
    {
        Ring9 dew12 = Ring9.NLinks;
        Ring9 warm27 = Ring9.D * Ring9.D * Ring9.D;
        Ring9 maillard68 = Ring9.KP + Ring9.DInSquare;
        Ring9 cutoff81 = warm27 * Ring9.D;
        Ring9 power180 = Ring9.NLinks * (Ring9.D * Ring9.DOutEdge);
        Ring9 min52 = Ring9.DInSquare * Ring9.NPhi;
        Ring9 min60 = Ring9.DInSquare * (Ring9.D * Ring9.DOutEdge);
        Ring9 energy156 = new((ushort)((power180.Bits9 * min52.Bits9) / min60.Bits9));
        return (dew12, warm27, maillard68, cutoff81, energy156);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 PyrazineAroma42, Ring9 PeptideBond22, Ring9 HydratedPeptide32, Bus36 AudibleSubharmonicClocks32768)
        ComputeMaillardAromaAndAudibleSubharmonic9()
    {
        Ring9 pyrazine42 = (Ring9.DInSquare * Ring9.NSq) + (Ring9.DInSquare * Ring9.MuVac) + (Ring9.C0OverCSq * Ring9.DOutVertex);
        Ring9 peptide22 = Ring9.NSq + Ring9.NTri + Ring9.DOutVertex + Ring9.MuVac;
        Ring9 hydrated32 = peptide22 + (Ring9.NTri + Ring9.C0OverCSq);
        Bus36 audible32768 = new(Ring9.Zero, Ring9.KP, Ring9.Zero, Ring9.Zero); // 0 + 64 * 512 = 32,768 = 4 * eps
        return (pyrazine42, peptide22, hydrated32, audible32768);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 WaterDielectric20C_80, Ring9 WaterDielectric68C_64, Ring9 UmamiGlutamicAcid78, Ring9 Cysteine64, Ring9 CystineDisulfide126, Ring9 Section36Count)
        ComputeBench36SectionFullClosure9()
    {
        Ring9 eps20C = Ring9.C0OverCSq * (Ring9.DPhiDisc - Ring9.MuVac); // 2 * 40 = 80
        Ring9 eps68C = eps20C - Ring9.ZOxygen2_16;                       // 80 - 16 = 64 = K_P
        Ring9 umami78 = Ring9.NSq * Ring9.NPhi;                          // 6 * 13 = 78 = Z(Pt)
        Ring9 cysteine64 = Ring9.KP;                                     // 64 = K_P
        Ring9 cystine126 = (Ring9.C0OverCSq * cysteine64) - Ring9.C0OverCSq; // 128 - 2 = 126 = 2 * 63
        Ring9 sec36 = Ring9.DInSquare * Ring9.ModBits;                   // 4 * 9 = 36 = dim(Bus36)
        return (eps20C, eps68C, umami78, cysteine64, cystine126, sec36);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 TheobromineChocolate94, Ring9 CitricAcidNectar100, Ring9 LacticAcidMilk48, Ring9 CalciumLactateCheese114, Ring9 WorstCaseNoise8, Ring9 HeisenbergMargin33, Ring9 FullSection42Count)
        ComputeBench42SecondSphereClosure9()
    {
        Ring9 glucose96 = Ring9.MHat - Ring9.MuVac;
        Ring9 theobromine94 = glucose96 - Ring9.C0OverCSq;
        Ring9 citric100 = Ring9.DInSquare * Ring9.BezoutInv41;
        Ring9 lactic48 = Ring9.DInSquare * Ring9.NLinks;
        Ring9 calciumLactate114 = Ring9.NSq * (Ring9.MuVac + Ring9.NLinks + Ring9.NSq);
        Ring9 worstNoise8 = Ring9.MuVac + Ring9.MuVac + Ring9.D + Ring9.C0OverCSq + Ring9.MuVac;
        Ring9 margin33 = Ring9.DPhiDisc - worstNoise8;
        Ring9 sec42 = (Ring9.DInSquare * Ring9.ModBits) + Ring9.NSq;
        return (theobromine94, citric100, lactic48, calciumLactate114, worstNoise8, margin33, sec42);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 MainsPeak312V, Ring9 CanonicalAminoAcids20, Ring9 AspartateHeme70, Ring9 HistidineLigand82, Ring9 TryptophanBus108, Ring9 MujocoGeoms49, Ring9 TwoShellBall55)
        ComputeBench55TwoShellBallClosure9()
    {
        Ring9 dc24 = Ring9.C0OverCSq * Ring9.NLinks;
        Ring9 mains312 = dc24 * Ring9.NPhi;
        Ring9 aa20 = Ring9.DInSquare * Ring9.DOutEdge;
        Ring9 asp70 = Ring9.KP + Ring9.NSq;
        Ring9 his82 = new Ring9(81) + Ring9.MuVac;
        Ring9 trp108 = Ring9.NLinks * Ring9.ModBits;
        Ring9 geoms49 = Ring9.DOutVertex * Ring9.DOutVertex;
        Ring9 ball55 = Ring9.DOutEdge * (Ring9.NLinks - Ring9.MuVac);
        return (mains312, aa20, asp70, his82, trp108, geoms49, ball55);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 HotCoilLag58, Ring9 MidCoilLag59, Ring9 EarlyCoilLag61, Ring9 MinCoilLag62, Ring9 HydrogenBond63, Ring9 StasisTemp65, Ring9 MaxCoilLag66, Ring9 Verlet72SectionCount)
        ComputeBench72VerletCellClosure9()
    {
        Ring9 lag58 = Ring9.KP - Ring9.NSq;
        Ring9 lag59 = Ring9.KP - Ring9.DOutEdge;
        Ring9 lag61 = Ring9.KP - Ring9.D;
        Ring9 lag62 = Ring9.KP - Ring9.C0OverCSq;
        Ring9 hb63 = Ring9.KP - Ring9.MuVac;
        Ring9 temp65 = Ring9.DOutEdge * Ring9.NPhi;
        Ring9 lag66 = Ring9.KP + Ring9.C0OverCSq;
        Ring9 sec72 = Ring9.NTri * Ring9.ModBits;
        return (lag58, lag59, lag61, lag62, hb63, temp65, lag66, sec72);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 FrameRadius73Mm, Ring9 Optocoupler74Ns, Ring9 CollagenGate75, Ring9 TetraOctahedral76, Ring9 GlutamateAnion77, Ring9 GoldPlating79, Ring9 Heisenberg81SectionCount)
        ComputeBench81HeisenbergClassClosure9()
    {
        Ring9 r73 = Ring9.KP + Ring9.ModBits;
        Ring9 ns74 = Ring9.C0OverCSq * ((Ring9.DOutVertex * Ring9.DOutEdge) + Ring9.C0OverCSq);
        Ring9 col75 = Ring9.D * Ring9.BezoutInv41;
        Ring9 tet76 = Ring9.DInSquare * (Ring9.MuVac + Ring9.NLinks + Ring9.NSq);
        Ring9 glu77 = (Ring9.NSq * Ring9.NPhi) - Ring9.MuVac;
        Ring9 gold79 = (Ring9.NSq * Ring9.NPhi) + Ring9.MuVac;
        Ring9 sec81 = (Ring9.NTri * Ring9.ModBits) + Ring9.ModBits;
        return (r73, ns74, col75, tet76, glu77, gold79, sec81);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 HotTotalLag83, Ring9 MidTotalLag84, Ring9 Lag20Min85, Ring9 EarlyTotalLag86, Ring9 SrfTripleCu87, Ring9 MsgLag88, Ring9 ColdNominalLag89, Ring9 MaxCoil11Lag91, Ring9 PolystyreneTotalLag96, Ring9 GaussianMhat97SectionCount)
        ComputeBench97GaussianModulusClosure9()
    {
        Ring9 gate25 = Ring9.BezoutInv41;
        Ring9 lag83 = (Ring9.KP - Ring9.NSq) + gate25;
        Ring9 lag84 = (Ring9.KP - Ring9.DOutEdge) + gate25;
        Ring9 lag85 = (Ring9.DOutEdge * Ring9.NLinks) + gate25;
        Ring9 lag86 = (Ring9.KP - Ring9.D) + gate25;
        Ring9 srf87 = (Ring9.KP - Ring9.C0OverCSq) + gate25;
        Ring9 msg88 = (Ring9.KP - Ring9.MuVac) + gate25;
        Ring9 cold89 = Ring9.KP + gate25;
        Ring9 max91 = (Ring9.KP + Ring9.C0OverCSq) + gate25;
        Ring9 poly96 = cold89 + Ring9.DOutVertex;
        Ring9 mhat97 = (Ring9.ModBits * Ring9.ModBits) + (Ring9.DInSquare * Ring9.DInSquare);
        return (lag83, lag84, lag85, lag86, srf87, msg88, cold89, max91, poly96, mhat97);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 WheelerPrime101, Ring9 KokOxoIodine106, Ring9 CollapseOddCurrent117, Ring9 MeniscusFresnel119, Ring9 KokS3Flash123, Ring9 PeltierCouples127, Ring9 QuarterRing128SectionCount)
        ComputeBench128QuarterRingSaddleClosure9()
    {
        Ring9 p101 = (Ring9.DInSquare * Ring9.BezoutInv41) + Ring9.MuVac;
        Ring9 kok106 = Ring9.C0OverCSq * new Ring9(53);
        Ring9 cur117 = Ring9.ModBits * Ring9.NPhi;
        Ring9 mf119 = Ring9.DOutVertex * (Ring9.ModBits + Ring9.NTri);
        Ring9 s3_123 = Ring9.D * Ring9.DPhiDisc;
        Ring9 pelt127 = Ring9.QSquare - Ring9.MuVac;
        Ring9 sec128 = Ring9.QSquare;
        return (p101, kok106, cur117, mf119, s3_123, pelt127, sec128);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 MeniscusGate129, Ring9 OrderedPairs132, Ring9 Cesium133Mass, Ring9 FineStructure137, Ring9 CollapsePhase139, Ring9 PalmiticLipid144, Ring9 ThreeShellBall147SectionCount)
        ComputeBench147ThreeShellBallClosure9()
    {
        Ring9 gate129 = Ring9.QSquare + Ring9.MuVac;
        Ring9 pairs132 = Ring9.NLinks * (Ring9.NLinks - Ring9.MuVac);
        Ring9 cs133 = Ring9.DOutVertex * (Ring9.MuVac + Ring9.NLinks + Ring9.NSq);
        Ring9 fs137 = Ring9.QSquare + Ring9.ModBits;
        Ring9 coh139 = Ring9.CollapsePhi;
        Ring9 lip144 = Ring9.NLinks * Ring9.NLinks;
        Ring9 sec147 = Ring9.D * (Ring9.DOutVertex * Ring9.DOutVertex);
        return (gate129, pairs132, cs133, fs137, coh139, lip144, sec147);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 KokS4Flash164, Ring9 CoreSquare169, Ring9 Power180SectionCount, Ring9 IvSeedBezout183, Ring9 Gold197SectionCount, Ring9 FpgaCold228SectionCount, Ring9 Uart234SectionCount, Ring9 Pauli256SectionCount, Ring9 IvSeed263SectionCount, Ring9 IcePoint273K, Ring9 Fresnel289SectionCount, Ring9 StandardTemp298K, Ring9 TwelveGateSum300, Ring9 Mains312SectionCount, Ring9 Maillard341K, Ring9 Branching343SectionCount, Ring9 TwelveCopper348, Ring9 Phase352SectionCount, Ring9 Cutoff354K, Ring9 FullCircle360, Ring9 Boiling373K, Ring9 Euler384SectionCount, Ring9 WireDiameter400Um, Ring9 GaussianBezout417, Ring9 Soret418Nm, Ring9 TwelveBus432, Ring9 BlueLed441SectionCount)
        ComputeBench441BlueLedMultiscaleClosure9()
    {
        Ring9 kok164 = Ring9.DInSquare * Ring9.DPhiDisc;
        Ring9 sq169 = Ring9.NPhi * Ring9.NPhi;
        Ring9 sec180 = Ring9.NLinks * (Ring9.D * Ring9.DOutEdge);
        Ring9 inv183 = Ring9.D * (Ring9.KP - Ring9.D);
        Ring9 au197 = (Ring9.C0OverCSq * Ring9.MHat) + Ring9.D;
        Ring9 fpga228 = Ring9.CollapsePhi + Ring9.KP + Ring9.BezoutInv41;
        Ring9 uart234 = (Ring9.C0OverCSq * Ring9.ModBits) * Ring9.NPhi;
        Ring9 pauli256 = Ring9.PhiPauli;
        Ring9 iv263 = Ring9.IvSeed;
        Ring9 ice273 = Ring9.SeamTicks * Ring9.NPhi;
        Ring9 fresnel17 = Ring9.ModBits + Ring9.NTri;
        Ring9 sec289 = fresnel17 * fresnel17;
        Ring9 s2_42 = (Ring9.DInSquare * Ring9.ModBits) + Ring9.NSq;
        Ring9 temp298 = pauli256 + s2_42;
        Ring9 gate300 = Ring9.NLinks * Ring9.BezoutInv41;
        Ring9 mains312 = (Ring9.C0OverCSq * Ring9.NLinks) * Ring9.NPhi;
        Ring9 maillard341 = (Ring9.KP + Ring9.DInSquare) + ice273;
        Ring9 branch343 = Ring9.DOutVertex * (Ring9.DOutVertex * Ring9.DOutVertex);
        Ring9 cu29 = (Ring9.DInSquare * Ring9.DOutVertex) + Ring9.MuVac;
        Ring9 cu348 = Ring9.NLinks * cu29;
        Ring9 ph352 = pauli256 + (Ring9.MHat - Ring9.MuVac);
        Ring9 cutoff354 = new Ring9(81) + ice273;
        Ring9 circle360 = Ring9.C0OverCSq * sec180;
        Ring9 boil373 = pauli256 + (Ring9.ModBits * Ring9.NPhi);
        Ring9 euler384 = Ring9.D * Ring9.QSquare;
        Ring9 aa20 = Ring9.DInSquare * Ring9.DOutEdge;
        Ring9 wire400 = aa20 * aa20;
        Ring9 bezout417 = Ring9.D * Ring9.CollapsePhi;
        Ring9 soret418 = bezout417 + Ring9.MuVac;
        Ring9 bus432 = Ring9.NLinks * (Ring9.DInSquare * Ring9.ModBits);
        Ring9 led441 = Ring9.SeamTicks * Ring9.SeamTicks;
        return (kok164, sq169, sec180, inv183, au197, fpga228, uart234, pauli256, iv263, ice273, sec289, temp298, gate300, mains312, maillard341, branch343, cu348, ph352, cutoff354, circle360, boil373, euler384, wire400, bezout417, soret418, bus432, led441);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 TwelveTickScan444Ns, Ring9 SevenMeniscusKp448, Ring9 DoubleFpgaDrive456, Ring9 MooreFresnelSrf459, Ring9 Co2SeamQuartz462, Ring9 Bus36CoreUart468, Ring9 UpperHeisenberg471, Ring9 ThirteenTickCore481Ns, Ring9 Co2Square484, Ring9 RingMinusGate487, Ring9 RingMinusSeam491, Ring9 TwelveHeisenberg492, Ring9 ThirdPerfect496, Ring9 SchottkyEnd499, Ring9 TwistedPair500Mm, Ring9 QuadPeltierFocal508, Ring9 QuartzFresnelRg510, Ring9 RingMax511, Bus36 FullRing512SectionCount)
        ComputeBench512FullRingMultiscaleClosure9()
    {
        Ring9 tick37 = (Ring9.DOutVertex * Ring9.DOutEdge) + Ring9.C0OverCSq;
        Ring9 scan444 = Ring9.NLinks * tick37;
        Ring9 kp448 = Ring9.DOutVertex * Ring9.KP;
        Ring9 fpga228 = Ring9.CollapsePhi + Ring9.KP + Ring9.BezoutInv41;
        Ring9 fpga456 = Ring9.C0OverCSq * fpga228;
        Ring9 fresnel17 = Ring9.ModBits + Ring9.NTri;
        Ring9 srf459 = (Ring9.D * (Ring9.D * Ring9.D)) * fresnel17;
        Ring9 co2_22 = (Ring9.C0OverCSq * Ring9.ModBits) + Ring9.DInSquare;
        Ring9 qz462 = co2_22 * Ring9.SeamTicks;
        Ring9 uart468 = (Ring9.DInSquare * Ring9.ModBits) * Ring9.NPhi;
        Ring9 heis471 = Ring9.Zero - Ring9.DPhiDisc;
        Ring9 core481 = tick37 * Ring9.NPhi;
        Ring9 co2Sq484 = co2_22 * co2_22;
        Ring9 gate487 = Ring9.Zero - Ring9.BezoutInv41;
        Ring9 seam491 = Ring9.Zero - Ring9.SeamTicks;
        Ring9 heis492 = Ring9.NLinks * Ring9.DPhiDisc;
        Ring9 ga31 = (Ring9.NTri * Ring9.DInSquare) - Ring9.MuVac;
        Ring9 perf496 = (Ring9.C0OverCSq * Ring9.NTri) * ga31;
        Ring9 schottky499 = Ring9.Zero - Ring9.NPhi;
        Ring9 cable500 = ((Ring9.DInSquare * Ring9.DOutEdge) * Ring9.BezoutInv41);
        Ring9 pelt127 = Ring9.QSquare - Ring9.MuVac;
        Ring9 focal508 = Ring9.DInSquare * pelt127;
        Ring9 sio2_30 = (Ring9.NPhi + Ring9.MuVac) + (Ring9.C0OverCSq * Ring9.NTri);
        Ring9 rg510 = sio2_30 * fresnel17;
        Ring9 max511 = Ring9.DOutVertex * (Ring9.KP + Ring9.ModBits);
        Bus36 full512 = Bus36.FromUnsignedRing9(max511) + Bus36.One;
        return (scan444, kp448, fpga456, srf459, qz462, uart468, heis471, core481, co2Sq484, gate487, seam491, heis492, perf496, schottky499, cable500, focal508, rg510, max511, full512);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Bus36 SecondTurnStart513, Bus36 FourteenFaceScan518Ns, Bus36 DoubleAtp520, Bus36 SeamGate525, Bus36 DoubleIvSeed526, Bus36 GalliumFresnelRipple527, Bus36 TwelveMosfetStacks528, Bus36 SilanolSquare529, Bus36 CoreHeisenberg533, Bus36 TotalWallCap540Pf, Bus36 FerroinGate544, Bus36 FourteenthGauss547, Bus36 QuadAlphaFs548, Bus36 SecondTurnHeisenberg553, Bus36 QuadCollapse556, Bus36 EightVerletVdcSq576, Bus36 OxymyoglobinCenter579Nm, Bus36 SixGaussianHeme582, Bus36 EightRadiusVerlet584, Bus36 RadiusCoreStasis585, Bus36 BarrierCollapse587SectionCount, Bus36 ReducedRedFerroin651Nm, Bus36 RedGaAsLed718Nm, Bus36 SpinorFullIdentity720Deg, Bus36 StasisMultiplier833, Bus36 MosfetGateDelay925Ns, Bus36 TwoTurnVerlet1024, Bus36 BezoutGateHeisenberg1025, Bus36 SiliconBandgap1107Nm, Bus36 CavityLoad1536Ohm, Bus36 EsdTotalCap1800Pf, Bus36 PauliDense2048, Bus36 HiddenMicrostates5184, Bus36 VacuumFloor8192)
        ComputeBench587BarrierCollapseMultiscaleClosureBus36()
    {
        Bus36 ring512 = Bus36.FromUnsignedRing9(Ring9.MinusOne) + Bus36.One;
        Bus36 start513 = ring512 + Bus36.One;
        Bus36 scan518 = ring512 + Bus36.FromUnsignedRing9(Ring9.NSq);
        Bus36 atp520 = ring512 + Bus36.FromUnsignedRing9(Ring9.NTri);
        Bus36 gate525 = ring512 + Bus36.FromUnsignedRing9(Ring9.NPhi);
        Bus36 iv526 = ring512 + Bus36.FromUnsignedRing9(Ring9.ZSilicon14);
        Bus36 ripple527 = ring512 + Bus36.FromUnsignedRing9(Ring9.D * Ring9.DOutEdge);
        Bus36 mos528 = ring512 + Bus36.FromUnsignedRing9(Ring9.ZOxygen2_16);
        Ring9 fresnel17 = Ring9.ModBits + Ring9.NTri;
        Bus36 sil529 = ring512 + Bus36.FromUnsignedRing9(fresnel17);
        Bus36 heis533 = ring512 + Bus36.FromUnsignedRing9(Ring9.SeamTicks);
        Bus36 cap540 = ring512 + Bus36.FromUnsignedRing9(Ring9.DInSquare * Ring9.DOutVertex);
        Bus36 fer544 = ring512 + Bus36.FromUnsignedRing9(Ring9.DInSquare * Ring9.NTri);
        Bus36 gauss547 = ring512 + Bus36.FromUnsignedRing9(Ring9.DOutEdge * Ring9.DOutVertex);
        Bus36 alpha548 = ring512 + Bus36.FromUnsignedRing9(Ring9.DInSquare * Ring9.ModBits);
        Bus36 heis553 = ring512 + Bus36.FromUnsignedRing9(Ring9.DPhiDisc);
        Bus36 col556 = ring512 + Bus36.FromUnsignedRing9(Ring9.ZSilicon14 + Ring9.ZQuartz30);
        Bus36 ver576 = ring512 + Bus36.FromUnsignedRing9(Ring9.KP);
        Bus36 oxy579 = ring512 + Bus36.FromUnsignedRing9(Ring9.KP + Ring9.D);
        Bus36 heme582 = ring512 + Bus36.FromUnsignedRing9(Ring9.KP + Ring9.NSq);
        Ring9 verlet72 = Ring9.NTri * Ring9.ModBits;
        Bus36 rad584 = ring512 + Bus36.FromUnsignedRing9(verlet72);
        Bus36 st585 = ring512 + Bus36.FromUnsignedRing9(Ring9.KP + Ring9.ModBits);
        Bus36 col587 = ring512 + Bus36.FromUnsignedRing9(Ring9.D * Ring9.BezoutInv41);
        Bus36 red651 = ring512 + Bus36.FromUnsignedRing9(Ring9.CollapsePhi);
        Bus36 led718 = oxy579 + Bus36.FromUnsignedRing9(Ring9.CollapsePhi);
        Bus36 spin720 = (Bus36.FromUnsignedRing9(verlet72) * (Ring9.NTri + Ring9.C0OverCSq));
        Bus36 st833 = (Bus36.FromUnsignedRing9(Ring9.DOutVertex * Ring9.DOutVertex) * fresnel17);
        Ring9 tick37 = (Ring9.DOutVertex * Ring9.DOutEdge) + Ring9.C0OverCSq;
        Bus36 gate925 = (Bus36.FromUnsignedRing9(Ring9.BezoutInv41) * tick37);
        Bus36 ver1024 = ring512 + ring512;
        Bus36 bez1025 = ver1024 + Bus36.One;
        Bus36 si1107 = (Bus36.FromUnsignedRing9(Ring9.D * Ring9.D * Ring9.D) * Ring9.DPhiDisc);
        Bus36 cav1536 = (ring512 * Ring9.D);
        Bus36 esd1800 = (Bus36.FromUnsignedRing9(Ring9.C0OverCSq * Ring9.ModBits) * (Ring9.DInSquare * Ring9.BezoutInv41));
        Bus36 pauli2048 = (ver1024 + ver1024);
        Bus36 hid5184 = (Bus36.FromUnsignedRing9(verlet72) * verlet72);
        Bus36 eps8192 = Bus36.EpsDenom;
        return (start513, scan518, atp520, gate525, iv526, ripple527, mos528, sil529, heis533, cap540, fer544, gauss547, alpha548, heis553, col556, ver576, oxy579, heme582, rad584, st585, col587, red651, led718, spin720, st833, gate925, ver1024, bez1025, si1107, cav1536, esd1800, pauli2048, hid5184, eps8192);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (
        Bus36 TwelveBeamsMujoco588,
        Bus36 TwelveBeamsTin600,
        Bus36 SecondTurnColdLag601,
        Bus36 RingPlusGlucose608,
        Bus36 Bus36Fresnel612,
        Bus36 TwelveRunMinutes624,
        Bus36 TenKp640,
        Bus36 ReducedRedFerroin651Nm,
        Bus36 BalmerHAlpha656Nm,
        Bus36 EighteenTickArgon666Ns,
        Bus36 IronSquare676,
        Bus36 FourteenFacesTin700,
        Bus36 RedGaAsLed718Nm,
        Bus36 SpinorIdentity720Deg,
        Bus36 SixteenthGauss721,
        Bus36 ModBitsCubed729,
        Bus36 TwentyTickFpga740Ns,
        Bus36 KleinJ744,
        Bus36 SecondTurnUart746,
        Bus36 AtmPressure760Torr,
        Bus36 ThreePauli768,
        Bus36 SecondTurnIvSeed775,
        Bus36 SecondTurnIce785K,
        Bus36 SecondTurnFresnel801,
        Bus36 SeventeenthGauss817,
        Bus36 AudibleFloor823Hz,
        Bus36 AudibleRounded824Hz,
        Bus36 GravityDipoleMult832,
        Bus36 StasisMultiplier833,
        Bus36 SixCollapse834,
        Bus36 CopperSquare841,
        Bus36 TwelveVerlet864,
        Bus36 TwentyFourTick888Ns,
        Bus36 QuartzSquare900,
        Bus36 EighteenthGauss919,
        Bus36 MosfetGateDelay925Ns,
        Bus36 CausalConeT1Bits936,
        Bus36 SecondTurnBlueLed953,
        Bus36 FifteenKp960,
        Bus36 TwelveHeisenbergClass972,
        Bus36 MooreCube27Tick999Ns,
        Bus36 DecimalCube1000,
        Bus36 Scheherazade1001,
        Bus36 FourteenVerlet1008,
        Bus36 EightPeltier1016,
        Bus36 StationaryFirstReturn1023,
        Bus36 TwoTurnVerlet1024,
        Bus36 BezoutGateHeisenberg1025)
        ComputeBench1025TwoTurnVerletBezoutConsecutiveClosureBus36()
    {
        Bus36 ring512 = Bus36.FromUnsignedRing9(Ring9.MinusOne) + Bus36.One;
        Ring9 geoms49 = Ring9.DOutVertex * Ring9.DOutVertex;
        Bus36 n588 = Bus36.FromUnsignedRing9(Ring9.NLinks) * geoms49;
        Ring9 tin50 = Ring9.C0OverCSq * Ring9.BezoutInv41;
        Bus36 n600 = Bus36.FromUnsignedRing9(Ring9.NLinks) * tin50;
        Bus36 n601 = ring512 + Bus36.FromUnsignedRing9(Ring9.KP + Ring9.BezoutInv41);
        Bus36 n608 = ring512 + Bus36.FromUnsignedRing9(Ring9.MHat - Ring9.MuVac);
        Ring9 bus36Bits = Ring9.DInSquare * Ring9.ModBits;
        Ring9 fresnel17 = Ring9.ModBits + Ring9.NTri;
        Bus36 n612 = Bus36.FromUnsignedRing9(bus36Bits) * fresnel17;
        Bus36 n624 = Bus36.FromUnsignedRing9(Ring9.NLinks) * (Ring9.DInSquare * Ring9.NPhi);
        Ring9 water10 = Ring9.NTri + Ring9.C0OverCSq;
        Bus36 n640 = Bus36.FromUnsignedRing9(water10) * Ring9.KP;
        Bus36 n651 = ring512 + Bus36.FromUnsignedRing9(Ring9.CollapsePhi);
        Bus36 n656 = Bus36.FromUnsignedRing9(Ring9.ZOxygen2_16) * Ring9.DPhiDisc;
        Ring9 tick37 = (Ring9.DOutVertex * Ring9.DOutEdge) + Ring9.C0OverCSq;
        Ring9 ar18 = Ring9.C0OverCSq * Ring9.ModBits;
        Bus36 n666 = Bus36.FromUnsignedRing9(ar18) * tick37;
        Ring9 fe26 = Ring9.C0OverCSq * Ring9.NPhi;
        Bus36 n676 = Bus36.FromUnsignedRing9(fe26) * fe26;
        Bus36 n700 = Bus36.FromUnsignedRing9(Ring9.ZSilicon14) * tin50;
        Bus36 n718 = ring512 + Bus36.FromUnsignedRing9(Ring9.QuartzSiO2Phi);
        Ring9 verlet72 = Ring9.NTri * Ring9.ModBits;
        Bus36 n720 = Bus36.FromUnsignedRing9(water10) * verlet72;
        Ring9 oct19 = Ring9.MuVac + Ring9.NLinks + Ring9.NSq;
        Bus36 n721 = ring512 + Bus36.FromUnsignedRing9((Ring9.NLinks - Ring9.MuVac) * oct19);
        Bus36 n729 = Bus36.FromUnsignedRing9(Ring9.ModBits * Ring9.ModBits) * Ring9.ModBits;
        Bus36 n740 = ring512 + Bus36.FromUnsignedRing9(Ring9.NLinks * oct19);
        Ring9 ga31 = (Ring9.NTri * Ring9.DInSquare) - Ring9.MuVac;
        Bus36 n744 = Bus36.FromUnsignedRing9(Ring9.C0OverCSq * Ring9.NLinks) * ga31;
        Bus36 n746 = ring512 + Bus36.FromUnsignedRing9(ar18 * Ring9.NPhi);
        Bus36 n760 = Bus36.FromUnsignedRing9(water10) * (Ring9.DInSquare * oct19);
        Bus36 n768 = Bus36.FromUnsignedRing9(Ring9.PhiPauli) * Ring9.D;
        Bus36 n775 = ring512 + Bus36.FromUnsignedRing9(Ring9.IvSeed);
        Bus36 n785 = ring512 + Bus36.FromUnsignedRing9(Ring9.SeamTicks * Ring9.NPhi);
        Bus36 n801 = ring512 + Bus36.FromUnsignedRing9(fresnel17 * fresnel17);
        Bus36 n817 = ring512 + Bus36.FromUnsignedRing9(Ring9.DOutEdge * (Ring9.KP - Ring9.D));
        Bus36 n823 = ring512 + Bus36.FromUnsignedRing9(Ring9.PhiPauli + (Ring9.DOutEdge * (Ring9.NLinks - Ring9.MuVac)));
        Bus36 n824 = n823 + Bus36.One;
        Bus36 n832 = Bus36.FromUnsignedRing9(Ring9.NPhi) * Ring9.KP;
        Bus36 n833 = n832 + Bus36.One;
        Bus36 n834 = Bus36.FromUnsignedRing9(Ring9.NSq) * Ring9.CollapsePhi;
        Ring9 cu29 = (Ring9.DInSquare * Ring9.DOutVertex) + Ring9.MuVac;
        Bus36 n841 = Bus36.FromUnsignedRing9(cu29) * cu29;
        Bus36 n864 = Bus36.FromUnsignedRing9(Ring9.NLinks) * verlet72;
        Bus36 n888 = Bus36.FromUnsignedRing9(Ring9.C0OverCSq * Ring9.NLinks) * tick37;
        Bus36 n900 = Bus36.FromUnsignedRing9(Ring9.ZQuartz30) * Ring9.ZQuartz30;
        Bus36 n919 = ring512 + Bus36.FromUnsignedRing9((Ring9.NLinks - Ring9.MuVac) * tick37);
        Bus36 n925 = Bus36.FromUnsignedRing9(Ring9.BezoutInv41) * tick37;
        Bus36 n936 = Bus36.FromUnsignedRing9(Ring9.NPhi) * verlet72;
        Bus36 n953 = ring512 + Bus36.FromUnsignedRing9(Ring9.SeamTicks * Ring9.SeamTicks);
        Bus36 n960 = Bus36.FromUnsignedRing9(Ring9.D * Ring9.DOutEdge) * Ring9.KP;
        Ring9 heis81 = (Ring9.NTri * Ring9.ModBits) + Ring9.ModBits;
        Bus36 n972 = Bus36.FromUnsignedRing9(Ring9.NLinks) * heis81;
        Bus36 n999 = Bus36.FromUnsignedRing9(Ring9.D * Ring9.D * Ring9.D) * tick37;
        Bus36 n1000 = Bus36.FromUnsignedRing9(water10) * (Ring9.DInSquare * Ring9.BezoutInv41);
        Bus36 n1001 = Bus36.FromUnsignedRing9(Ring9.DOutVertex * (Ring9.NLinks - Ring9.MuVac)) * Ring9.NPhi;
        Bus36 n1008 = Bus36.FromUnsignedRing9(Ring9.ZSilicon14) * verlet72;
        Bus36 n1016 = Bus36.FromUnsignedRing9(Ring9.NTri) * (Ring9.QSquare - Ring9.MuVac);
        Bus36 n1024 = ring512 + ring512;
        Bus36 n1023 = n1024 - Bus36.One;
        Bus36 n1025 = n1024 + Bus36.One;
        return (n588, n600, n601, n608, n612, n624, n640, n651, n656, n666, n676, n700, n718, n720, n721, n729, n740, n744, n746, n760, n768, n775, n785, n801, n817, n823, n824, n832, n833, n834, n841, n864, n888, n900, n919, n925, n936, n953, n960, n972, n999, n1000, n1001, n1008, n1016, n1023, n1024, n1025);
    }
}
