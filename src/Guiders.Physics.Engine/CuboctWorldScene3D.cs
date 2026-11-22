using System.Runtime.CompilerServices;
using Guiders.Physics.Core;

namespace Guiders.Physics.Engine;

public enum WorldClockMode : byte
{
    FreeRunning = 0,
    AtomSynced  = 1
}

public sealed class CuboctWorldScene3D
{
    public WorldClockMode ClockMode { get; set; }
    public Ring9 CurrentTick { get; private set; }
    public Ring9 AtomicSlopeV0 { get; private set; }
    public Ring9 MasterAtomPhase { get; private set; }
    public Ring9 SpatialAnchorX { get; private set; }

    public CuboctWorldScene3D(WorldClockMode clockMode = WorldClockMode.FreeRunning)
    {
        ClockMode = clockMode;
        CurrentTick = Ring9.Zero;
        AtomicSlopeV0 = Ring9.Zero;
        MasterAtomPhase = Ring9.IvSeed;
        SpatialAnchorX = Ring9.Zero;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void SyncWithPhysicalAtom(Ring9 phiAtomNow, Ring9 deltaPhiVerletNow = default)
    {
        MasterAtomPhase = phiAtomNow;
        SpatialAnchorX = CuboctahedronAutomatonKernel.DecodeAtomicSpatialCoord(phiAtomNow - Ring9.IvSeed);
        CurrentTick = CuboctahedronAutomatonKernel.DecodeAtomicTimeTick(deltaPhiVerletNow);
        AtomicSlopeV0 = CuboctahedronAutomatonKernel.TwoPointWorldlineSlopeRecovery(Ring9.IvSeed, phiAtomNow);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void SyncWithSiliconThermalTrng()
    {
        var (phiNow, dphiVerlet) = CuboctahedronAutomatonKernel.SampleSiliconThermalTwoTickPair9();
        SyncWithPhysicalAtom(phiNow, dphiVerlet);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void AdvanceClock()
    {
        ushort isAtomMode = (ushort)ClockMode;
        var (trngPhi, trngDphi) = CuboctahedronAutomatonKernel.SampleSiliconThermalTwoTickPair9();
        Ring9 syncedTick = CuboctahedronAutomatonKernel.DecodeAtomicTimeTick(trngDphi);
        Ring9 slopeV0 = CuboctahedronAutomatonKernel.TwoPointWorldlineSlopeRecovery(Ring9.IvSeed, trngPhi);
        Ring9 freeNext = CurrentTick + Ring9.MuVac;

        Ring9 selAtom = new(isAtomMode);
        Ring9 selFree = Ring9.MuVac - selAtom;

        CurrentTick = selFree * freeNext + selAtom * syncedTick;
        AtomicSlopeV0 = selFree * AtomicSlopeV0 + selAtom * slopeV0;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void AdvanceClock(Ring9 externalAtomPhase, Ring9 externalDeltaPhiVerlet)
    {
        ushort isAtomMode = (ushort)ClockMode;
        Ring9 syncedTick = CuboctahedronAutomatonKernel.DecodeAtomicTimeTick(externalDeltaPhiVerlet);
        Ring9 slopeV0 = CuboctahedronAutomatonKernel.TwoPointWorldlineSlopeRecovery(Ring9.IvSeed, externalAtomPhase);
        Ring9 freeNext = CurrentTick + Ring9.MuVac;

        Ring9 selAtom = new(isAtomMode);
        Ring9 selFree = Ring9.MuVac - selAtom;

        CurrentTick = selFree * freeNext + selAtom * syncedTick;
        AtomicSlopeV0 = selFree * AtomicSlopeV0 + selAtom * slopeV0;
    }

    public void JumpToTick(Ring9 tick) => CurrentTick = tick;

    public void StepForward() => CurrentTick = CurrentTick + Ring9.MuVac;

    public void StepBackward() => CurrentTick = CurrentTick - Ring9.MuVac;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public VerletCell72 EvaluateWorldlineAt(Ring9 dz, Ring9 dy, Ring9 dx, Ring9 t) =>
        CuboctahedronAutomatonKernel.EvaluateClosedFormWorldlineCellAt(dz, dy, dx, t, AtomicSlopeV0);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public VerletCell72 GetCellAtCurrentTick(Ring9 dz, Ring9 dy, Ring9 dx) =>
        EvaluateWorldlineAt(dz, dy, dx, CurrentTick);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public MultimodalReadout9 SampleCameraVoxel(Ring9 dz, Ring9 dy, Ring9 dx) =>
        CuboctahedronAutomatonKernel.ReadCellSensors(GetCellAtCurrentTick(dz, dy, dx));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public (Ring9 LeftEarWave, Ring9 RightEarWave) SampleBinauralAudioListener(
        (Ring9 Dz, Ring9 Dy, Ring9 Dx) leftEar, (Ring9 Dz, Ring9 Dy, Ring9 Dx) rightEar)
    {
        var sL = SampleCameraVoxel(leftEar.Dz, leftEar.Dy, leftEar.Dx);
        var sR = SampleCameraVoxel(rightEar.Dz, rightEar.Dy, rightEar.Dx);
        return (sL.AcousticWave, sR.AcousticWave);
    }
}

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
/// Full Digital Twin of the 12-Beam Cuboctahedron O_2(3) Kitchen Bench (`Tang Nano 9K` 27.000 MHz)
/// inside `Guiders.Physics.Engine`:
/// - Models all 12 hand-wound copper coils (`Z(Cu) = 29`, `N = 60` turns, `MMF = 42 = S_2(3)` At, Wheeler denom `202`),
/// - Models MOSFET gate stack (`Z(Si+SiO2) = 44`, Miller plateau lag `25 = d_out(E)^2` ticks),
/// - Models 1N5819 Schottky flyback relaxation (`3 * 81 + 13 = 256` ticks),
/// - Models quartz Petri dish wall (`Z(SiO2) = 30 == Z(N2+O2)`, Fresnel 9:8 transmission `288/289`),
/// - Generates and verifies 16-byte XOR-CRC UART frames (`0x06`, `0x08`, `0x02`),
/// - Executes the 52-minute (`52 = 4 * 13 = Z(Te)`) 6-stage synthesis schedule (`Stage 0..5`) with 0-tick phase error!
/// </summary>
public sealed class KitchenBenchDigitalTwin
{
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

    public static bool VerifyUart16ByteFrameCrc(byte[] frame)
    {
        if (frame is null || frame.Length != 16 || frame[0] != 0xA5 || frame[15] != 0x5A)
            return false;
        byte crc = 0;
        for (int i = 1; i <= 13; i++)
            crc ^= frame[i];
        return crc == frame[14];
    }

    public KitchenBenchStageSnapshot9 EvaluateMinute(int minute)
    {
        int m = minute < 0 ? 0 : (minute > 52 ? 52 : minute);
        Ring9 stageIdx;
        Ring9 targetDphi;
        Ring9 targetZ;

        if (m < 4)
        {
            stageIdx = Ring9.Zero;
            targetDphi = Ring9.CollapsePhi; // 139
            targetZ = Ring9.NTri + Ring9.C0OverCSq; // 10 (H2O condensation on TEC1-12706)
        }
        else if (m < 18)
        {
            stageIdx = Ring9.MuVac;
            targetDphi = Ring9.CollapsePhi; // 139
            targetZ = Ring9.MHat - Ring9.MuVac; // 96 (Glucose + 6 O2 photosynthesis)
        }
        else if (m < 34)
        {
            stageIdx = Ring9.C0OverCSq;
            targetDphi = Ring9.CollapsePhi; // 139
            targetZ = Ring9.KP + Ring9.NSq; // 70 (Chiral L-amino & Heme FeN4O2)
        }
        else if (m < 46)
        {
            stageIdx = Ring9.D;
            targetDphi = Ring9.QSquare; // 128 (12-beam 52.734 kHz ultrasonic folding)
            targetZ = Ring9.KP + Ring9.NSq; // 70
        }
        else if (m < 52)
        {
            stageIdx = Ring9.DInSquare;
            targetDphi = Ring9.CollapsePhi; // 139 (68 C Maillard golden crust)
            targetZ = Ring9.KP + Ring9.NSq; // 70
        }
        else
        {
            stageIdx = Ring9.DOutEdge;
            targetDphi = Ring9.Zero; // 0 (Phase stasis lock D = 6,823,936)
            targetZ = Ring9.KP + Ring9.NSq; // 70
        }

        // Thermal coil lag shifts from 64 ticks (25 C) to 58 ticks (63.8 C at 52 min):
        ushort thermalDrop = (ushort)((m * 6) / 52);
        Ring9 coilLag = Ring9.KP - new Ring9(thermalDrop);
        Ring9 totalLag = coilLag + GateLagTicks + DishWallLagTicks;
        Ring9 fpgaDrive = CuboctahedronAutomatonKernel.ComputeKitchenBenchPrecompensatedBeamPhase9(
            targetDphi, GateLagTicks, coilLag, Ring9.Zero, DishWallLagTicks, Ring9.Zero);
        Ring9 arrived = CuboctahedronAutomatonKernel.ComputeKitchenBenchArrivedDishPhase9(
            fpgaDrive, GateLagTicks, coilLag, Ring9.Zero, DishWallLagTicks, Ring9.Zero);

        var dishCell = VerletCell72.Create(Ring9.MuVac, arrived, Ring9.Zero, Ring9.MuVac, Ring9.Zero, Ring9.Zero);
        var telemetry = CuboctahedronAutomatonKernel.ReadCellSensors(dishCell);
        Ring9 recovered = CuboctahedronAutomatonKernel.RecoverRingPhaseFromHydrophoneAndChiralCurrent9(
            Ring9.MuVac, Ring9.MuVac, telemetry.AcousticWave, telemetry.ReactionCurrent);

        uint lagInt = 0;
        for (ushort cand = 0; cand < 512; cand++)
        {
            if (new Ring9(cand) == totalLag)
            {
                lagInt = cand;
                break;
            }
        }
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
        Ring9 recDish0 = CuboctahedronAutomatonKernel.RecoverRingPhaseFromHydrophoneAndChiralCurrent9(
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
}
