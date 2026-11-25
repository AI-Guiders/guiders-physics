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
        var (phiNow, dphiVerlet) = KitchenBenchDigitalTwin.SampleSiliconThermalTwoTickPair9();
        SyncWithPhysicalAtom(phiNow, dphiVerlet);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void AdvanceClock()
    {
        ushort isAtomMode = (ushort)ClockMode;
        var (trngPhi, trngDphi) = KitchenBenchDigitalTwin.SampleSiliconThermalTwoTickPair9();
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
