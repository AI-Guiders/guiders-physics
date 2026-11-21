using System.Diagnostics;
using System.Runtime.CompilerServices;
using Guiders.Physics.Core;

namespace Guiders.Physics.Engine;

public sealed class CuboctWorldScene3D
{
    public Ring9 CurrentTick { get; private set; }
    public Ring9 AtomicSlopeV0 { get; private set; }
    public Ring9 MasterAtomPhase { get; private set; }

    public CuboctWorldScene3D()
    {
        CurrentTick = Ring9.Zero;
        AtomicSlopeV0 = Ring9.Zero;
        MasterAtomPhase = Ring9.IvSeed;
    }

    /// <summary>
    /// Locks the Two-Point Worldline Parabola between Origin `P_0 = IV(3) = 263` and a live physical atom/oscillator phase `atomPhase`!
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void SyncWithPhysicalAtom(Ring9 atomPhase)
    {
        MasterAtomPhase = atomPhase;
        var (syncedTick, slopeV0) = CuboctahedronAutomatonKernel.DecodeTwoPointAtomicTimeLock(atomPhase, Ring9.IvSeed);
        CurrentTick = syncedTick;
        AtomicSlopeV0 = slopeV0;
    }

    /// <summary>
    /// Reads the physical quartz crystal oscillator of the host hardware (`Stopwatch.GetTimestamp() &amp; 0x1FF`)
    /// as a live physical atomic phase reference and locks the worldline parabola to it.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void SyncWithHostQuartzCrystal()
    {
        ushort crystalTicks9 = (ushort)(Stopwatch.GetTimestamp() & 0x1FF);
        SyncWithPhysicalAtom(new Ring9(crystalTicks9));
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
