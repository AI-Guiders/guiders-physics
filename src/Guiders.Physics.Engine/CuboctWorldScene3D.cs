using System.Collections.Generic;
using System.Diagnostics;
using System.Runtime.CompilerServices;
using Guiders.Physics.Core;

namespace Guiders.Physics.Engine;

public enum WorldClockMode : byte
{
    FreeRunning = 0,
    AtomSynced  = 1
}

public readonly record struct WorldEvent9(
    Ring9 Z,
    Ring9 Y,
    Ring9 X,
    Ring9 BirthTick,
    Ring9 Mu,
    Ring9 GateImpulse);

public sealed class CuboctWorldScene3D
{
    private readonly List<WorldEvent9> _events = [];

    public WorldClockMode ClockMode { get; set; }
    public Ring9 CurrentTick { get; private set; }
    public Ring9 AtomicSlopeV0 { get; private set; }
    public Ring9 MasterAtomPhase { get; private set; }
    public IReadOnlyList<WorldEvent9> Events => _events;

    public CuboctWorldScene3D(WorldClockMode clockMode = WorldClockMode.FreeRunning)
    {
        ClockMode = clockMode;
        CurrentTick = Ring9.Zero;
        AtomicSlopeV0 = Ring9.Zero;
        MasterAtomPhase = Ring9.IvSeed;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void SyncWithPhysicalAtom(Ring9 atomPhase)
    {
        MasterAtomPhase = atomPhase;
        var (syncedTick, slopeV0) = CuboctahedronAutomatonKernel.DecodeTwoPointAtomicTimeLock(atomPhase, Ring9.IvSeed);
        CurrentTick = syncedTick;
        AtomicSlopeV0 = slopeV0;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void SyncWithHostQuartzCrystal()
    {
        ushort crystalTicks9 = (ushort)(Stopwatch.GetTimestamp() & 0x1FF);
        SyncWithPhysicalAtom(new Ring9(crystalTicks9));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void AdvanceClock()
    {
        ushort isAtomMode = (ushort)ClockMode;
        ushort crystalTicks9 = (ushort)(Stopwatch.GetTimestamp() & 0x1FF);
        var (syncedTick, slopeV0) = CuboctahedronAutomatonKernel.DecodeTwoPointAtomicTimeLock(new Ring9(crystalTicks9), Ring9.IvSeed);
        Ring9 freeNext = CurrentTick + Ring9.MuVac;

        Ring9 selAtom = new(isAtomMode);
        Ring9 selFree = Ring9.MuVac - selAtom;

        CurrentTick = selFree * freeNext + selAtom * syncedTick;
        AtomicSlopeV0 = selFree * AtomicSlopeV0 + selAtom * slopeV0;
    }

    public void JumpToTick(Ring9 tick) => CurrentTick = tick;

    public void StepForward() => CurrentTick = CurrentTick + Ring9.MuVac;

    public void StepBackward() => CurrentTick = CurrentTick - Ring9.MuVac;

    public void EmitGameEvent(Ring9 z, Ring9 y, Ring9 x, Ring9 mu, Ring9 gateImpulse) =>
        _events.Add(new WorldEvent9(z, y, x, CurrentTick, mu, gateImpulse));

    public void ClearGameEvents() => _events.Clear();

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public VerletCell72 EvaluateWorldlineAt(Ring9 dz, Ring9 dy, Ring9 dx, Ring9 t)
    {
        var cell = CuboctahedronAutomatonKernel.EvaluateClosedFormWorldlineCellAt(dz, dy, dx, t, AtomicSlopeV0);
        Ring9 p0 = cell.Phi0Curr;
        Ring9 p1 = cell.Phi1Curr;
        Ring9 mu0 = cell.Mu0;
        Ring9 mu1 = cell.Mu1;

        foreach (var ev in _events)
        {
            Ring9 dt = t - ev.BirthTick;
            Ring9 dist = (dz - ev.Z).AbsWrap + (dy - ev.Y).AbsWrap + (dx - ev.X).AbsWrap;
            Ring9 wavePhase = dist.JumpVerletParabola(ev.Mu, ev.GateImpulse, dt);
            Ring9 gatedWave = CuboctahedronAutomatonKernel.ApplyHeisenbergDeadband(wavePhase);
            p0 = p0 + gatedWave;
            p1 = p1 - gatedWave;
        }

        var outCell = VerletCell72.Create(mu0, p0, cell.Phi0Past, mu1, p1, cell.Phi1Past);
        outCell.Gate0 = cell.Gate0;
        outCell.Gate1 = cell.Gate1;
        return outCell;
    }

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
