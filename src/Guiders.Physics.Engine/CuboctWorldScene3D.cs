using System.Runtime.CompilerServices;
using Guiders.Physics.Core;

namespace Guiders.Physics.Engine;

/// <summary>
/// Unified 3D World Engine on Pure 9-Bit `Ring9` and 72-Bit `VerletCell72`.
/// Zero `int`, zero `if`, zero 3D grid loops — pure closed-form O(1) Verlet worldline parabola on `Ring9`!
/// </summary>
public sealed class CuboctWorldScene3D
{
    public Ring9 CurrentTick { get; private set; }

    public CuboctWorldScene3D() => CurrentTick = Ring9.Zero;

    public void JumpToTick(Ring9 tick) => CurrentTick = tick;

    public void StepForward() => CurrentTick = CurrentTick + Ring9.MuVac;

    public void StepBackward() => CurrentTick = CurrentTick - Ring9.MuVac;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public VerletCell72 EvaluateWorldlineAt(Ring9 dz, Ring9 dy, Ring9 dx, Ring9 t) =>
        CuboctahedronAutomatonKernel.EvaluateClosedFormWorldlineCellAt(dz, dy, dx, t);

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
