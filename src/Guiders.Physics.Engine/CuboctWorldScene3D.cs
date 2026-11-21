using System;
using System.Runtime.CompilerServices;
using Guiders.Physics.Core;

namespace Guiders.Physics.Engine;

/// <summary>
/// Unified 3D World Engine on the 9-Bit `Ring9` / 72-Bit `VerletCell72` Cuboctahedron O_2(3) Carrier.
/// Features:
/// - Coordinate-Free Closed-Form O(1) Worldline Trajectory `EvaluateWorldlineAt(dz, dy, dx, t)`:
///   evaluates ANY cell `(dz, dy, dx) in Z_512^3` at ANY past or future time tick `t` in O(1) via the
///   Verlet parabola `phi(t) = phi(0) + t*v0 + (t*(t+1)/2)*Gate` with ZERO 3D grid memory allocation!
/// - Live Multimodal Sensor Readout at any `(dz, dy, dx, t)`: optical wavelength [440..718 nm],
///   RGB triad (R+G+B = 128), binaural stereo acoustic pressure (LeftEar, RightEar), and thermal quanta n_E.
/// - Exact 3D Topological Vortex Plaquette Birth Census (`86 = 38 + 16 + 32` vortex particles born at `t = 2` from `IV(3) = 263`).
/// </summary>
public sealed class CuboctWorldScene3D
{
    public int CurrentTick { get; private set; }

    // Optional local perturbation impulse at center (dz=0, dy=0, dx=0)
    private Ring9 _impulseMu;
    private Ring9 _impulsePhase;
    private int _hasImpulseFlag;

    public CuboctWorldScene3D()
    {
        CurrentTick = 0;
        _hasImpulseFlag = 0;
    }

    public void JumpToTick(int tick) => CurrentTick = tick;

    public void StepForward() => CurrentTick++;

    public void StepBackward() => CurrentTick--;

    public void TriggerCenterImpulse(Ring9 amplitudeMu, Ring9 phaseImpulse)
    {
        _impulseMu = amplitudeMu;
        _impulsePhase = phaseImpulse;
        _hasImpulseFlag = 1;
    }

    /// <summary>
    /// Closed-form O(1) evaluation of the 72-bit `VerletCell72` at any 3D coordinate `(dz, dy, dx)` in `Ring9^3`
    /// and any time tick `t` (past, present, or future!) along the Verlet worldline parabola!
    /// Zero 3D array allocation, zero `if` branches, zero time-stepping loops!
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public VerletCell72 EvaluateWorldlineAt(Ring9 dz, Ring9 dy, Ring9 dx, int t)
    {
        var cell = CuboctahedronAutomatonKernel.EvaluateClosedFormWorldlineCellAt(dz, dy, dx, t);
        int isOrigin = CuboctahedronAutomatonKernel.BranchlessIsZeroNat(
            CuboctahedronAutomatonKernel.BranchlessAbs(dz.SignedWrap) +
            CuboctahedronAutomatonKernel.BranchlessAbs(dy.SignedWrap) +
            CuboctahedronAutomatonKernel.BranchlessAbs(dx.SignedWrap));
        int applyImp = _hasImpulseFlag * isOrigin;
        Ring9 p0 = cell.Phi0Curr + _impulsePhase * applyImp;
        Ring9 p1 = cell.Phi1Curr - _impulsePhase * applyImp;
        Ring9 mu = new(CuboctahedronAutomatonKernel.BranchlessBlend(applyImp, _impulseMu.Bits9, cell.Mu0.Bits9));
        var res = VerletCell72.Create(mu, p0, cell.Phi0Past, mu, p1, cell.Phi1Past);
        res.Gate0 = cell.Gate0;
        res.Gate1 = cell.Gate1;
        return res;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public VerletCell72 GetCellAtCurrentTick(Ring9 dz, Ring9 dy, Ring9 dx) =>
        EvaluateWorldlineAt(dz, dy, dx, CurrentTick);

    /// <summary>
    /// Samples live optical wavelength [440..718 nm], RGB triad (R+G+B=128), acoustic pressure,
    /// and thermal quanta at any `(dz, dy, dx)` at the current world tick in O(1)!
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public MultimodalReadout SampleCameraVoxel(Ring9 dz, Ring9 dy, Ring9 dx) =>
        CuboctahedronAutomatonKernel.ReadCellSensors(GetCellAtCurrentTick(dz, dy, dx));

    /// <summary>
    /// Samples live binaural stereo 16-bit PCM audio `(LeftEar, RightEar)` directly from the closed-form worldline at `CurrentTick`!
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public (short LeftPcm16, short RightPcm16) SampleBinauralAudioListener(
        (Ring9 Dz, Ring9 Dy, Ring9 Dx) leftEar, (Ring9 Dz, Ring9 Dy, Ring9 Dx) rightEar)
    {
        var sL = SampleCameraVoxel(leftEar.Dz, leftEar.Dy, leftEar.Dx);
        var sR = SampleCameraVoxel(rightEar.Dz, rightEar.Dy, rightEar.Dx);
        return ((short)(sL.AcousticPressure - 128), (short)(sR.AcousticPressure - 128));
    }
}
