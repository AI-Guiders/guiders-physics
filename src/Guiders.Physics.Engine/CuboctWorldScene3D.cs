using System;
using System.Threading.Tasks;
using Guiders.Physics.Core;

namespace Guiders.Physics.Engine;

/// <summary>
/// Unified 3D World Scene on a Cuboctahedral O_2(3) lattice of 72-bit `VerletCell72` cells.
/// Simulates matter, shockwaves/explosions, light/color, heat, live binaural stereo audio,
/// and exact zero-memory time reversal.
/// </summary>
public sealed class CuboctWorldScene3D
{
    public int SizeZ { get; }
    public int SizeY { get; }
    public int SizeX { get; }
    public long Tick { get; private set; }

    private VerletCell72[,,] _curr;
    private VerletCell72[,,] _next;

    public CuboctWorldScene3D(int sizeZ, int sizeY, int sizeX)
    {
        if (sizeZ < 5 || sizeY < 5 || sizeX < 5)
            throw new ArgumentOutOfRangeException(nameof(sizeZ), "Lattice dimensions must be at least 5x5x5.");
        SizeZ = sizeZ;
        SizeY = sizeY;
        SizeX = sizeX;
        _curr = new VerletCell72[sizeZ, sizeY, sizeX];
        _next = new VerletCell72[sizeZ, sizeY, sizeX];
        InitializeVacuum();
    }

    public long TotalMemoryFootprintBytes => (long)SizeZ * SizeY * SizeX * 9L; // 9 bytes = 72 bits per cell!

    public VerletCell72 GetCell(int z, int y, int x) => _curr[z, y, x];

    public void SetCell(int z, int y, int x, VerletCell72 cell) => _curr[z, y, x] = cell;

    public void InitializeVacuum()
    {
        for (int z = 0; z < SizeZ; z++)
            for (int y = 0; y < SizeY; y++)
                for (int x = 0; x < SizeX; x++)
                    _curr[z, y, x] = VerletCell72.Create(CarrierSpecD3.MuVac, 0, 0, CarrierSpecD3.MuVac, 0, 0);
        Tick = 0;
    }

    /// <summary>
    /// Initializes the 3D scene with the canonical IV(3) = 263 = 0o407 cosmological/nucleation field.
    /// </summary>
    public void InitializeIv263Field()
    {
        int cz = SizeZ / 2, cy = SizeY / 2, cx = SizeX / 2;
        for (int z = 0; z < SizeZ; z++)
            for (int y = 0; y < SizeY; y++)
                for (int x = 0; x < SizeX; x++)
                {
                    int dz = z - cz, dy = y - cy, dx = x - cx;
                    int mu = CuboctahedronAutomatonKernel.Iv263MuAt(dz, dy, dx);
                    int p0 = CuboctahedronAutomatonKernel.Iv263PhiAt(dz, dy, dx, 0).UnsignedTicks;
                    int p1 = CuboctahedronAutomatonKernel.Iv263PhiAt(dz, dy, dx, 1).UnsignedTicks;
                    _curr[z, y, x] = VerletCell72.Create(mu, p0, p0, mu, p1, p1);
                }
        Tick = 0;
    }

    /// <summary>
    /// Triggers a physical detonation / solar emitter impulse at `(cz, cy, cx)` by shifting phase across the Janus saddle `Q_square = 128`.
    /// </summary>
    public void TriggerImpulse(int cz, int cy, int cx, int amplitudeMu = 64, int phaseImpulse = 164)
    {
        var c = _curr[cz, cy, cx];
        _curr[cz, cy, cx] = VerletCell72.Create(
            amplitudeMu,
            (c.Phi0Curr.UnsignedTicks + phaseImpulse) & CarrierSpecD3.RingMask,
            c.Phi0Past.UnsignedTicks,
            amplitudeMu,
            (c.Phi1Curr.UnsignedTicks - phaseImpulse) & CarrierSpecD3.RingMask,
            c.Phi1Past.UnsignedTicks);
    }

    /// <summary>
    /// Advances the entire 3D world forward by 1 discrete Verlet tick (`t -> t + 1`) across all CPU cores.
    /// </summary>
    public void StepForward()
    {
        int sz = SizeZ, sy = SizeY, sx = SizeX;
        Parallel.For(0, sz, z =>
        {
            Span<int> locMu0 = stackalloc int[12];
            Span<int> locMu1 = stackalloc int[12];
            Span<Ring9> locP0 = stackalloc Ring9[12];
            Span<Ring9> locP1 = stackalloc Ring9[12];

            for (int y = 0; y < sy; y++)
                for (int x = 0; x < sx; x++)
                {
                    var c = _curr[z, y, x];
                    for (int k = 0; k < 12; k++)
                    {
                        var (dz, dy, dx) = CarrierSpecD3.CuboctOffsets12[k];
                        int nz = (z - dz + sz) % sz;
                        int ny = (y - dy + sy) % sy;
                        int nx = (x - dx + sx) % sx;
                        var nb = _curr[nz, ny, nx];
                        locMu0[k] = nb.Mu0.UnsignedTicks;
                        locMu1[k] = nb.Mu1.UnsignedTicks;
                        locP0[k] = nb.Phi0Curr;
                        locP1[k] = nb.Phi1Curr;
                    }
                    var g0 = CuboctahedronAutomatonKernel.Compute12NeighborGate(c.Mu0.UnsignedTicks, c.Phi0Curr, locMu0, locP0, true);
                    var g1 = CuboctahedronAutomatonKernel.Compute12NeighborGate(c.Mu1.UnsignedTicks, c.Phi1Curr, locMu1, locP1, true);
                    Ring9 gLead = Math.Abs(g1.GatePhi.SignedWrap) > Math.Abs(g0.GatePhi.SignedWrap) ? g1.GatePhi : g0.GatePhi;
                    _next[z, y, x] = CuboctahedronAutomatonKernel.StepCellForward(c, gLead, gLead);
                }
        });
        (_curr, _next) = (_next, _curr);
        Tick++;
    }

    /// <summary>
    /// Reverses the entire 3D world backward in time by 1 discrete Verlet tick (`t + 1 -> t`) with ZERO history stack!
    /// Uses only the current two-tick cell state `VerletCell72` and recomputes the holonomy gate at `t`.
    /// </summary>
    public void StepBackward()
    {
        int sz = SizeZ, sy = SizeY, sx = SizeX;
        Parallel.For(0, sz, z =>
        {
            Span<int> locMu0 = stackalloc int[12];
            Span<int> locMu1 = stackalloc int[12];
            Span<Ring9> locP0AtT = stackalloc Ring9[12];
            Span<Ring9> locP1AtT = stackalloc Ring9[12];

            for (int y = 0; y < sy; y++)
                for (int x = 0; x < sx; x++)
                {
                    var c = _curr[z, y, x];
                    for (int k = 0; k < 12; k++)
                    {
                        var (dz, dy, dx) = CarrierSpecD3.CuboctOffsets12[k];
                        int nz = (z - dz + sz) % sz;
                        int ny = (y - dy + sy) % sy;
                        int nx = (x - dx + sx) % sx;
                        var nb = _curr[nz, ny, nx];
                        locMu0[k] = nb.Mu0.UnsignedTicks;
                        locMu1[k] = nb.Mu1.UnsignedTicks;
                        // At tick t+1, PhiPast holds the exact phase at tick t!
                        locP0AtT[k] = nb.Phi0Past;
                        locP1AtT[k] = nb.Phi1Past;
                    }
                    var g0 = CuboctahedronAutomatonKernel.Compute12NeighborGate(c.Mu0.UnsignedTicks, c.Phi0Past, locMu0, locP0AtT, true);
                    var g1 = CuboctahedronAutomatonKernel.Compute12NeighborGate(c.Mu1.UnsignedTicks, c.Phi1Past, locMu1, locP1AtT, true);
                    Ring9 gLead = Math.Abs(g1.GatePhi.SignedWrap) > Math.Abs(g0.GatePhi.SignedWrap) ? g1.GatePhi : g0.GatePhi;
                    _next[z, y, x] = CuboctahedronAutomatonKernel.StepCellReverse(c, gLead, gLead);
                }
        });
        (_curr, _next) = (_next, _curr);
        Tick--;
    }

    /// <summary>
    /// Reads live binaural stereo 16-bit PCM audio sample `(LeftEar, RightEar)` directly from two listener ear cells in the 3D scene
    /// (for direct streaming to Windows sound card / headphones without pre-recorded audio files).
    /// </summary>
    public (short LeftPcm16, short RightPcm16) SampleBinauralAudioListener(
        (int Z, int Y, int X) leftEarPos, (int Z, int Y, int X) rightEarPos)
    {
        var sL = CuboctahedronAutomatonKernel.ReadCellSensors(_curr[leftEarPos.Z, leftEarPos.Y, leftEarPos.X]);
        var sR = CuboctahedronAutomatonKernel.ReadCellSensors(_curr[rightEarPos.Z, rightEarPos.Y, rightEarPos.X]);
        short pcmL = (short)Math.Clamp(sL.AcousticPressure - 128, short.MinValue, short.MaxValue);
        short pcmR = (short)Math.Clamp(sR.AcousticPressure - 128, short.MinValue, short.MaxValue);
        return (pcmL, pcmR);
    }

    /// <summary>
    /// Reads live optical, thermal, and acoustic sensor state at a camera/listener voxel `(z, y, x)`.
    /// </summary>
    public MultimodalReadout SampleCameraVoxel(int z, int y, int x) =>
        CuboctahedronAutomatonKernel.ReadCellSensors(_curr[z, y, x]);
}
