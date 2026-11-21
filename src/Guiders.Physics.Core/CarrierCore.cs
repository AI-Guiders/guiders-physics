using System;
using System.Diagnostics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace Guiders.Physics.Core;

/// <summary>
/// Complete passport of the M-layer carrier derived purely from spatial dimension D = 3.
/// All ring constants live natively in 9 bits (Z_512).
/// </summary>
public static class CarrierSpecD3
{
    public const int D = 3;
    public const int MuVac = (D * D) - (1 << D);              // 1
    public const int C0OverCSq = D - MuVac;                   // 2
    public const int NLinks = 2 * D * (D - MuVac);            // 12
    public const int NTri = 1 << D;                           // 8
    public const int NSq = 2 * D;                             // 6
    public const int DInSquare = NLinks / D;                  // 4
    public const int DOutEdge = NSq - MuVac;                  // 5
    public const int DOutVertex = NTri - MuVac;               // 7
    public const int ModBits = MuVac + NTri;                  // 9
    public const int FracBits = ModBits - D;                  // 6
    public const int NRing = 1 << ModBits;                    // 512
    public const int RingMask = NRing - MuVac;                // 511 = 0x1FF
    public const int KP = 1 << FracBits;                      // 64
    public const int QSquare = NRing / DInSquare;             // 128
    public const int PhiPauli = NRing / C0OverCSq;            // 256
    public const int EpsDenom = KP * QSquare;                 // 8192
    public const int NPhi = MuVac + NLinks;                   // 13
    public const int DPhiDisc = (NRing + NPhi + NTri) / NPhi; // 41
    public const int BezoutInv41 = 25;                        // 25 * 41 === 1 (mod 512)
    public const int SeamTicks = NPhi * DPhiDisc - NRing;     // 21
    public const int MHat = MuVac + NLinks * NTri;            // 97
    public const int CollapsePhi = QSquare + NLinks - MuVac;  // 139
    public const int IvSeed = PhiPauli + DOutVertex;          // 263 = 0o407
    public const int B1SingleTickBits = DInSquare * ModBits;  // 36 bits
    public const int B2TwoTickBits = NTri * ModBits;          // 72 bits = 9 bytes!

    // Canonical 12 FCC Cuboctahedron O_2(3) neighbor offsets in Ring9
    public static readonly (Ring9 Dz, Ring9 Dy, Ring9 Dx)[] CuboctOffsets12 =
    [
        (new(0),   new(1),   new(1)),   (new(0),   new(1),   new(511)),
        (new(0),   new(511), new(1)),   (new(0),   new(511), new(511)),
        (new(1),   new(0),   new(1)),   (new(1),   new(0),   new(511)),
        (new(511), new(0),   new(1)),   (new(511), new(0),   new(511)),
        (new(1),   new(1),   new(0)),   (new(1),   new(511), new(0)),
        (new(511), new(1),   new(0)),   (new(511), new(511), new(0))
    ];
}

/// <summary>
/// Native 9-Bit Modular Ring Element in Z_512 (bits [8:0]):
/// - 100% branchless Gray-Quadrant Butterfly Rotator (C_3, S_3)
/// - Closed-form O(1) Verlet Worldline Parabola `JumpVerletParabola(v0, gatePhi, t)`
/// - Bezout Coordinate Decoder `* 25 mod 512` (`41^-1 === 25 mod 512`)
/// </summary>
public readonly record struct Ring9
{
    private readonly ushort _bits9; // Strictly 9 bits [0..511]

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public Ring9(int value) => _bits9 = (ushort)(value & CarrierSpecD3.RingMask);

    public ushort Bits9
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _bits9;
    }

    public int UnsignedTicks
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _bits9;
    }

    /// <summary>Branchless signed representative in (-256, +256].</summary>
    public int SignedWrap
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => CarrierSpecD3.PhiPauli - ((CarrierSpecD3.PhiPauli - _bits9) & CarrierSpecD3.RingMask);
    }

    /// <summary>
    /// 100% branchless Gray-Quadrant Butterfly projections (C_3(phi), S_3(phi)) on 9-bit wires:
    /// bit 8 = q1 (half-turn 256), bit 7 = q0 (quadrant 128), 7-bit mantissa m7 = phi[6:0].
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public (int ProjC, int ProjS) GrayButterflyRotator()
    {
        int q1 = (_bits9 >> 8) & 1;
        int q0 = (_bits9 >> 7) & 1;
        int m7 = _bits9 & 127;
        int comp7 = CarrierSpecD3.QSquare - m7;
        int magC = (1 - q0) * comp7 + q0 * m7;
        int magS = (1 - q0) * m7 + q0 * comp7;
        int sgnC = 1 - 2 * (q1 ^ q0);
        int sgnS = 1 - 2 * q1;
        return (sgnC * magC, sgnS * magS);
    }

    public int ProjC
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => GrayButterflyRotator().ProjC;
    }

    public int ProjS
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => GrayButterflyRotator().ProjS;
    }

    /// <summary>
    /// Closed-form O(1) Verlet Worldline Parabola Jump on Z_512 (`verletParabolaClosedFormJumpOp`):
    /// Jumps directly from initial phase `this = phi(0)` to ANY arbitrary time tick `t` (forward or backward!)
    /// in O(1) without stepping through intermediate ticks:
    /// `phi(t) = (phi(0) + t * v0 + (t * (t + 1) / 2) * gatePhi) mod 512`.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public Ring9 JumpVerletParabola(Ring9 v0, Ring9 gatePhi, int t)
    {
        int triT = (t * (t + 1)) >> 1;
        return new Ring9(_bits9 + t * v0.SignedWrap + triT * gatePhi.SignedWrap);
    }

    /// <summary>
    /// Coordinate-free Bezout inverse of Heisenberg quantum 41 (`25 * 41 === 1 mod 512`).
    /// Decodes spatial coordinate or time tick in 1 modular multiplication on Ring9!
    /// </summary>
    public Ring9 DecodeBezout41
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new(_bits9 * CarrierSpecD3.BezoutInv41);
    }

    public (int O2, int O1, int O0) OctalTriad => ((_bits9 >> 6) & 7, (_bits9 >> 3) & 7, _bits9 & 7);

    public Ring9 PauliPartner
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new(_bits9 ^ CarrierSpecD3.PhiPauli);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 operator +(Ring9 a, Ring9 b) => new(a._bits9 + b._bits9);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 operator -(Ring9 a, Ring9 b) => new(a._bits9 - b._bits9);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 operator *(Ring9 a, Ring9 b) => new(a._bits9 * b._bits9);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 operator *(Ring9 a, int k) => new(a._bits9 * k);

    public override string ToString()
    {
        var (o2, o1, o0) = OctalTriad;
        return $"{_bits9} (0o{o2}{o1}{o0}, wrap={SignedWrap:+#;-#;0})";
    }
}

/// <summary>
/// Exact 72-Bit (9-Byte, Pack = 1) Two-Tick Verlet Polar Cell:
/// packs eight 9-bit `Ring9` lanes (Mu0, Phi0Curr, Phi0Past, Mu1, Phi1Curr, Phi1Past, Gate0, Gate1)
/// into literally 9 bytes = 72 bits (`B_2(3) = 8 * 9 = 72 bits`) with ZERO padding and ZERO branches!
/// </summary>
[StructLayout(LayoutKind.Sequential, Pack = 1, Size = 9)]
public struct VerletCell72
{
    private ulong _lo64;
    private byte _hi8;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public readonly Ring9 GetLane9(int laneIdx)
    {
        UInt128 w = (UInt128)_lo64 | ((UInt128)_hi8 << 64);
        int bits = (int)((w >> (laneIdx * CarrierSpecD3.ModBits)) & (UInt128)CarrierSpecD3.RingMask);
        return new Ring9(bits);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void SetLane9(int laneIdx, Ring9 val)
    {
        int shift = laneIdx * CarrierSpecD3.ModBits;
        UInt128 w = (UInt128)_lo64 | ((UInt128)_hi8 << 64);
        UInt128 mask = ~((UInt128)CarrierSpecD3.RingMask << shift);
        w = (w & mask) | ((UInt128)val.Bits9 << shift);
        _lo64 = (ulong)w;
        _hi8 = (byte)(w >> 64);
    }

    public Ring9 Mu0      { readonly get => GetLane9(0); set => SetLane9(0, value); }
    public Ring9 Phi0Curr { readonly get => GetLane9(1); set => SetLane9(1, value); }
    public Ring9 Phi0Past { readonly get => GetLane9(2); set => SetLane9(2, value); }
    public Ring9 Mu1      { readonly get => GetLane9(3); set => SetLane9(3, value); }
    public Ring9 Phi1Curr { readonly get => GetLane9(4); set => SetLane9(4, value); }
    public Ring9 Phi1Past { readonly get => GetLane9(5); set => SetLane9(5, value); }
    public Ring9 Gate0    { readonly get => GetLane9(6); set => SetLane9(6, value); }
    public Ring9 Gate1    { readonly get => GetLane9(7); set => SetLane9(7, value); }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static VerletCell72 Create(Ring9 mu0, Ring9 p0Curr, Ring9 p0Past, Ring9 mu1, Ring9 p1Curr, Ring9 p1Past)
    {
        VerletCell72 c = default;
        c.Mu0 = mu0;
        c.Phi0Curr = p0Curr;
        c.Phi0Past = p0Past;
        c.Mu1 = mu1;
        c.Phi1Curr = p1Curr;
        c.Phi1Past = p1Past;
        return c;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static VerletCell72 Create(int mu0, int p0Curr, int p0Past, int mu1, int p1Curr, int p1Past) =>
        Create(new Ring9(mu0), new Ring9(p0Curr), new Ring9(p0Past), new Ring9(mu1), new Ring9(p1Curr), new Ring9(p1Past));

    /// <summary>
    /// Closed-form O(1) Worldline Evaluation of this 72-bit cell at ANY arbitrary time tick `t`!
    /// Uses the stationary-gate Verlet parabola `phi(t) = phi(0) + (t*(t+1)/2)*Gate` directly on `Ring9`!
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public readonly VerletCell72 EvaluateClosedFormAtTick(int t)
    {
        Ring9 v0_0 = Phi0Curr - Phi0Past;
        Ring9 v0_1 = Phi1Curr - Phi1Past;
        Ring9 p0AtT = Phi0Curr.JumpVerletParabola(v0_0, Gate0, t);
        Ring9 p0AtTMinus1 = Phi0Curr.JumpVerletParabola(v0_0, Gate0, t - 1);
        Ring9 p1AtT = Phi1Curr.JumpVerletParabola(v0_1, Gate1, t);
        Ring9 p1AtTMinus1 = Phi1Curr.JumpVerletParabola(v0_1, Gate1, t - 1);
        var c = Create(Mu0, p0AtT, p0AtTMinus1, Mu1, p1AtT, p1AtTMinus1);
        c.Gate0 = Gate0;
        c.Gate1 = Gate1;
        return c;
    }
}

/// <summary>
/// Closed-form Multimodal Sensor Readout from a 72-bit Verlet cell:
/// Optical wavelength [440..718 nm], RGB triad (R+G+B = 128), acoustic pressure P_sound, concentration, and thermal quanta n_E.
/// </summary>
public readonly record struct MultimodalReadout(
    int ConcentrationNorm,
    int WavelengthNm,
    byte ColorR,
    byte ColorG,
    byte ColorB,
    int AcousticPressure,
    int ReactionCurrent,
    int ThermalQuantaNe);

/// <summary>
/// 100% Branchless 9-Bit Ring9 &amp; 12-Neighbor Cuboctahedron O_2(3) Closed-Form Worldline Kernel.
/// Zero `if` statements and zero `?:` ternary operators in the M-layer data path!
/// </summary>
public static class CuboctahedronAutomatonKernel
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int BranchlessAbs(int x) => (x ^ (x >> 31)) - (x >> 31);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static long BranchlessAbs64(long x) => (x ^ (x >> 63)) - (x >> 63);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int BranchlessSignNonNeg64(long x) => 1 + (int)((x >> 63) << 1);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int BranchlessIsNonNeg(int x) => 1 - ((x >> 31) & 1);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int BranchlessIsNonNeg64(long x) => 1 - (int)((x >> 63) & 1);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int BranchlessIsZeroNat(int k) => 1 / (BranchlessMax(k, 0) + 1);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int BranchlessIsPosNat(int k) => 1 - (1 / (BranchlessMax(k, 0) + 1));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int BranchlessBlend(int sel, int a, int b) => sel * a + (1 - sel) * b;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static long BranchlessBlend64(long sel, long a, long b) => sel * a + (1L - sel) * b;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int BranchlessMin(int a, int b) => (int)(((long)a + b - BranchlessAbs64((long)a - b)) >> 1);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int BranchlessMax(int a, int b) => (int)(((long)a + b + BranchlessAbs64((long)a - b)) >> 1);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int FloorDivInt(long num, long den)
    {
        long q = num / den;
        long r = num - q * den;
        return (int)(q - ((r >> 63) & 1L));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int BranchQuotient512(int dphi) =>
        FloorDivInt(dphi + (CarrierSpecD3.PhiPauli - 1), CarrierSpecD3.NRing);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int FaceWindingFromEdges(int k01, int k12, int k32, int k03) =>
        -(k01 + k12 - k32 - k03);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 SelectLeadingGate(Ring9 g0, Ring9 g1)
    {
        int a0 = BranchlessAbs(g0.SignedWrap);
        int a1 = BranchlessAbs(g1.SignedWrap);
        int sel1 = BranchlessIsPosNat(BranchlessMax(a1 - a0, 0));
        return new Ring9(BranchlessBlend(sel1, g1.SignedWrap, g0.SignedWrap));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 GatePhi, Ring9 RawWrapped, int DenomD, int EnergyNe) Compute12NeighborGate(
        Ring9 muC, Ring9 phiC, ReadOnlySpan<Ring9> nbMu, ReadOnlySpan<Ring9> nbPhi, int includeEvenDeficitFlag = 1)
    {
        long zrRaw = 0;
        long ziRaw = 0;
        long sumMu = 0;
        for (int k = 0; k < 12; k++)
        {
            Ring9 dphi = nbPhi[k] - phiC;
            var (cProj, sProj) = dphi.GrayButterflyRotator();
            int m = nbMu[k].Bits9;
            zrRaw += (long)m * cProj;
            ziRaw += (long)m * sProj;
            sumMu += m;
        }
        int mc = muC.Bits9;
        long zr = (long)mc * zrRaw;
        long zi = (long)mc * ziRaw;
        long deficit = (long)mc * sumMu * CarrierSpecD3.QSquare - zr;
        int sgnZi = BranchlessSignNonNeg64(zi);
        long zDrive = zi + includeEvenDeficitFlag * sgnZi * deficit;
        long rawDen = zr + (long)mc * mc * CarrierSpecD3.QSquare + CarrierSpecD3.EpsDenom;
        long isDenPos = BranchlessIsNonNeg64(rawDen - 1);
        long safeDen = BranchlessBlend64(isDenPos, rawDen, CarrierSpecD3.EpsDenom);

        Ring9 rawRing = new(FloorDivInt(CarrierSpecD3.KP * zDrive, safeDen));
        int rawWrap = rawRing.SignedWrap;
        int absRaw = BranchlessAbs(rawWrap);
        int ne = absRaw / CarrierSpecD3.DPhiDisc;
        int gateSigned = BranchlessIsPosNat(ne) * rawWrap;
        return (new Ring9(gateSigned), rawRing, (int)safeDen, ne);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static VerletCell72 StepCellForward(VerletCell72 cell, Ring9 gate0, Ring9 gate1)
    {
        Ring9 p0Next = cell.Phi0Curr + (cell.Phi0Curr - cell.Phi0Past) + gate0;
        Ring9 p1Next = cell.Phi1Curr + (cell.Phi1Curr - cell.Phi1Past) + gate1;
        var next = VerletCell72.Create(cell.Mu0, p0Next, cell.Phi0Curr, cell.Mu1, p1Next, cell.Phi1Curr);
        next.Gate0 = gate0;
        next.Gate1 = gate1;
        return next;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static VerletCell72 StepCellReverse(VerletCell72 steppedCell, Ring9 gate0AtT, Ring9 gate1AtT)
    {
        Ring9 p0Curr = steppedCell.Phi0Past;
        Ring9 p1Curr = steppedCell.Phi1Past;
        Ring9 p0PastRec = p0Curr + (p0Curr - steppedCell.Phi0Curr) + gate0AtT;
        Ring9 p1PastRec = p1Curr + (p1Curr - steppedCell.Phi1Curr) + gate1AtT;
        return VerletCell72.Create(steppedCell.Mu0, p0Curr, p0PastRec, steppedCell.Mu1, p1Curr, p1PastRec);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static MultimodalReadout ReadCellSensors(VerletCell72 cell)
    {
        int mu0 = cell.Mu0.Bits9;
        int mu1 = cell.Mu1.Bits9;
        Ring9 dphi = cell.Phi1Curr - cell.Phi0Curr;
        var (c3, s3) = dphi.GrayButterflyRotator();

        int conc = mu0 * mu0 + mu1 * mu1;
        int wavelengthNm = 579 + FloorDivInt(CarrierSpecD3.CollapsePhi * c3, CarrierSpecD3.QSquare);
        byte r = (byte)BranchlessMax(0, c3);
        byte g = (byte)BranchlessAbs(s3);
        byte b = (byte)BranchlessMax(0, -c3);
        int soundPressure = mu0 * mu1 * (CarrierSpecD3.QSquare - c3);
        int rxnCurrent = mu0 * mu1 * s3;
        int ne = BranchlessAbs(cell.Gate0.SignedWrap) / CarrierSpecD3.DPhiDisc;

        return new MultimodalReadout(conc, wavelengthNm, r, g, b, soundPressure, rxnCurrent, ne);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 Iv263MuAt(Ring9 dz, Ring9 dy, Ring9 dx)
    {
        int z = dz.SignedWrap, y = dy.SignedWrap, x = dx.SignedWrap;
        int r2 = z * z + y * y + x * x;
        int inBall = BranchlessIsNonNeg(CarrierSpecD3.NSq * CarrierSpecD3.NSq - r2);
        return new Ring9(CarrierSpecD3.MuVac + CarrierSpecD3.NLinks * inBall);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 Iv263PhiAt(Ring9 dz, Ring9 dy, Ring9 dx, int lane)
    {
        int z = dz.SignedWrap, y = dy.SignedWrap, x = dx.SignedWrap;
        int r2 = z * z + y * y + x * x;
        int inBall = BranchlessIsNonNeg(CarrierSpecD3.NSq * CarrierSpecD3.NSq - r2);
        int oct = (BranchlessIsNonNeg(z) << 2) + (BranchlessIsNonNeg(y) << 1) + BranchlessIsNonNeg(x);
        int shiftOct = BranchlessBlend(BranchlessIsZeroNat(lane), oct, oct ^ CarrierSpecD3.D);
        int bit = ((CarrierSpecD3.IvSeed % CarrierSpecD3.PhiPauli) >> shiftOct) & 1;
        int rem13 = ((x + CarrierSpecD3.NLinks) % CarrierSpecD3.NPhi + CarrierSpecD3.NPhi) % CarrierSpecD3.NPhi;
        int bg = (rem13 * CarrierSpecD3.DPhiDisc) & CarrierSpecD3.RingMask;
        return new Ring9(bg + (CarrierSpecD3.D * CarrierSpecD3.DPhiDisc) * inBall * bit);
    }

    /// <summary>
    /// Coordinate-Free O(1) Closed-Form Worldline Cell Evaluator at ANY 4D spacetime point `(dz, dy, dx, t)`!
    /// Needs ZERO 3D grid allocation in RAM: evaluates the 12-neighbor holonomy gate at `(dz, dy, dx)` on `Ring9`
    /// and jumps directly along the Verlet parabola to tick `t` in O(1)!
    /// </summary>
    public static VerletCell72 EvaluateClosedFormWorldlineCellAt(Ring9 dz, Ring9 dy, Ring9 dx, int t)
    {
        Span<Ring9> nbMu = stackalloc Ring9[12];
        Span<Ring9> nbP0 = stackalloc Ring9[12];
        Span<Ring9> nbP1 = stackalloc Ring9[12];

        Ring9 muC = Iv263MuAt(dz, dy, dx);
        Ring9 p0_0 = Iv263PhiAt(dz, dy, dx, 0);
        Ring9 p1_0 = Iv263PhiAt(dz, dy, dx, 1);

        for (int k = 0; k < 12; k++)
        {
            var (odz, ody, odx) = CarrierSpecD3.CuboctOffsets12[k];
            Ring9 nz = dz - odz, ny = dy - ody, nx = dx - odx;
            nbMu[k] = Iv263MuAt(nz, ny, nx);
            nbP0[k] = Iv263PhiAt(nz, ny, nx, 0);
            nbP1[k] = Iv263PhiAt(nz, ny, nx, 1);
        }

        Ring9 g0 = Compute12NeighborGate(muC, p0_0, nbMu, nbP0, 1).GatePhi;
        Ring9 g1 = Compute12NeighborGate(muC, p1_0, nbMu, nbP1, 1).GatePhi;
        Ring9 gLead = SelectLeadingGate(g0, g1);

        var cell0 = VerletCell72.Create(muC, p0_0, p0_0, muC, p1_0, p1_0);
        cell0.Gate0 = gLead;
        cell0.Gate1 = gLead;
        return cell0.EvaluateClosedFormAtTick(t);
    }
}

/// <summary>
/// Direct P/Invoke wrapper for NVIDIA CUDA Driver API (`nvcuda.dll`) on RTX 4060 (`sm_89`).
/// </summary>
public static partial class CudaRtx4060Driver
{
    [LibraryImport("nvcuda.dll")]
    private static partial int cuInit(uint flags);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuDeviceGet(out int device, int ordinal);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuDeviceGetName(Span<byte> name, int len, int dev);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuCtxCreate_v2(out IntPtr pctx, uint flags, int dev);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuCtxDestroy_v2(IntPtr ctx);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuModuleLoadDataEx(out IntPtr module, IntPtr image, uint numOptions, IntPtr options, IntPtr optionValues);

    [LibraryImport("nvcuda.dll", StringMarshalling = StringMarshalling.Utf8)]
    private static partial int cuModuleGetFunction(out IntPtr hfunc, IntPtr hmod, string name);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuMemAlloc_v2(out ulong dptr, nuint bytesize);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuMemFree_v2(ulong dptr);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuMemcpyHtoD_v2(ulong dstDevice, IntPtr srcHost, nuint byteCount);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuMemcpyDtoH_v2(IntPtr dstHost, ulong srcDevice, nuint byteCount);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuLaunchKernel(
        IntPtr f, uint gridDimX, uint gridDimY, uint gridDimZ,
        uint blockDimX, uint blockDimY, uint blockDimZ,
        uint sharedMemBytes, IntPtr hStream, IntPtr kernelParams, IntPtr extra);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuCtxSynchronize();

    private const string RotatorVerletPtx = """
        .version 6.5
        .target sm_70
        .address_size 64

        .visible .entry cuboct_rotator_verlet_ptx(
            .param .u64 p_in,
            .param .u64 p_out,
            .param .u32 n_elems
        )
        {
            .reg .pred %p1;
            .reg .b32 %r_tid, %r_ntid, %r_ctaid, %idx, %n, %phi, %w, %tmp, %c_proj;
            .reg .b64 %base_in, %base_out, %addr_in, %addr_out, %off;

            ld.param.u64 %base_in, [p_in];
            ld.param.u64 %base_out, [p_out];
            ld.param.u32 %n, [n_elems];

            mov.u32 %r_tid, %tid.x;
            mov.u32 %r_ntid, %ntid.x;
            mov.u32 %r_ctaid, %ctaid.x;
            mad.lo.u32 %idx, %r_ctaid, %r_ntid, %r_tid;

            setp.ge.u32 %p1, %idx, %n;
            @%p1 bra DONE;

            mul.wide.u32 %off, %idx, 4;
            add.u64 %addr_in, %base_in, %off;
            add.u64 %addr_out, %base_out, %off;

            ld.global.u32 %phi, [%addr_in];
            and.b32 %phi, %phi, 511;

            mov.s32 %tmp, 256;
            sub.s32 %tmp, %tmp, %phi;
            and.b32 %tmp, %tmp, 511;
            mov.s32 %w, 256;
            sub.s32 %w, %w, %tmp;

            abs.s32 %tmp, %w;
            mov.s32 %c_proj, 128;
            sub.s32 %c_proj, %c_proj, %tmp;

            st.global.s32 [%addr_out], %c_proj;
        DONE:
            ret;
        }
        """;

    public static unsafe (bool Ok, string GpuName, double GbpSec, int VerifiedCount) RunRtx4060PtxSelfTest()
    {
        if (cuInit(0) != 0) return (false, "cuInit failed", 0, 0);
        if (cuDeviceGet(out int dev, 0) != 0) return (false, "cuDeviceGet failed", 0, 0);

        Span<byte> nameBuf = stackalloc byte[128];
        cuDeviceGetName(nameBuf, 128, dev);
        int nullIdx = nameBuf.IndexOf((byte)0);
        string gpuName = System.Text.Encoding.UTF8.GetString(nameBuf[..(nullIdx > 0 ? nullIdx : 128)]);

        if (cuCtxCreate_v2(out IntPtr ctx, 0, dev) != 0) return (false, gpuName, 0, 0);
        try
        {
            byte[] ptxBytes = System.Text.Encoding.ASCII.GetBytes(RotatorVerletPtx.Replace("\r\n", "\n") + "\0");
            byte[] errLog = new byte[4096];
            int loadRc;
            IntPtr mod;
            fixed (byte* pPtx = ptxBytes)
            fixed (byte* pErr = errLog)
            {
                int* opts = stackalloc int[2] { 5, 6 };
                void** vals = stackalloc void*[2] { pErr, (void*)4096 };
                loadRc = cuModuleLoadDataEx(out mod, (IntPtr)pPtx, 2, (IntPtr)opts, (IntPtr)vals);
            }
            if (loadRc != 0)
            {
                int errEnd = Array.IndexOf(errLog, (byte)0);
                string jitMsg = System.Text.Encoding.UTF8.GetString(errLog, 0, errEnd > 0 ? errEnd : 256).Trim();
                return (false, $"{gpuName} (PTX err={loadRc}: {jitMsg})", 0, 0);
            }
            if (cuModuleGetFunction(out IntPtr fn, mod, "cuboct_rotator_verlet_ptx") != 0) return (false, gpuName, 0, 0);

            const int N = 16 * 1024 * 1024;
            int[] hostIn = new int[N];
            int[] hostOut = new int[N];
            for (int i = 0; i < N; i++) hostIn[i] = i & 511;

            nuint byteSize = (nuint)(N * sizeof(int));
            cuMemAlloc_v2(out ulong dIn, byteSize);
            cuMemAlloc_v2(out ulong dOut, byteSize);
            try
            {
                fixed (int* pIn = hostIn)
                    cuMemcpyHtoD_v2(dIn, (IntPtr)pIn, byteSize);

                uint nElems = (uint)N;
                ulong argIn = dIn;
                ulong argOut = dOut;
                void** kParams = stackalloc void*[3];
                kParams[0] = &argIn;
                kParams[1] = &argOut;
                kParams[2] = &nElems;

                cuLaunchKernel(fn, (uint)(N / 256), 1, 1, 256, 1, 1, 0, IntPtr.Zero, (IntPtr)kParams, IntPtr.Zero);
                cuCtxSynchronize();

                var sw = Stopwatch.StartNew();
                const int Iters = 20;
                for (int it = 0; it < Iters; it++)
                {
                    cuLaunchKernel(fn, (uint)(N / 256), 1, 1, 256, 1, 1, 0, IntPtr.Zero, (IntPtr)kParams, IntPtr.Zero);
                }
                cuCtxSynchronize();
                sw.Stop();

                fixed (int* pOut = hostOut)
                    cuMemcpyDtoH_v2((IntPtr)pOut, dOut, byteSize);

                int verified = 0;
                for (int i = 0; i < 512; i++)
                {
                    if (hostOut[i] == new Ring9(i).ProjC) verified++;
                }

                double totalBytes = (double)N * sizeof(int) * 2.0 * Iters;
                double gbpSec = (totalBytes / (1024.0 * 1024.0 * 1024.0)) / sw.Elapsed.TotalSeconds;
                return (verified == 512, gpuName, gbpSec, verified);
            }
            finally
            {
                cuMemFree_v2(dIn);
                cuMemFree_v2(dOut);
            }
        }
        finally
        {
            cuCtxDestroy_v2(ctx);
        }
    }
}
