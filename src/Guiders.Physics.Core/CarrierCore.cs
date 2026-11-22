using System;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace Guiders.Physics.Core;

public readonly record struct Ring9
{
    private readonly ushort _bits9;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public Ring9(ushort bits9) => _bits9 = (ushort)(bits9 & 0x1FF);

    public ushort Bits9
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _bits9;
    }

    public static Ring9 Zero           => new(0);
    public static Ring9 MuVac          => new(1);
    public static Ring9 C0OverCSq      => new(2);
    public static Ring9 D              => new(3);
    public static Ring9 DInSquare      => new(4);
    public static Ring9 DOutEdge       => new(5);
    public static Ring9 NSq            => new(6);
    public static Ring9 DOutVertex     => new(7);
    public static Ring9 NTri           => new(8);
    public static Ring9 ModBits        => new(9);
    public static Ring9 NLinks         => new(12);
    public static Ring9 NPhi           => new(13);
    public static Ring9 ZSilicon14     => new(14);
    public static Ring9 ZOxygen2_16    => new(16);
    public static Ring9 SeamTicks      => new(21);
    public static Ring9 BezoutInv41    => new(25);
    public static Ring9 ZQuartz30      => new(30);
    public static Ring9 DPhiDisc       => new(41);
    public static Ring9 KP             => new(64);
    public static Ring9 ZAlphaQuartz90 => new(90);
    public static Ring9 MHat           => new(97);
    public static Ring9 QSquare        => new(128);
    public static Ring9 CollapsePhi    => new(139);
    public static Ring9 BezoutInv139   => new(183);
    public static Ring9 BezoutInv13    => new(197);
    public static Ring9 QuartzSiO2Phi  => new(206);
    public static Ring9 PhiPauli       => new(256);
    public static Ring9 IvSeed         => new(263);
    public static Ring9 BezoutInv97    => new(417);
    public static Ring9 MinusOne       => new(511);

    public static readonly (Ring9 Dz, Ring9 Dy, Ring9 Dx)[] CuboctOffsets12 =
    [
        (Zero,     MuVac,    MuVac),    (Zero,     MuVac,    MinusOne),
        (Zero,     MinusOne, MuVac),    (Zero,     MinusOne, MinusOne),
        (MuVac,    Zero,     MuVac),    (MuVac,    Zero,     MinusOne),
        (MinusOne, Zero,     MuVac),    (MinusOne, Zero,     MinusOne),
        (MuVac,    MuVac,    Zero),     (MuVac,    MinusOne, Zero),
        (MinusOne, MuVac,    Zero),     (MinusOne, MinusOne, Zero)
    ];

    public Ring9 BitQ1
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new((ushort)((_bits9 >> 8) & 1));
    }

    public Ring9 BitQ0
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new((ushort)((_bits9 >> 7) & 1));
    }

    public Ring9 Mantissa7
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new((ushort)(_bits9 & 127));
    }

    public Ring9 PauliPartner
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new((ushort)(_bits9 ^ 256));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public (Ring9 ProjC, Ring9 ProjS) GrayButterflyRotator() => (ProjC, ProjS);

    public Ring9 ProjC
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get
        {
            ushort q1 = (ushort)((_bits9 >> 8) & 1);
            ushort q0 = (ushort)((_bits9 >> 7) & 1);
            ushort m7 = (ushort)(_bits9 & 127);
            ushort comp7 = (ushort)(128 - m7);
            ushort magC = (ushort)((1 - q0) * comp7 + q0 * m7);
            ushort cNeg = (ushort)(q1 ^ q0);
            return new Ring9((ushort)((1 - 2 * cNeg) * magC));
        }
    }

    public Ring9 ProjS
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get
        {
            ushort q1 = (ushort)((_bits9 >> 8) & 1);
            ushort q0 = (ushort)((_bits9 >> 7) & 1);
            ushort m7 = (ushort)(_bits9 & 127);
            ushort comp7 = (ushort)(128 - m7);
            ushort magS = (ushort)((1 - q0) * m7 + q0 * comp7);
            return new Ring9((ushort)((1 - 2 * q1) * magS));
        }
    }

    public Ring9 AbsWrap
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get
        {
            ushort q1 = (ushort)((_bits9 >> 8) & 1);
            ushort q0 = (ushort)((_bits9 >> 7) & 1);
            ushort m7 = (ushort)(_bits9 & 127);
            ushort comp7 = (ushort)(128 - m7);
            ushort magC = (ushort)((1 - q0) * comp7 + q0 * m7);
            ushort cNeg = (ushort)(q1 ^ q0);
            ushort absW = (ushort)(128 - (1 - 2 * cNeg) * magC);
            return new Ring9(absW);
        }
    }

    public Ring9 DecodeBezout41
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => this * BezoutInv41;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public Ring9 JumpVerletParabola(Ring9 v0, Ring9 gatePhi, Ring9 t)
    {
        ushort t10 = t._bits9;
        ushort tri9 = (ushort)(((t10 * (t10 + 1)) >> 1) & 0x1FF);
        Ring9 triT = new(tri9);
        return this + (t * v0) + (triT * gatePhi);
    }

    public (Ring9 O2, Ring9 O1, Ring9 O0) OctalTriad =>
        (new((ushort)((_bits9 >> 6) & 7)), new((ushort)((_bits9 >> 3) & 7)), new((ushort)(_bits9 & 7)));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 operator +(Ring9 a, Ring9 b) => new((ushort)(a._bits9 + b._bits9));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 operator -(Ring9 a, Ring9 b) => new((ushort)(a._bits9 - b._bits9));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 operator -(Ring9 a) => new((ushort)(512 - a._bits9));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 operator *(Ring9 a, Ring9 b) => new((ushort)(a._bits9 * b._bits9));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 operator ^(Ring9 a, Ring9 b) => new((ushort)(a._bits9 ^ b._bits9));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 operator &(Ring9 a, Ring9 b) => new((ushort)(a._bits9 & b._bits9));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 operator |(Ring9 a, Ring9 b) => new((ushort)(a._bits9 | b._bits9));

    public override string ToString()
    {
        var (o2, o1, o0) = OctalTriad;
        return $"{_bits9} (0o{o2._bits9}{o1._bits9}{o0._bits9})";
    }
}

/// <summary>
/// Hardware 9-Bit Full Adder &amp; Multiplier Object (`Adder9`):
/// Splits any 9-bit addition or multiplication into:
///   - `Sum9` (low 9 bits [8:0] = Channel 1 wave phase on Z_512),
///   - `Carry9` (high carry-out bits [17:9] = Channel 2 corpuscle quantum!).
/// </summary>
public readonly record struct Adder9
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 Sum9, Ring9 Carry9) AddWithCarry(Ring9 a, Ring9 b, Ring9 carryIn = default)
    {
        ushort full = (ushort)(a.Bits9 + b.Bits9 + carryIn.Bits9);
        return (new Ring9(full), new Ring9((ushort)(full >> 9)));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 Lo9, Ring9 Hi9) MultiplyFull(Ring9 a, Ring9 b)
    {
        uint prod = (uint)a.Bits9 * b.Bits9;
        return (new Ring9((ushort)prod), new Ring9((ushort)(prod >> 9)));
    }
}

/// <summary>
/// 4-Lane Cascaded 9-Bit Bus (`Bus36`: `L0`, `L1`, `L2`, `L3` = 4 x 9 bits = 36 bits = `B_1(3)`!).
/// Represents a 36-bit two's-complement wire bundle built purely from four `Ring9` lanes
/// chained through `Adder9.AddWithCarry`!
/// </summary>
public readonly record struct Bus36(Ring9 L0, Ring9 L1, Ring9 L2, Ring9 L3)
{
    public static Bus36 Zero     => new(Ring9.Zero, Ring9.Zero, Ring9.Zero, Ring9.Zero);
    public static Bus36 One      => new(Ring9.MuVac, Ring9.Zero, Ring9.Zero, Ring9.Zero);
    public static Bus36 MinusOne => new(Ring9.MinusOne, Ring9.MinusOne, Ring9.MinusOne, Ring9.MinusOne);
    public static Bus36 QSquare  => new(Ring9.QSquare, Ring9.Zero, Ring9.Zero, Ring9.Zero);
    public static Bus36 KP       => new(Ring9.KP, Ring9.Zero, Ring9.Zero, Ring9.Zero);
    public static Bus36 EpsDenom => new(Ring9.Zero, Ring9.ZOxygen2_16, Ring9.Zero, Ring9.Zero); // 16 * 512 = 8192!

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Bus36 FromUnsignedRing9(Ring9 u) =>
        new(u, Ring9.Zero, Ring9.Zero, Ring9.Zero);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Bus36 FromSignedRing9(Ring9 s)
    {
        ushort isNeg = s.BitQ1.Bits9;
        Ring9 ext = new((ushort)(511 * isNeg));
        return new(s, ext, ext, ext);
    }

    public Ring9 SignBit
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => L3.BitQ1;
    }

    public Ring9 IsNonNegative
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => Ring9.MuVac - L3.BitQ1;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Bus36 operator +(Bus36 a, Bus36 b)
    {
        var (s0, c0) = Adder9.AddWithCarry(a.L0, b.L0);
        var (s1, c1) = Adder9.AddWithCarry(a.L1, b.L1, c0);
        var (s2, c2) = Adder9.AddWithCarry(a.L2, b.L2, c1);
        var (s3, _)  = Adder9.AddWithCarry(a.L3, b.L3, c2);
        return new Bus36(s0, s1, s2, s3);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Bus36 operator ~(Bus36 a) =>
        new(a.L0 ^ Ring9.MinusOne, a.L1 ^ Ring9.MinusOne, a.L2 ^ Ring9.MinusOne, a.L3 ^ Ring9.MinusOne);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Bus36 operator -(Bus36 a) => (~a) + One;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Bus36 operator -(Bus36 a, Bus36 b) => a + (-b);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Bus36 operator *(Bus36 a, Ring9 k)
    {
        var (p0Lo, p0Hi) = Adder9.MultiplyFull(a.L0, k);
        var (p1Lo, p1Hi) = Adder9.MultiplyFull(a.L1, k);
        var (p2Lo, p2Hi) = Adder9.MultiplyFull(a.L2, k);
        var (p3Lo, _)    = Adder9.MultiplyFull(a.L3, k);

        var (s1, c1) = Adder9.AddWithCarry(p1Lo, p0Hi);
        var (s2, c2) = Adder9.AddWithCarry(p2Lo, p1Hi, c1);
        var (s3, _)  = Adder9.AddWithCarry(p3Lo, p2Hi, c2);
        return new Bus36(p0Lo, s1, s2, s3);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Bus36 operator *(Bus36 a, Bus36 b)
    {
        ulong ua = (ulong)a.L0.Bits9 | ((ulong)a.L1.Bits9 << 9) | ((ulong)a.L2.Bits9 << 18) | ((ulong)a.L3.Bits9 << 27);
        ulong ub = (ulong)b.L0.Bits9 | ((ulong)b.L1.Bits9 << 9) | ((ulong)b.L2.Bits9 << 18) | ((ulong)b.L3.Bits9 << 27);
        ulong prod = (ua * ub) & 0xF_FFFF_FFFFUL;
        return new Bus36(
            new Ring9((ushort)prod),
            new Ring9((ushort)(prod >> 9)),
            new Ring9((ushort)(prod >> 18)),
            new Ring9((ushort)(prod >> 27)));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 FloorDivideToRing9(Bus36 num, Bus36 den)
    {
        ulong uNum = (ulong)num.L0.Bits9 | ((ulong)num.L1.Bits9 << 9) | ((ulong)num.L2.Bits9 << 18) | ((ulong)num.L3.Bits9 << 27);
        ulong uDen = (ulong)den.L0.Bits9 | ((ulong)den.L1.Bits9 << 9) | ((ulong)den.L2.Bits9 << 18) | ((ulong)den.L3.Bits9 << 27);

        ulong numNeg = (uNum >> 35) & 1UL;
        ulong magNum = ((1UL - numNeg) * uNum + numNeg * ((0x10_0000_0000UL - uNum) & 0xF_FFFF_FFFFUL)) & 0xF_FFFF_FFFFUL;
        ulong qMag = magNum / uDen;
        ulong rMag = magNum - qMag * uDen;
        ulong hasRem = 1UL - (1UL / (rMag + 1UL));
        ulong qFloor = (1UL - numNeg) * qMag + numNeg * (0x10_0000_0000UL - qMag - hasRem);
        return new Ring9((ushort)(qFloor & 0x1FFUL));
    }
}

[StructLayout(LayoutKind.Sequential, Pack = 1, Size = 9)]
public struct VerletCell72
{
    private ulong _lo64;
    private byte _hi8;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public readonly Ring9 GetLane9(Ring9 laneIdx)
    {
        UInt128 w = (UInt128)_lo64 | ((UInt128)_hi8 << 64);
        ushort bits = (ushort)((w >> ((laneIdx.Bits9 & 7) * 9)) & 0x1FFUL);
        return new Ring9(bits);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void SetLane9(Ring9 laneIdx, Ring9 val)
    {
        ushort shift = (ushort)((laneIdx.Bits9 & 7) * 9);
        UInt128 w = (UInt128)_lo64 | ((UInt128)_hi8 << 64);
        UInt128 mask = ~((UInt128)0x1FFUL << shift);
        w = (w & mask) | ((UInt128)val.Bits9 << shift);
        _lo64 = (ulong)w;
        _hi8 = (byte)(w >> 64);
    }

    public Ring9 Mu0      { readonly get => GetLane9(new(0)); set => SetLane9(new(0), value); }
    public Ring9 Phi0Curr { readonly get => GetLane9(new(1)); set => SetLane9(new(1), value); }
    public Ring9 Phi0Past { readonly get => GetLane9(new(2)); set => SetLane9(new(2), value); }
    public Ring9 Mu1      { readonly get => GetLane9(new(3)); set => SetLane9(new(3), value); }
    public Ring9 Phi1Curr { readonly get => GetLane9(new(4)); set => SetLane9(new(4), value); }
    public Ring9 Phi1Past { readonly get => GetLane9(new(5)); set => SetLane9(new(5), value); }
    public Ring9 Gate0    { readonly get => GetLane9(new(6)); set => SetLane9(new(6), value); }
    public Ring9 Gate1    { readonly get => GetLane9(new(7)); set => SetLane9(new(7), value); }

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
    public static VerletCell72 CreatePhysicalAlphaQuartzCell(Ring9 crystalOscillatorA, Ring9 crystalOscillatorB)
    {
        Ring9 phiSi = Ring9.QuartzSiO2Phi + crystalOscillatorA;
        Ring9 phiO2 = Ring9.QuartzSiO2Phi + crystalOscillatorB;
        Ring9 dphiBeat = phiO2 - phiSi;
        var cell = Create(Ring9.ZSilicon14, phiSi, phiSi - dphiBeat, Ring9.ZOxygen2_16, phiO2, phiO2 - dphiBeat);
        cell.Gate0 = CuboctahedronAutomatonKernel.ApplyHeisenbergDeadband(dphiBeat);
        cell.Gate1 = cell.Gate0;
        return cell;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public readonly VerletCell72 EvaluateClosedFormAtTick(Ring9 t)
    {
        Ring9 v0_0 = Phi0Curr - Phi0Past;
        Ring9 v0_1 = Phi1Curr - Phi1Past;
        Ring9 tPrev = t - Ring9.MuVac;
        Ring9 p0AtT = Phi0Curr.JumpVerletParabola(v0_0, Gate0, t);
        Ring9 p0AtTMinus1 = Phi0Curr.JumpVerletParabola(v0_0, Gate0, tPrev);
        Ring9 p1AtT = Phi1Curr.JumpVerletParabola(v0_1, Gate1, t);
        Ring9 p1AtTMinus1 = Phi1Curr.JumpVerletParabola(v0_1, Gate1, tPrev);
        var c = Create(Mu0, p0AtT, p0AtTMinus1, Mu1, p1AtT, p1AtTMinus1);
        c.Gate0 = Gate0;
        c.Gate1 = Gate1;
        return c;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public readonly VerletCell72 StepForward(Ring9 gate) =>
        CuboctahedronAutomatonKernel.StepCellForward(this, gate, gate);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public readonly VerletCell72 StepReverse(Ring9 gateAtT) =>
        CuboctahedronAutomatonKernel.StepCellReverse(this, gateAtT, gateAtT);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public readonly MultimodalReadout9 ReadSensors() =>
        CuboctahedronAutomatonKernel.ReadCellSensors(this);
}

public sealed class CuboctahedronCarrier
{
    private VerletCell72[,,] _grid;
    private VerletCell72[,,] _scratch;
    private readonly ushort _n;

    public Ring9 Size => new(_n);
    public Ring9 Tick { get; private set; }

    public CuboctahedronCarrier(Ring9 size)
    {
        _n = size.Bits9;
        _grid = new VerletCell72[_n, _n, _n];
        _scratch = new VerletCell72[_n, _n, _n];
        Tick = Ring9.Zero;
        ushort c = (ushort)(_n >> 1);
        Ring9 center = new(c);
        for (ushort iz = 0; iz < _n; iz++)
            for (ushort iy = 0; iy < _n; iy++)
                for (ushort ix = 0; ix < _n; ix++)
                {
                    Ring9 dz = new Ring9(iz) - center;
                    Ring9 dy = new Ring9(iy) - center;
                    Ring9 dx = new Ring9(ix) - center;
                    Ring9 mu = CuboctahedronAutomatonKernel.Iv263MuAt(dz, dy, dx);
                    Ring9 p0 = CuboctahedronAutomatonKernel.Iv263PhiAt(dz, dy, dx, Ring9.Zero, center);
                    Ring9 p1 = CuboctahedronAutomatonKernel.Iv263PhiAt(dz, dy, dx, Ring9.MuVac, center);
                    _grid[iz, iy, ix] = VerletCell72.Create(mu, p0, p0, mu, p1, p1);
                }
    }

    public VerletCell72 this[Ring9 z, Ring9 y, Ring9 x]
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _grid[(z.Bits9 + _n) % _n, (y.Bits9 + _n) % _n, (x.Bits9 + _n) % _n];
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set => _grid[(z.Bits9 + _n) % _n, (y.Bits9 + _n) % _n, (x.Bits9 + _n) % _n] = value;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public MultimodalReadout9 ReadSensors(Ring9 z, Ring9 y, Ring9 x) =>
        this[z, y, x].ReadSensors();

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private VerletCell72 CellAt(ushort z, ushort y, ushort x) =>
        _grid[z % _n, y % _n, x % _n];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private Ring9 ComputeGateAtNode(ushort z, ushort y, ushort x, Ring9 usePastFlag)
    {
        ushort n = _n;
        ushort zm = (ushort)(z + n - 1), zp = (ushort)(z + 1);
        ushort ym = (ushort)(y + n - 1), yp = (ushort)(y + 1);
        ushort xm = (ushort)(x + n - 1), xp = (ushort)(x + 1);

        var c   = CellAt(z,  y,  x);
        var n0  = CellAt(z,  ym, xm);
        var n1  = CellAt(z,  ym, xp);
        var n2  = CellAt(z,  yp, xm);
        var n3  = CellAt(z,  yp, xp);
        var n4  = CellAt(zm, y,  xm);
        var n5  = CellAt(zm, y,  xp);
        var n6  = CellAt(zp, y,  xm);
        var n7  = CellAt(zp, y,  xp);
        var n8  = CellAt(zm, ym, x);
        var n9  = CellAt(zm, yp, x);
        var n10 = CellAt(zp, ym, x);
        var n11 = CellAt(zp, yp, x);

        Ring9 selPast = usePastFlag;
        Ring9 selCurr = Ring9.MuVac - selPast;

        Ring9 P0(VerletCell72 v) => selCurr * v.Phi0Curr + selPast * v.Phi0Past;
        Ring9 P1(VerletCell72 v) => selCurr * v.Phi1Curr + selPast * v.Phi1Past;

        return CuboctahedronAutomatonKernel.ComputeTwoLaneCellGate(
            c.Mu0, P0(c), c.Mu1, P1(c),
            n0.Mu0,  P0(n0),  n1.Mu0,  P0(n1),  n2.Mu0,  P0(n2),  n3.Mu0,  P0(n3),
            n4.Mu0,  P0(n4),  n5.Mu0,  P0(n5),  n6.Mu0,  P0(n6),  n7.Mu0,  P0(n7),
            n8.Mu0,  P0(n8),  n9.Mu0,  P0(n9),  n10.Mu0, P0(n10), n11.Mu0, P0(n11),
            n0.Mu1,  P1(n0),  n1.Mu1,  P1(n1),  n2.Mu1,  P1(n2),  n3.Mu1,  P1(n3),
            n4.Mu1,  P1(n4),  n5.Mu1,  P1(n5),  n6.Mu1,  P1(n6),  n7.Mu1,  P1(n7),
            n8.Mu1,  P1(n8),  n9.Mu1,  P1(n9),  n10.Mu1, P1(n10), n11.Mu1, P1(n11),
            Ring9.MuVac);
    }

    public void StepForward()
    {
        ushort n = _n;
        for (ushort z = 0; z < n; z++)
            for (ushort y = 0; y < n; y++)
                for (ushort x = 0; x < n; x++)
                {
                    Ring9 gate = ComputeGateAtNode(z, y, x, Ring9.Zero);
                    _scratch[z, y, x] = _grid[z, y, x].StepForward(gate);
                }
        (_grid, _scratch) = (_scratch, _grid);
        Tick = Tick + Ring9.MuVac;
    }

    public void StepBackward()
    {
        ushort n = _n;
        for (ushort z = 0; z < n; z++)
            for (ushort y = 0; y < n; y++)
                for (ushort x = 0; x < n; x++)
                {
                    Ring9 gateAtT = ComputeGateAtNode(z, y, x, Ring9.MuVac);
                    _scratch[z, y, x] = _grid[z, y, x].StepReverse(gateAtT);
                }
        (_grid, _scratch) = (_scratch, _grid);
        Tick = Tick - Ring9.MuVac;
    }
}

public readonly record struct MultimodalReadout9(
    Ring9 ColorR,
    Ring9 ColorG,
    Ring9 ColorB,
    Ring9 AcousticWave,
    Ring9 ReactionCurrent,
    Ring9 EnergyQuantaNe,
    Ring9 OlfactoryShellZ,
    Ring9 ChiralOdorProjectionS);

public static unsafe class CuboctahedronAutomatonKernel
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
    public static Ring9 HeisenbergEnergyQuanta(Ring9 gateRaw)
    {
        ushort abs9 = gateRaw.AbsWrap.Bits9;
        return new Ring9((ushort)(abs9 / 41));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 ApplyHeisenbergDeadband(Ring9 gateRaw)
    {
        ushort ne = (ushort)(gateRaw.AbsWrap.Bits9 / 41);
        ushort pass = (ushort)(1 - (1 / (ne + 1)));
        return new Ring9((ushort)(pass * gateRaw.Bits9));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 AtomicTimeServerPfdError(Ring9 phiAtom, Ring9 phiMachine) =>
        (phiAtom - phiMachine).ProjS;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 AtomicSelfClockedDrivePhase(Ring9 phiAtom, Ring9 deltaPhiMode) =>
        phiAtom + deltaPhiMode;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 AtomicRelativePhaseAtTarget(Ring9 phiAtom, Ring9 deltaPhiMode) =>
        AtomicSelfClockedDrivePhase(phiAtom, deltaPhiMode) - phiAtom;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 DecodeAtomicSpatialCoord(Ring9 phiAtomBg) =>
        phiAtomBg.DecodeBezout41;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 DecodeAtomicTimeTick(Ring9 deltaPhiVerlet) =>
        deltaPhiVerlet.DecodeBezout41;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 TwoPointWorldlineSlopeRecovery(Ring9 phiOrigin0, Ring9 phiNowAt1) =>
        phiNowAt1 - phiOrigin0;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 TwoPointWorldlineInterpolate(Ring9 phiOrigin0, Ring9 phiNowAt1, Ring9 phiImp, Ring9 tau)
    {
        Ring9 v0 = TwoPointWorldlineSlopeRecovery(phiOrigin0, phiNowAt1) - phiImp;
        return phiOrigin0.JumpVerletParabola(v0, phiImp, tau);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 ChronovisorArbitraryTargetPhase(
        Ring9 phiAnchorNow, Ring9 deltaPhiAnchorNow, Ring9 xTarget, Ring9 tauTarget)
    {
        Ring9 xNow = DecodeAtomicSpatialCoord(phiAnchorNow);
        Ring9 tauNow = DecodeAtomicTimeTick(deltaPhiAnchorNow);
        Ring9 deltaX = xTarget - xNow;
        Ring9 deltaTau = tauTarget - tauNow;
        Ring9 phiAtTargetX = phiAnchorNow + (deltaX * Ring9.DPhiDisc);
        return phiAtTargetX.JumpVerletParabola(deltaPhiAnchorNow, Ring9.DPhiDisc, deltaTau);
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
    public static MultimodalReadout9 ReadCellSensors(VerletCell72 cell)
    {
        Ring9 phi0 = cell.Phi0Curr;
        Ring9 q1 = phi0.BitQ1;
        Ring9 q0 = phi0.BitQ0;
        Ring9 m7 = phi0.Mantissa7;
        Ring9 comp7 = Ring9.QSquare - m7;

        Ring9 magC = (Ring9.MuVac - q0) * comp7 + q0 * m7;
        Ring9 magS = (Ring9.MuVac - q0) * m7 + q0 * comp7;
        Ring9 cNeg = q1 ^ q0;

        Ring9 r = (Ring9.MuVac - cNeg) * magC;
        Ring9 g = magS;
        Ring9 b = cNeg * magC;

        Ring9 c3 = phi0.ProjC;
        Ring9 s3 = phi0.ProjS;
        Ring9 sound = cell.Mu0 * cell.Mu1 * (Ring9.QSquare - c3);
        Ring9 rxn = cell.Mu0 * cell.Mu1 * s3;
        Ring9 ne = HeisenbergEnergyQuanta(cell.Gate0);
        Ring9 odorShellZ = cell.Mu0 + cell.Mu1;

        return new MultimodalReadout9(r, g, b, sound, rxn, ne, odorShellZ, s3);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 Iv263InBall9(Ring9 dz, Ring9 dy, Ring9 dx)
    {
        ushort az = dz.AbsWrap.Bits9;
        ushort ay = dy.AbsWrap.Bits9;
        ushort ax = dx.AbsWrap.Bits9;
        ushort r2 = (ushort)(az * az + ay * ay + ax * ax);
        ushort inBall = (ushort)(1 / ((r2 / 37) + 1));
        return new Ring9(inBall);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 Iv263MuAt(Ring9 dz, Ring9 dy, Ring9 dx) =>
        Ring9.MuVac + Ring9.NLinks * Iv263InBall9(dz, dy, dx);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 Iv263PhiAt(Ring9 dz, Ring9 dy, Ring9 dx, Ring9 lane, Ring9 centerX = default)
    {
        Ring9 inBall = Iv263InBall9(dz, dy, dx);
        ushort zPos = (ushort)(1 - dz.BitQ1.Bits9);
        ushort yPos = (ushort)(1 - dy.BitQ1.Bits9);
        ushort xPos = (ushort)(1 - dx.BitQ1.Bits9);
        ushort oct = (ushort)((zPos << 2) | (yPos << 1) | xPos);
        ushort isLane1 = (ushort)(lane.Bits9 & 1);
        ushort shift = (ushort)(oct ^ (3 * isLane1));
        ushort bit = (ushort)(((263 % 256) >> shift) & 1);

        ushort cx = (ushort)((1 - (1 / (centerX.Bits9 + 1))) * centerX.Bits9 + (1 / (centerX.Bits9 + 1)) * 12);
        ushort xShifted = (ushort)(dx.Bits9 + cx + 260 - 512 * dx.BitQ1.Bits9);
        ushort rem13 = (ushort)(xShifted % 13);
        ushort bg = (ushort)((rem13 * 41) & 0x1FF);
        return new Ring9((ushort)(bg + 123 * inBall.Bits9 * bit));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Bus36 NeighborTerm(Ring9 m, Ring9 proj)
    {
        ushort isNeg = proj.BitQ1.Bits9;
        Ring9 ext = new((ushort)(511 * isNeg));
        Bus36 sProj = new(proj, ext, ext, ext);
        return sProj * m;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 Isqrt4096(ushort q)
    {
        ushort r = 0;
        ushort c6 = (ushort)(r | 64); ushort b6 = (ushort)(1 / (((c6 * c6) / (q + 1)) + 1)); r = (ushort)(r | (b6 * 64));
        ushort c5 = (ushort)(r | 32); ushort b5 = (ushort)(1 / (((c5 * c5) / (q + 1)) + 1)); r = (ushort)(r | (b5 * 32));
        ushort c4 = (ushort)(r | 16); ushort b4 = (ushort)(1 / (((c4 * c4) / (q + 1)) + 1)); r = (ushort)(r | (b4 * 16));
        ushort c3 = (ushort)(r |  8); ushort b3 = (ushort)(1 / (((c3 * c3) / (q + 1)) + 1)); r = (ushort)(r | (b3 *  8));
        ushort c2 = (ushort)(r |  4); ushort b2 = (ushort)(1 / (((c2 * c2) / (q + 1)) + 1)); r = (ushort)(r | (b2 *  4));
        ushort c1 = (ushort)(r |  2); ushort b1 = (ushort)(1 / (((c1 * c1) / (q + 1)) + 1)); r = (ushort)(r | (b1 *  2));
        ushort c0 = (ushort)(r |  1); ushort b0 = (ushort)(1 / (((c0 * c0) / (q + 1)) + 1)); r = (ushort)(r | (b0 *  1));
        return new Ring9(r);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Bus36 Px, Bus36 Py, Bus36 Pr) PythagoreanRotator36(Ring9 phi)
    {
        Ring9 c = phi.ProjC;
        Ring9 s = phi.ProjS;
        Ring9 absC = (Ring9.QSquare - phi.AbsWrap).AbsWrap;
        Ring9 absS = (Ring9.QSquare - (phi - Ring9.QSquare).AbsWrap).AbsWrap;

        ushort sMag = absS.Bits9;
        ushort cMag = absC.Bits9;
        ushort cNonZero = (ushort)(1 - (1 / (cMag + 1)));
        ushort cNeg = c.BitQ1.Bits9;
        ushort sNeg = s.BitQ1.Bits9;

        uint q2 = 128u * 128u;
        uint s2 = (uint)sMag * sMag;
        uint pxMag = cNonZero * (q2 - s2);
        uint pyMag = 256u * sMag;
        uint prVal = q2 + s2;

        Bus36 pxPos = new(new Ring9((ushort)pxMag), new Ring9((ushort)(pxMag >> 9)), Ring9.Zero, Ring9.Zero);
        Bus36 px = (pxPos * new Ring9((ushort)(1 - cNeg))) + ((-pxPos) * new Ring9(cNeg));

        Bus36 pyPos = new(new Ring9((ushort)pyMag), new Ring9((ushort)(pyMag >> 9)), Ring9.Zero, Ring9.Zero);
        Bus36 py = (pyPos * new Ring9((ushort)(1 - sNeg))) + ((-pyPos) * new Ring9(sNeg));

        Bus36 pr = new(new Ring9((ushort)prVal), new Ring9((ushort)(prVal >> 9)), Ring9.Zero, Ring9.Zero);
        return (px, py, pr);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 EnforcePlanckCellFloor(Ring9 mu)
    {
        ushort isZero = (ushort)(1 / (mu.Bits9 + 1));
        return new Ring9((ushort)(mu.Bits9 + isZero));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Bus36 ComputeSingleDirectedEdgeFlux36(Ring9 muX, Ring9 phiX, Ring9 muY, Ring9 phiY)
    {
        Ring9 dphiXY = phiX - phiY;
        Bus36 sVal = NeighborTerm(muY, dphiXY.ProjS) * muX;
        Bus36 denom1536 = new(Ring9.Zero, Ring9.D, Ring9.Zero, Ring9.Zero); // 3 * 512 = 1536 = 12 * 128!

        ulong uNum = (ulong)sVal.L0.Bits9 | ((ulong)sVal.L1.Bits9 << 9) | ((ulong)sVal.L2.Bits9 << 18) | ((ulong)sVal.L3.Bits9 << 27);
        ulong numNeg = (uNum >> 35) & 1UL;
        ulong magNum = ((1UL - numNeg) * uNum + numNeg * ((0x10_0000_0000UL - uNum) & 0xF_FFFF_FFFFUL)) & 0xF_FFFF_FFFFUL;
        ulong qMag = magNum / 1536UL;
        ulong rMag = magNum - qMag * 1536UL;
        ulong hasRem = 1UL - (1UL / (rMag + 1UL));

        ushort qX = (ushort)(muX.Bits9 * muX.Bits9);
        ushort qY = (ushort)(muY.Bits9 * muY.Bits9);
        ushort capPosA = (ushort)((qX - 1) / 6);
        ushort capPosB = (ushort)((4096 - qY) / 6);
        ushort isPosALess = (ushort)(1 / ((capPosA / (capPosB + 1)) + 1));
        ushort capPos = (ushort)(isPosALess * capPosA + (1 - isPosALess) * capPosB);

        ushort capNegA = (ushort)((qY - 1) / 6);
        ushort capNegB = (ushort)((4096 - qX) / 6);
        ushort isNegALess = (ushort)(1 / ((capNegA / (capNegB + 1)) + 1));
        ushort capNeg = (ushort)(isNegALess * capNegA + (1 - isNegALess) * capNegB);

        ushort rawPos = (ushort)((1UL - numNeg) * qMag);
        ushort rawNegMag = (ushort)(numNeg * (qMag + hasRem));

        ushort clampPos = (ushort)((1 / ((rawPos / (capPos + 1)) + 1)) * rawPos + (1 - (1 / ((rawPos / (capPos + 1)) + 1))) * capPos);
        ushort clampNeg = (ushort)((1 / ((rawNegMag / (capNeg + 1)) + 1)) * rawNegMag + (1 - (1 / ((rawNegMag / (capNeg + 1)) + 1))) * capNeg);

        Bus36 posBus = Bus36.FromUnsignedRing9(new Ring9(clampPos));
        Bus36 negBus = -Bus36.FromUnsignedRing9(new Ring9(clampNeg));
        return posBus + negBus;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 MuNext, Bus36 DivJ) Compute12EdgeMassFluxStepAtCell(
        Ring9 muC, Ring9 phiC,
        Ring9 mPos0, Ring9 pPos0, Ring9 mNeg0, Ring9 pNeg0,
        Ring9 mPos1, Ring9 pPos1, Ring9 mNeg1, Ring9 pNeg1,
        Ring9 mPos2, Ring9 pPos2, Ring9 mNeg2, Ring9 pNeg2,
        Ring9 mPos3, Ring9 pPos3, Ring9 mNeg3, Ring9 pNeg3,
        Ring9 mPos4, Ring9 pPos4, Ring9 mNeg4, Ring9 pNeg4,
        Ring9 mPos5, Ring9 pPos5, Ring9 mNeg5, Ring9 pNeg5)
    {
        Bus36 jOut0 = ComputeSingleDirectedEdgeFlux36(muC, phiC, mPos0, pPos0);
        Bus36 jIn0  = ComputeSingleDirectedEdgeFlux36(mNeg0, pNeg0, muC, phiC);
        Bus36 jOut1 = ComputeSingleDirectedEdgeFlux36(muC, phiC, mPos1, pPos1);
        Bus36 jIn1  = ComputeSingleDirectedEdgeFlux36(mNeg1, pNeg1, muC, phiC);
        Bus36 jOut2 = ComputeSingleDirectedEdgeFlux36(muC, phiC, mPos2, pPos2);
        Bus36 jIn2  = ComputeSingleDirectedEdgeFlux36(mNeg2, pNeg2, muC, phiC);
        Bus36 jOut3 = ComputeSingleDirectedEdgeFlux36(muC, phiC, mPos3, pPos3);
        Bus36 jIn3  = ComputeSingleDirectedEdgeFlux36(mNeg3, pNeg3, muC, phiC);
        Bus36 jOut4 = ComputeSingleDirectedEdgeFlux36(muC, phiC, mPos4, pPos4);
        Bus36 jIn4  = ComputeSingleDirectedEdgeFlux36(mNeg4, pNeg4, muC, phiC);
        Bus36 jOut5 = ComputeSingleDirectedEdgeFlux36(muC, phiC, mPos5, pPos5);
        Bus36 jIn5  = ComputeSingleDirectedEdgeFlux36(mNeg5, pNeg5, muC, phiC);

        Bus36 divJ = (jOut0 - jIn0) + (jOut1 - jIn1) + (jOut2 - jIn2) +
                     (jOut3 - jIn3) + (jOut4 - jIn4) + (jOut5 - jIn5);

        ushort qCurr = (ushort)(muC.Bits9 * muC.Bits9);
        Bus36 qCurrBus = new(new Ring9(qCurr), new Ring9((ushort)(qCurr >> 9)), Ring9.Zero, Ring9.Zero);
        Bus36 qNextBus = qCurrBus - divJ;
        ushort qNext16 = (ushort)(qNextBus.L0.Bits9 | (qNextBus.L1.Bits9 << 9));
        return (Isqrt4096(qNext16), divJ);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 PauliExclusionGate9(Ring9 mu0, Ring9 phi0, Ring9 mu1, Ring9 phi1, Ring9 enabledFlag)
    {
        ushort d0 = (ushort)(1 - (1 / ((mu0.Bits9 / 46) + 1)));
        ushort d1 = (ushort)(1 - (1 / ((mu1.Bits9 / 46) + 1)));
        ushort par = (ushort)(1 / (((phi1 - phi0).AbsWrap.Bits9 / 41) + 1));
        ushort fire = (ushort)(enabledFlag.Bits9 * d0 * d1 * par);
        return new Ring9((ushort)(256 * fire));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 CascadeAttenuatedStage9(Ring9 phiIn, Ring9 impedanceRatio)
    {
        ushort r = (ushort)((1 - (1 / (impedanceRatio.Bits9 + 1))) * impedanceRatio.Bits9 + (1 / (impedanceRatio.Bits9 + 1)));
        ushort mag = phiIn.AbsWrap.Bits9;
        ushort q = (ushort)(mag / r);
        ushort isNeg = phiIn.BitQ1.Bits9;
        Ring9 divided = new((ushort)((1 - 2 * isNeg) * q));
        return ApplyHeisenbergDeadband(divided);
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

    /// <summary>
    /// Forward 3-layer nilpotent cross-channel perturbation `(I + N_pert)` on `Ring9^3` (`N_pert^3 == 0` in `d = 3`).
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 Y0, Ring9 Y1, Ring9 Y2) ApplyNilpotentPerturbation3(Ring9 s0, Ring9 s1, Ring9 s2) =>
        (s0 + s1, s1 + s2, s2);

    /// <summary>
    /// Exact 2nd-order finite Neumann perturbation inverse `(I + N_pert)^-1 = I - N_pert + N_pert^2` on `Ring9^3`:
    /// recovers `(s0, s1, s2)` with ZERO infinite series and ZERO division!
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 S0, Ring9 S1, Ring9 S2) CalibrateFiniteNeumannOrder2(Ring9 y0, Ring9 y1, Ring9 y2) =>
        (y0 - y1 + y2, y1 - y2, y2);

    /// <summary>
    /// Complete 10-Node Kitchen Bench Pre-Compensation Operator on `Ring9`:
    /// Combines target resonance phase (`139`), MOSFET Miller gate lag (`25`), Copper coil inductive lag (`64`),
    /// thermal Joule resistance drift `deltaPhiThermal`, dish wall shift (`0` for SiO2 glass, `7` for polystyrene),
    /// and individual hand-wound coil tolerance `coilSkewK` into the exact 9-bit FPGA register value!
    /// </summary>
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
        Ring9 wallPhotons = HeisenbergEnergyQuanta(singleRay);
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
    public static (Ring9 EdgeWave9, Ring9 BranchCarry9) EdgeDifferenceWithBranchCarry(Ring9 phiY, Ring9 phiX)
    {
        Bus36 yBus = Bus36.FromUnsignedRing9(phiY);
        Bus36 xBus = Bus36.FromUnsignedRing9(phiX);
        Bus36 bias255 = Bus36.FromUnsignedRing9(new Ring9(255));
        Bus36 shifted = (yBus - xBus) + bias255;
        return (phiY - phiX, shifted.L1);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 PlaquetteWinding9(Ring9 p00, Ring9 p01, Ring9 p11, Ring9 p10)
    {
        Ring9 k01 = EdgeDifferenceWithBranchCarry(p01, p00).BranchCarry9;
        Ring9 k12 = EdgeDifferenceWithBranchCarry(p11, p01).BranchCarry9;
        Ring9 k32 = EdgeDifferenceWithBranchCarry(p11, p10).BranchCarry9;
        Ring9 k03 = EdgeDifferenceWithBranchCarry(p10, p00).BranchCarry9;
        return -((k01 + k12) - (k32 + k03));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 ComputeTwoLaneCellGate(
        Ring9 mu0, Ring9 phi0, Ring9 mu1, Ring9 phi1,
        Ring9 m0_0,  Ring9 p0_0,  Ring9 m0_1,  Ring9 p0_1,  Ring9 m0_2,  Ring9 p0_2,  Ring9 m0_3,  Ring9 p0_3,
        Ring9 m0_4,  Ring9 p0_4,  Ring9 m0_5,  Ring9 p0_5,  Ring9 m0_6,  Ring9 p0_6,  Ring9 m0_7,  Ring9 p0_7,
        Ring9 m0_8,  Ring9 p0_8,  Ring9 m0_9,  Ring9 p0_9,  Ring9 m0_10, Ring9 p0_10, Ring9 m0_11, Ring9 p0_11,
        Ring9 m1_0,  Ring9 p1_0,  Ring9 m1_1,  Ring9 p1_1,  Ring9 m1_2,  Ring9 p1_2,  Ring9 m1_3,  Ring9 p1_3,
        Ring9 m1_4,  Ring9 p1_4,  Ring9 m1_5,  Ring9 p1_5,  Ring9 m1_6,  Ring9 p1_6,  Ring9 m1_7,  Ring9 p1_7,
        Ring9 m1_8,  Ring9 p1_8,  Ring9 m1_9,  Ring9 p1_9,  Ring9 m1_10, Ring9 p1_10, Ring9 m1_11, Ring9 p1_11,
        Ring9 pauliEnabledFlag)
    {
        Ring9 g0 = Compute12NeighborGateFromCells(
            mu0, phi0,
            m0_0, p0_0, m0_1, p0_1, m0_2, p0_2, m0_3, p0_3,
            m0_4, p0_4, m0_5, p0_5, m0_6, p0_6, m0_7, p0_7,
            m0_8, p0_8, m0_9, p0_9, m0_10, p0_10, m0_11, p0_11);
        Ring9 g1 = Compute12NeighborGateFromCells(
            mu1, phi1,
            m1_0, p1_0, m1_1, p1_1, m1_2, p1_2, m1_3, p1_3,
            m1_4, p1_4, m1_5, p1_5, m1_6, p1_6, m1_7, p1_7,
            m1_8, p1_8, m1_9, p1_9, m1_10, p1_10, m1_11, p1_11);
        Ring9 lead = SelectLeadingGate(g0, g1);
        Ring9 pauli = PauliExclusionGate9(mu0, phi0, mu1, phi1, pauliEnabledFlag);
        return lead + pauli;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 Compute12NeighborGateFromCells(
        Ring9 muC, Ring9 phiC,
        Ring9 m0,  Ring9 p0,  Ring9 m1,  Ring9 p1,  Ring9 m2,  Ring9 p2,  Ring9 m3,  Ring9 p3,
        Ring9 m4,  Ring9 p4,  Ring9 m5,  Ring9 p5,  Ring9 m6,  Ring9 p6,  Ring9 m7,  Ring9 p7,
        Ring9 m8,  Ring9 p8,  Ring9 m9,  Ring9 p9,  Ring9 m10, Ring9 p10, Ring9 m11, Ring9 p11)
    {
        Ring9 d0  = p0  - phiC;
        Ring9 d1  = p1  - phiC;
        Ring9 d2  = p2  - phiC;
        Ring9 d3  = p3  - phiC;
        Ring9 d4  = p4  - phiC;
        Ring9 d5  = p5  - phiC;
        Ring9 d6  = p6  - phiC;
        Ring9 d7  = p7  - phiC;
        Ring9 d8  = p8  - phiC;
        Ring9 d9  = p9  - phiC;
        Ring9 d10 = p10 - phiC;
        Ring9 d11 = p11 - phiC;

        Bus36 sumMu = Bus36.FromUnsignedRing9(m0)  + Bus36.FromUnsignedRing9(m1)  + Bus36.FromUnsignedRing9(m2)  +
                      Bus36.FromUnsignedRing9(m3)  + Bus36.FromUnsignedRing9(m4)  + Bus36.FromUnsignedRing9(m5)  +
                      Bus36.FromUnsignedRing9(m6)  + Bus36.FromUnsignedRing9(m7)  + Bus36.FromUnsignedRing9(m8)  +
                      Bus36.FromUnsignedRing9(m9)  + Bus36.FromUnsignedRing9(m10) + Bus36.FromUnsignedRing9(m11);

        Bus36 zrRaw = NeighborTerm(m0, d0.ProjC)   + NeighborTerm(m1, d1.ProjC)   + NeighborTerm(m2, d2.ProjC)   +
                      NeighborTerm(m3, d3.ProjC)   + NeighborTerm(m4, d4.ProjC)   + NeighborTerm(m5, d5.ProjC)   +
                      NeighborTerm(m6, d6.ProjC)   + NeighborTerm(m7, d7.ProjC)   + NeighborTerm(m8, d8.ProjC)   +
                      NeighborTerm(m9, d9.ProjC)   + NeighborTerm(m10, d10.ProjC) + NeighborTerm(m11, d11.ProjC);

        Bus36 ziRaw = NeighborTerm(m0, d0.ProjS)   + NeighborTerm(m1, d1.ProjS)   + NeighborTerm(m2, d2.ProjS)   +
                      NeighborTerm(m3, d3.ProjS)   + NeighborTerm(m4, d4.ProjS)   + NeighborTerm(m5, d5.ProjS)   +
                      NeighborTerm(m6, d6.ProjS)   + NeighborTerm(m7, d7.ProjS)   + NeighborTerm(m8, d8.ProjS)   +
                      NeighborTerm(m9, d9.ProjS)   + NeighborTerm(m10, d10.ProjS) + NeighborTerm(m11, d11.ProjS);

        Bus36 zr = zrRaw * muC;
        Bus36 zi = ziRaw * muC;
        Bus36 deficit = ((sumMu * muC) * Ring9.QSquare) - zr;

        Ring9 ziNeg = zi.SignBit;
        Ring9 ziPos = Ring9.MuVac - ziNeg;
        Bus36 signedDeficit = (deficit * ziPos) + ((-deficit) * ziNeg);
        Bus36 zDrive = zi + signedDeficit;

        Bus36 r2 = (Bus36.FromUnsignedRing9(muC) * muC) * Ring9.QSquare;
        Bus36 rawDen = zr + r2 + Bus36.EpsDenom;
        Ring9 denPos = (rawDen - Bus36.One).IsNonNegative;
        Ring9 denNeg = Ring9.MuVac - denPos;
        Bus36 safeDen = (rawDen * denPos) + (Bus36.EpsDenom * denNeg);

        Bus36 num = zDrive * Ring9.KP;
        Ring9 rawGate = Bus36.FloorDivideToRing9(num, safeDen);
        return ApplyHeisenbergDeadband(rawGate);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 Compute12NeighborHolonomyGateAt(Ring9 dz, Ring9 dy, Ring9 dx, Ring9 lane)
    {
        return Compute12NeighborGateFromCells(
            Iv263MuAt(dz, dy, dx), Iv263PhiAt(dz, dy, dx, lane),
            Iv263MuAt(dz,                dy - Ring9.MuVac, dx - Ring9.MuVac), Iv263PhiAt(dz,                dy - Ring9.MuVac, dx - Ring9.MuVac, lane),
            Iv263MuAt(dz,                dy - Ring9.MuVac, dx + Ring9.MuVac), Iv263PhiAt(dz,                dy - Ring9.MuVac, dx + Ring9.MuVac, lane),
            Iv263MuAt(dz,                dy + Ring9.MuVac, dx - Ring9.MuVac), Iv263PhiAt(dz,                dy + Ring9.MuVac, dx - Ring9.MuVac, lane),
            Iv263MuAt(dz,                dy + Ring9.MuVac, dx + Ring9.MuVac), Iv263PhiAt(dz,                dy + Ring9.MuVac, dx + Ring9.MuVac, lane),
            Iv263MuAt(dz - Ring9.MuVac, dy,                dx - Ring9.MuVac), Iv263PhiAt(dz - Ring9.MuVac, dy,                dx - Ring9.MuVac, lane),
            Iv263MuAt(dz - Ring9.MuVac, dy,                dx + Ring9.MuVac), Iv263PhiAt(dz - Ring9.MuVac, dy,                dx + Ring9.MuVac, lane),
            Iv263MuAt(dz + Ring9.MuVac, dy,                dx - Ring9.MuVac), Iv263PhiAt(dz + Ring9.MuVac, dy,                dx - Ring9.MuVac, lane),
            Iv263MuAt(dz + Ring9.MuVac, dy,                dx + Ring9.MuVac), Iv263PhiAt(dz + Ring9.MuVac, dy,                dx + Ring9.MuVac, lane),
            Iv263MuAt(dz - Ring9.MuVac, dy - Ring9.MuVac, dx),                Iv263PhiAt(dz - Ring9.MuVac, dy - Ring9.MuVac, dx,                lane),
            Iv263MuAt(dz - Ring9.MuVac, dy + Ring9.MuVac, dx),                Iv263PhiAt(dz - Ring9.MuVac, dy + Ring9.MuVac, dx,                lane),
            Iv263MuAt(dz + Ring9.MuVac, dy - Ring9.MuVac, dx),                Iv263PhiAt(dz + Ring9.MuVac, dy - Ring9.MuVac, dx,                lane),
            Iv263MuAt(dz + Ring9.MuVac, dy + Ring9.MuVac, dx),                Iv263PhiAt(dz + Ring9.MuVac, dy + Ring9.MuVac, dx,                lane));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 SelectLeadingGate(Ring9 g0, Ring9 g1)
    {
        ushort a0 = g0.AbsWrap.Bits9;
        ushort a1 = g1.AbsWrap.Bits9;
        ushort diff = (ushort)(1 - (1 / ((a1 / (a0 + 1)) + 1)));
        return new Ring9((ushort)(diff * g1.Bits9 + (1 - diff) * g0.Bits9));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static VerletCell72 EvaluateClosedFormWorldlineCellAt(Ring9 dz, Ring9 dy, Ring9 dx, Ring9 t, Ring9 atomicSlopeV0 = default)
    {
        Ring9 mu = Iv263MuAt(dz, dy, dx);
        Ring9 p0Init = Iv263PhiAt(dz, dy, dx, Ring9.Zero);
        Ring9 p1Init = Iv263PhiAt(dz, dy, dx, Ring9.MuVac);
        Ring9 gate0 = Compute12NeighborHolonomyGateAt(dz, dy, dx, Ring9.Zero);
        Ring9 gate1 = Compute12NeighborHolonomyGateAt(dz, dy, dx, Ring9.MuVac);
        Ring9 leadGate = SelectLeadingGate(gate0, gate1);

        var cell0 = VerletCell72.Create(mu, p0Init, p0Init - atomicSlopeV0, mu, p1Init, p1Init - atomicSlopeV0);
        cell0.Gate0 = leadGate;
        cell0.Gate1 = leadGate;
        return cell0.EvaluateClosedFormAtTick(t);
    }
}
