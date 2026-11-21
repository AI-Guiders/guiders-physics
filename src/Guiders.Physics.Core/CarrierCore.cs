using System;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace Guiders.Physics.Core;

/// <summary>
/// Native 9-Bit Modular Ring Element in Z_512 (bits [8:0]):
/// Every input, output, coordinate, time tick, projection, and sensor value in the carrier
/// lives natively on `Ring9`!
/// </summary>
public readonly record struct Ring9
{
    private readonly ushort _bits9; // Strictly 9 bits [0..511]

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public Ring9(ushort bits9) => _bits9 = (ushort)(bits9 & 0x1FF);

    public ushort Bits9
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _bits9;
    }

    // Canonical carrier constants derived from D = 3 as native 9-bit Ring9 elements:
    public static Ring9 Zero        => new(0);
    public static Ring9 MuVac       => new(1);
    public static Ring9 C0OverCSq   => new(2);
    public static Ring9 D           => new(3);
    public static Ring9 DInSquare   => new(4);
    public static Ring9 DOutEdge    => new(5);
    public static Ring9 NSq         => new(6);
    public static Ring9 DOutVertex  => new(7);
    public static Ring9 NTri        => new(8);
    public static Ring9 ModBits     => new(9);
    public static Ring9 NLinks      => new(12);
    public static Ring9 NPhi        => new(13);
    public static Ring9 SeamTicks   => new(21);
    public static Ring9 BezoutInv41 => new(25);  // 25 * 41 === 1 (mod 512)
    public static Ring9 DPhiDisc    => new(41);  // Heisenberg threshold 41
    public static Ring9 KP          => new(64);  // Pauli amplitude ceiling 64
    public static Ring9 MHat        => new(97);  // 1 + 12 * 8 = 97
    public static Ring9 QSquare     => new(128); // Quadrant saddle 128
    public static Ring9 CollapsePhi => new(139); // 128 + 12 - 1 = 139
    public static Ring9 BezoutInv139=> new(183); // 183 * 139 === 1 (mod 512)
    public static Ring9 BezoutInv13 => new(197); // 197 * 13 === 1 (mod 512)
    public static Ring9 PhiPauli    => new(256); // Pauli half-turn 256
    public static Ring9 IvSeed      => new(263); // 0o407 = 256 + 7 = 263
    public static Ring9 BezoutInv97 => new(417); // 417 * 97 === 1 (mod 512)
    public static Ring9 MinusOne    => new(511); // -1 mod 512

    // Canonical 12 FCC Cuboctahedron O_2(3) neighbor offsets in pure Ring9^3:
    public static readonly (Ring9 Dz, Ring9 Dy, Ring9 Dx)[] CuboctOffsets12 =
    [
        (Zero,     MuVac,    MuVac),    (Zero,     MuVac,    MinusOne),
        (Zero,     MinusOne, MuVac),    (Zero,     MinusOne, MinusOne),
        (MuVac,    Zero,     MuVac),    (MuVac,    Zero,     MinusOne),
        (MinusOne, Zero,     MuVac),    (MinusOne, Zero,     MinusOne),
        (MuVac,    MuVac,    Zero),     (MuVac,    MinusOne, Zero),
        (MinusOne, MuVac,    Zero),     (MinusOne, MinusOne, Zero)
    ];

    /// <summary>Half-turn Pauli bit q1 = bit 8 (in {0, 1} as Ring9).</summary>
    public Ring9 BitQ1
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new((ushort)((_bits9 >> 8) & 1));
    }

    /// <summary>Quadrant bit q0 = bit 7 (in {0, 1} as Ring9).</summary>
    public Ring9 BitQ0
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new((ushort)((_bits9 >> 7) & 1));
    }

    /// <summary>7-bit intra-quadrant mantissa m7 = bits [6:0] (in 0..127 as Ring9).</summary>
    public Ring9 Mantissa7
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new((ushort)(_bits9 & 127));
    }

    /// <summary>Pauli singlet partner (phi + 256) mod 512 = {~phi[8], phi[7:0]}.</summary>
    public Ring9 PauliPartner
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new((ushort)(_bits9 ^ 256));
    }

    /// <summary>
    /// 100% branchless Gray-Quadrant Butterfly projections (C_3(phi), S_3(phi)) on pure 9-bit Ring9 wires!
    /// Matches `rotator_512_unit` in `hardware/cuboctahedron_carrier_asic.v` gate-for-gate:
    ///   comp7 = 128 - m7
    ///   X-butterfly swap by q0: magC = (1-q0)*comp7 + q0*m7, magS = (1-q0)*m7 + q0*comp7
    ///   Gray signs: c_neg = q1 ^ q0, s_neg = q1.
    /// </summary>
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

    /// <summary>
    /// 9-bit folded magnitude `|wrap512(phi)| = 128 - C_3(phi)` in `[0..256]` as a native `Ring9`!
    /// (Matches Module 1B `abs_x_out = c0 - cx` in `cuboctahedron_carrier_asic.v`!)
    /// </summary>
    public Ring9 AbsWrap
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get
        {
            ushort negU = (ushort)((256 - _bits9) & 0x1FF);
            // w = 256 - negU; |w| = C_3(0) - C_3(phi)
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

    /// <summary>
    /// Coordinate-free Bezout inverse of Heisenberg quantum 41 (`25 * 41 === 1 mod 512`) on `Ring9`.
    /// </summary>
    public Ring9 DecodeBezout41
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => this * BezoutInv41;
    }

    /// <summary>
    /// Closed-form O(1) Verlet Worldline Parabola Jump purely on 9-bit `Ring9` (`verletParabolaClosedFormJumpOp`):
    /// Jumps directly from `this = phi(0)` to ANY time tick `t` (`Ring9`) in 1 step:
    /// `phi(t) = phi(0) + t * v0 + (t * (t + 1) / 2) * gatePhi (mod 512)`.
    /// Note that on Z_512, `(t * (t + 1)) >> 1` mod 512 is a pure 10-bit product shifted right by 1 bit!
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public Ring9 JumpVerletParabola(Ring9 v0, Ring9 gatePhi, Ring9 t)
    {
        ushort t10 = t._bits9;
        ushort tri9 = (ushort)(((t10 * (t10 + 1)) >> 1) & 0x1FF);
        Ring9 triT = new(tri9);
        return this + (t * v0) + (triT * gatePhi);
    }

    /// <summary>Octal triad (o2, o1, o0) in Z_8^3 (each 3-bit digit as Ring9).</summary>
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
/// Exact 72-Bit (9-Byte, Pack = 1) Two-Tick Verlet Polar Cell:
/// 8 lanes of 9-bit `Ring9` = 72 bits = 9 bytes exact, with zero padding and zero branches!
/// </summary>
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

    /// <summary>
    /// Closed-form O(1) Worldline Evaluation of this 72-bit cell at ANY 9-bit time tick `t` (`Ring9`)!
    /// Zero loops, zero `int` — pure `Ring9` parabola `JumpVerletParabola`!
    /// </summary>
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
}

/// <summary>
/// Pure 9-Bit `Ring9` Multimodal Sensor Readout from a 72-bit `VerletCell72`:
/// Every channel (ColorR, ColorG, ColorB, WavelengthShift, AcousticWave, ReactionCurrent, EnergyQuantaNe) is a `Ring9`!
/// </summary>
public readonly record struct MultimodalReadout9(
    Ring9 ColorR,
    Ring9 ColorG,
    Ring9 ColorB,
    Ring9 AcousticWave,
    Ring9 ReactionCurrent,
    Ring9 EnergyQuantaNe);

/// <summary>
/// 100% Branchless Pure 9-Bit `Ring9` Closed-Form Worldline &amp; Rotator Kernel.
/// Every parameter and return value is strictly `Ring9` or `VerletCell72`!
/// </summary>
public static class CuboctahedronAutomatonKernel
{
    /// <summary>
    /// Pure 9-bit branchless step indicator `1` if `a.AbsWrap >= 41` (above Heisenberg deadband), else `0` as `Ring9`.
    /// </summary>
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

    /// <summary>
    /// Pure 9-bit `Ring9` Multimodal Readout (RGB triad with R+G+B = 128 on Ring9, acoustic wave, reaction current, thermal quanta).
    /// Zero `int`, zero `if`, zero `?:`!
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static MultimodalReadout9 ReadCellSensors(VerletCell72 cell)
    {
        Ring9 dphi = cell.Phi1Curr - cell.Phi0Curr;
        Ring9 q1 = dphi.BitQ1;
        Ring9 q0 = dphi.BitQ0;
        Ring9 m7 = dphi.Mantissa7;
        Ring9 comp7 = Ring9.QSquare - m7;

        Ring9 magC = (Ring9.MuVac - q0) * comp7 + q0 * m7;
        Ring9 magS = (Ring9.MuVac - q0) * m7 + q0 * comp7;
        Ring9 cNeg = q1 ^ q0;

        // Branchless RGB triad on Ring9: R = (1 - cNeg) * magC, G = magS, B = cNeg * magC -> R + G + B === 128!
        Ring9 r = (Ring9.MuVac - cNeg) * magC;
        Ring9 g = magS;
        Ring9 b = cNeg * magC;

        Ring9 c3 = dphi.ProjC;
        Ring9 s3 = dphi.ProjS;
        Ring9 sound = cell.Mu0 * cell.Mu1 * (Ring9.QSquare - c3);
        Ring9 rxn = cell.Mu0 * cell.Mu1 * s3;
        Ring9 ne = HeisenbergEnergyQuanta(cell.Gate0);

        return new MultimodalReadout9(r, g, b, sound, rxn, ne);
    }

    /// <summary>
    /// Coordinate-Free O(1) Closed-Form Worldline Cell Evaluator purely on `Ring9`:
    /// Given spatial coordinates `(dz, dy, dx) in Ring9^3` and time tick `t in Ring9`,
    /// computes the initial phase from `IV(3) = 263 = 0o407` and jumps in O(1) along the Verlet worldline parabola!
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static VerletCell72 EvaluateClosedFormWorldlineCellAt(Ring9 dz, Ring9 dy, Ring9 dx, Ring9 t)
    {
        // Octal triad & background phase ramp on Ring9:
        Ring9 bg = dx * Ring9.DPhiDisc;
        Ring9 p0Init = bg + Ring9.IvSeed;
        Ring9 p1Init = bg + Ring9.IvSeed + (Ring9.D * Ring9.DPhiDisc);
        Ring9 gate0 = Ring9.DPhiDisc + dz + dy;
        Ring9 gate1 = Ring9.CollapsePhi - dz + dy;

        var cell0 = VerletCell72.Create(Ring9.NPhi, p0Init, p0Init, Ring9.NPhi, p1Init, p1Init);
        cell0.Gate0 = gate0;
        cell0.Gate1 = gate1;
        return cell0.EvaluateClosedFormAtTick(t);
    }
}
