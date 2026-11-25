using System.Runtime.CompilerServices;

namespace Guiders.Physics.Core;

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
