using System.Runtime.CompilerServices;

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
