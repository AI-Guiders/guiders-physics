using System;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Security.Cryptography;

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
        // Direct x86_64 hardware silicon thermal TRNG instruction:
        //   48 0F C7 F8    rdseed rax
    // Direct x86_64 RDSEED instruction bytes (rdseed rax; jnc -6; and eax, 0x1FF; ret)
        //   25 FF 01 00 00 and eax, 0x1FF (native 9-bit Ring9 mask!)
        //   C3             ret
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
    public static Ring9 Iv263PhiAt(Ring9 dz, Ring9 dy, Ring9 dx, Ring9 lane)
    {
        Ring9 inBall = Iv263InBall9(dz, dy, dx);
        ushort zPos = (ushort)(1 - dz.BitQ1.Bits9);
        ushort yPos = (ushort)(1 - dy.BitQ1.Bits9);
        ushort xPos = (ushort)(1 - dx.BitQ1.Bits9);
        ushort oct = (ushort)((zPos << 2) | (yPos << 1) | xPos);
        ushort isLane1 = (ushort)(lane.Bits9 & 1);
        ushort shift = (ushort)(oct ^ (3 * isLane1));
        ushort bit = (ushort)(((263 % 256) >> shift) & 1);
        Ring9 bg = dx * Ring9.DPhiDisc;
        return bg + new Ring9((ushort)(123 * inBall.Bits9 * bit));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 Compute12NeighborHolonomyGateAt(Ring9 dz, Ring9 dy, Ring9 dx, Ring9 lane)
    {
        Ring9 muC = Iv263MuAt(dz, dy, dx);
        Ring9 phiC = Iv263PhiAt(dz, dy, dx, lane);
        Ring9 d0  = Iv263PhiAt(dz,                dy - Ring9.MuVac, dx - Ring9.MuVac, lane) - phiC;
        Ring9 d1  = Iv263PhiAt(dz,                dy - Ring9.MuVac, dx + Ring9.MuVac, lane) - phiC;
        Ring9 d2  = Iv263PhiAt(dz,                dy + Ring9.MuVac, dx - Ring9.MuVac, lane) - phiC;
        Ring9 d3  = Iv263PhiAt(dz,                dy + Ring9.MuVac, dx + Ring9.MuVac, lane) - phiC;
        Ring9 d4  = Iv263PhiAt(dz - Ring9.MuVac, dy,                dx - Ring9.MuVac, lane) - phiC;
        Ring9 d5  = Iv263PhiAt(dz - Ring9.MuVac, dy,                dx + Ring9.MuVac, lane) - phiC;
        Ring9 d6  = Iv263PhiAt(dz + Ring9.MuVac, dy,                dx - Ring9.MuVac, lane) - phiC;
        Ring9 d7  = Iv263PhiAt(dz + Ring9.MuVac, dy,                dx + Ring9.MuVac, lane) - phiC;
        Ring9 d8  = Iv263PhiAt(dz - Ring9.MuVac, dy - Ring9.MuVac, dx,                lane) - phiC;
        Ring9 d9  = Iv263PhiAt(dz - Ring9.MuVac, dy + Ring9.MuVac, dx,                lane) - phiC;
        Ring9 d10 = Iv263PhiAt(dz + Ring9.MuVac, dy - Ring9.MuVac, dx,                lane) - phiC;
        Ring9 d11 = Iv263PhiAt(dz + Ring9.MuVac, dy + Ring9.MuVac, dx,                lane) - phiC;

        Ring9 sumS = d0.ProjS + d1.ProjS + d2.ProjS + d3.ProjS + d4.ProjS + d5.ProjS +
                     d6.ProjS + d7.ProjS + d8.ProjS + d9.ProjS + d10.ProjS + d11.ProjS;
        Ring9 sumC = d0.ProjC + d1.ProjC + d2.ProjC + d3.ProjC + d4.ProjC + d5.ProjC +
                     d6.ProjC + d7.ProjC + d8.ProjC + d9.ProjC + d10.ProjC + d11.ProjC;

        Ring9 rawGate = (muC * (sumS + Ring9.QSquare - sumC)) * Ring9.BezoutInv97;
        return ApplyHeisenbergDeadband(rawGate);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static VerletCell72 EvaluateClosedFormWorldlineCellAt(Ring9 dz, Ring9 dy, Ring9 dx, Ring9 t, Ring9 atomicSlopeV0 = default)
    {
        Ring9 mu = Iv263MuAt(dz, dy, dx);
        Ring9 p0Init = Iv263PhiAt(dz, dy, dx, Ring9.Zero);
        Ring9 p1Init = Iv263PhiAt(dz, dy, dx, Ring9.MuVac);
        Ring9 gate0 = Compute12NeighborHolonomyGateAt(dz, dy, dx, Ring9.Zero);
        Ring9 gate1 = Compute12NeighborHolonomyGateAt(dz, dy, dx, Ring9.MuVac);

        var cell0 = VerletCell72.Create(mu, p0Init, p0Init - atomicSlopeV0, mu, p1Init, p1Init - atomicSlopeV0);
        cell0.Gate0 = gate0;
        cell0.Gate1 = gate1;
        return cell0.EvaluateClosedFormAtTick(t);
    }
}
