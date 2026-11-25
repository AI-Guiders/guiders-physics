using System.Runtime.CompilerServices;

namespace Guiders.Physics.Core;

public static class CuboctahedronAutomatonKernel
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 HeisenbergEnergyQuanta(Ring9 gateRaw)
    {
        ushort abs9 = gateRaw.AbsWrap.Bits9;
        return new Ring9((ushort)(abs9 / 41));
    }

    /// <summary>
    /// Branchlessly scores a 9-bit phase residual against the Heisenberg-41 deadband (`DPhiDisc = 41`):
    /// returns `4` (`3 + 1`) for exact zero stasis (`delta == 0`), `1` for sub-threshold deadband (`1 &lt;= |delta| &lt;= 40`), and `0` otherwise.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 ScoreHeisenbergStasis9(Ring9 delta9)
    {
        ushort abs9 = delta9.AbsWrap.Bits9;
        ushort ne = HeisenbergEnergyQuanta(delta9).Bits9;
        ushort inDeadband = (ushort)(1 / (ne + 1));
        ushort isZero = (ushort)(1 / (abs9 + 1));
        return new Ring9((ushort)(3 * isZero + inDeadband));
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
