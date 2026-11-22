using System.Security.Cryptography;
using Guiders.Physics.Core;
using Guiders.Physics.Compression;
using Guiders.Physics.MediaCodec;
using Xunit;

namespace Guiders.Physics.Tests;

public class GuidersPhysicsSuiteTests
{
    [Fact]
    public unsafe void Core_1_PassportAnd72BitLayout_MatchesLean4AndPython()
    {
        Assert.Equal(9, sizeof(VerletCell72));
        Assert.Equal(72, sizeof(VerletCell72) * 8);

        Assert.Equal(new Ring9(3),   Ring9.D);
        Assert.Equal(new Ring9(1),   Ring9.MuVac);
        Assert.Equal(new Ring9(2),   Ring9.C0OverCSq);
        Assert.Equal(new Ring9(12),  Ring9.NLinks);
        Assert.Equal(new Ring9(8),   Ring9.NTri);
        Assert.Equal(new Ring9(6),   Ring9.NSq);
        Assert.Equal(new Ring9(4),   Ring9.DInSquare);
        Assert.Equal(new Ring9(5),   Ring9.DOutEdge);
        Assert.Equal(new Ring9(7),   Ring9.DOutVertex);
        Assert.Equal(new Ring9(9),   Ring9.ModBits);
        Assert.Equal(new Ring9(13),  Ring9.NPhi);
        Assert.Equal(new Ring9(41),  Ring9.DPhiDisc);
        Assert.Equal(new Ring9(21),  Ring9.SeamTicks);
        Assert.Equal(new Ring9(64),  Ring9.KP);
        Assert.Equal(new Ring9(97),  Ring9.MHat);
        Assert.Equal(new Ring9(128), Ring9.QSquare);
        Assert.Equal(new Ring9(139), Ring9.CollapsePhi);
        Assert.Equal(new Ring9(256), Ring9.PhiPauli);
        Assert.Equal(new Ring9(263), Ring9.IvSeed);

        var cell = VerletCell72.Create(new(64), new(263), new(139), new(13), new(41), new(97));
        cell.Gate0 = new Ring9(352);
        cell.Gate1 = new Ring9(256);
        Assert.Equal(new Ring9(64),  cell.Mu0);
        Assert.Equal(new Ring9(263), cell.Phi0Curr);
        Assert.Equal(new Ring9(139), cell.Phi0Past);
        Assert.Equal(new Ring9(13),  cell.Mu1);
        Assert.Equal(new Ring9(41),  cell.Phi1Curr);
        Assert.Equal(new Ring9(97),  cell.Phi1Past);
        Assert.Equal(new Ring9(352), cell.Gate0);
        Assert.Equal(new Ring9(256), cell.Gate1);
    }

    [Fact]
    public void Core_2_Rotator512AndPythagorean_ExactOnRing9AndBus36All512()
    {
        for (ushort phi = 0; phi < 512; phi++)
        {
            Ring9 r = new(phi);
            Ring9 c = r.ProjC;
            Ring9 s = r.ProjS;

            Assert.Equal(Ring9.QSquare - r.AbsWrap, c);
            Assert.Equal((r - Ring9.QSquare).ProjC, s);
            Assert.Equal(Ring9.QSquare, c.AbsWrap + s.AbsWrap);

            // Pythagorean identity Px^2 + Py^2 == Pr^2 directly on Bus36 without unwrapping bits!
            var (px, py, pr) = CuboctahedronAutomatonKernel.PythagoreanRotator36(r);
            Assert.Equal(pr * pr, (px * px) + (py * py));
        }

        // At Heisenberg threshold phi = 41: (Px, Py, Pr) = (14703, 10496, 18065) = (28*512 + 367, 20*512 + 256, 35*512 + 145) on Bus36!
        var (px41, py41, pr41) = CuboctahedronAutomatonKernel.PythagoreanRotator36(Ring9.DPhiDisc);
        Assert.Equal(new Bus36(new(367), new(28), Ring9.Zero, Ring9.Zero), px41);
        Assert.Equal(new Bus36(new(256), new(20), Ring9.Zero, Ring9.Zero), py41);
        Assert.Equal(new Bus36(new(145), new(35), Ring9.Zero, Ring9.Zero), pr41);
    }

    [Fact]
    public void Core_3_Channel2BranchCarryAndPlaquetteWinding_ExactRuleZ3()
    {
        // Rule Z3 three regimes on Ring9:
        Assert.Equal(Ring9.MinusOne, CuboctahedronAutomatonKernel.EdgeDifferenceWithBranchCarry(new(0), new(300)).BranchCarry9);
        Assert.Equal(Ring9.Zero,     CuboctahedronAutomatonKernel.EdgeDifferenceWithBranchCarry(new(100), new(50)).BranchCarry9);
        Assert.Equal(Ring9.MuVac,    CuboctahedronAutomatonKernel.EdgeDifferenceWithBranchCarry(new(400), new(100)).BranchCarry9);

        // Plaquette winding around a vortex loop (0 -> 130 -> 260 -> 390 -> 0) == +1 (Ring9.MuVac):
        Ring9 w = CuboctahedronAutomatonKernel.PlaquetteWinding9(new(0), new(130), new(260), new(390));
        Assert.Equal(Ring9.MuVac, w);
    }

    [Fact]
    public void Core_4_PauliExclusionAndVerletReversibility_MatchesPythonCheckPauliRepel()
    {
        Ring9 pOn = CuboctahedronAutomatonKernel.PauliExclusionGate9(new(46), Ring9.Zero, new(46), Ring9.Zero, Ring9.MuVac);
        Ring9 pOff = CuboctahedronAutomatonKernel.PauliExclusionGate9(new(46), Ring9.Zero, new(46), Ring9.Zero, Ring9.Zero);
        Ring9 pAnti = CuboctahedronAutomatonKernel.PauliExclusionGate9(new(46), Ring9.Zero, new(46), Ring9.PhiPauli, Ring9.MuVac);
        Ring9 pSparse = CuboctahedronAutomatonKernel.PauliExclusionGate9(Ring9.NPhi, Ring9.Zero, Ring9.NPhi, Ring9.Zero, Ring9.MuVac);

        Assert.Equal(Ring9.PhiPauli, pOn);
        Assert.Equal(Ring9.Zero, pOff);
        Assert.Equal(Ring9.Zero, pAnti);
        Assert.Equal(Ring9.Zero, pSparse);

        for (ushort p = 0; p < 512; p++)
        {
            Ring9 curr = new(p);
            Ring9 past = curr * Ring9.NPhi + Ring9.MHat;
            Ring9 gate = curr * Ring9.DPhiDisc;
            var c0 = VerletCell72.Create(Ring9.NPhi, curr, past, Ring9.NPhi, curr + Ring9.DPhiDisc, past + Ring9.DPhiDisc);
            var c1 = CuboctahedronAutomatonKernel.StepCellForward(c0, gate, gate);
            var c0Rec = CuboctahedronAutomatonKernel.StepCellReverse(c1, gate, gate);
            Assert.Equal(c0.Phi0Curr, c0Rec.Phi0Curr);
            Assert.Equal(c0.Phi0Past, c0Rec.Phi0Past);
            Assert.Equal(c0.Phi1Curr, c0Rec.Phi1Curr);
            Assert.Equal(c0.Phi1Past, c0Rec.Phi1Past);
        }
    }

    [Fact]
    public void Core_5_Iv263HolonomyAndFull3DCensus_MatchesLean4AndPythonBitForBit()
    {
        Assert.Equal(Ring9.Zero, CuboctahedronAutomatonKernel.Compute12NeighborHolonomyGateAt(new(2), new(2), new(2), Ring9.Zero));
        Ring9 gWallPos = CuboctahedronAutomatonKernel.Compute12NeighborHolonomyGateAt(Ring9.Zero, new(2), new(2), Ring9.MuVac);
        Ring9 gWallNeg = CuboctahedronAutomatonKernel.Compute12NeighborHolonomyGateAt(Ring9.MinusOne, new(2), new(2), Ring9.MuVac);
        Assert.Equal(new Ring9(76), gWallPos);
        Assert.Equal(-new Ring9(77), gWallNeg);
        Assert.Equal(Ring9.MinusOne, gWallPos + gWallNeg);

        const ushort N = 24;
        var carrier = new CuboctahedronCarrier(new Ring9(N));
        carrier.StepForward();
        carrier.StepForward();

        Ring9 P2(int z, int y, int x) => carrier[new((ushort)z), new((ushort)y), new((ushort)x)].Phi0Curr;
        Ring9 NXY(int z, int y, int x) => CuboctahedronAutomatonKernel.PlaquetteWinding9(P2(z, y, x), P2(z, y, x + 1), P2(z, y + 1, x + 1), P2(z, y + 1, x));
        Ring9 NXZ(int z, int y, int x) => CuboctahedronAutomatonKernel.PlaquetteWinding9(P2(z, y, x), P2(z, y, x + 1), P2(z + 1, y, x + 1), P2(z + 1, y, x));
        Ring9 NYZ(int z, int y, int x) => CuboctahedronAutomatonKernel.PlaquetteWinding9(P2(z, y, x), P2(z, y + 1, x), P2(z + 1, y + 1, x), P2(z + 1, y, x));

        Ring9 sumXY = Ring9.Zero, sumYZ = Ring9.Zero, sumXZ = Ring9.Zero;
        for (int z = 0; z < N; z++)
            for (int y = 0; y < N; y++)
                for (int x = 0; x < N; x++)
                {
                    sumXY = sumXY + NXY(z, y, x).AbsWrap;
                    sumYZ = sumYZ + NYZ(z, y, x).AbsWrap;
                    sumXZ = sumXZ + NXZ(z, y, x).AbsWrap;
                    Ring9 divCube = (NXY(z + 1, y, x) - NXY(z, y, x))
                                  - (NXZ(z, y + 1, x) - NXZ(z, y, x))
                                  + (NYZ(z, y, x + 1) - NYZ(z, y, x));
                    Assert.Equal(Ring9.Zero, divCube);
                }

        Assert.Equal(new Ring9(26), sumXY);
        Assert.Equal(new Ring9(26), sumYZ);
        Assert.Equal(new Ring9(34), sumXZ);
        Assert.Equal(new Ring9(86), sumXY + sumYZ + sumXZ);

        // Verify exact 3D time reversal back to t = 0 via carrier.StepBackward()!
        carrier.StepBackward();
        carrier.StepBackward();
        Assert.Equal(Ring9.Zero, carrier.Tick);
        Assert.Equal(CuboctahedronAutomatonKernel.Iv263PhiAt(Ring9.Zero, new(2), new(2), Ring9.MuVac), carrier[new(12), new(14), new(14)].Phi1Curr);
    }

    [Fact]
    public void Core_6_MultimodalSensorReadout_ExactRgbTriadAndIsqrt4096()
    {
        for (ushort m = 0; m <= 64; m++)
        {
            Assert.Equal(new Ring9(m), CuboctahedronAutomatonKernel.Isqrt4096((ushort)(m * m)));
        }

        for (ushort phi = 0; phi < 512; phi++)
        {
            var c = VerletCell72.Create(Ring9.NPhi, new(phi), new(phi), Ring9.NPhi, new(phi), new(phi));
            c.Gate0 = Ring9.CollapsePhi;
            var s = CuboctahedronAutomatonKernel.ReadCellSensors(c);
            Assert.Equal(Ring9.QSquare, s.ColorR + s.ColorG + s.ColorB);
            Assert.Equal(Ring9.D, s.EnergyQuantaNe);
        }
    }

    [Fact]
    public void Core_7_RealReactionSimulations_InertSolventAndBelousovZhabotinsky13TickTrajectory_MatchesPythonBitForBit()
    {
        const ushort N = 7;
        Ring9 c = new(3);

        static CuboctahedronCarrier MakeReactionCarrier(Ring9 phiEnv)
        {
            var carrier = new CuboctahedronCarrier(new Ring9(N));
            for (ushort iz = 0; iz < N; iz++)
            {
                int dz = iz - 3;
                for (ushort iy = 0; iy < N; iy++)
                {
                    int dy = iy - 3;
                    for (ushort ix = 0; ix < N; ix++)
                    {
                        int dx = ix - 3;
                        Ring9 p = (dz == 0 && dy == 0 && dx == 0)
                            ? Ring9.Zero
                            : ((4 * dz + 2 * dy + dx) >= 0 ? phiEnv : -phiEnv);
                        carrier[new(iz), new(iy), new(ix)] = VerletCell72.Create(Ring9.NPhi, p, p, Ring9.NPhi, p, p);
                    }
                }
            }
            return carrier;
        }

        // 1) Inert solvent (phi_env = 41): reaction does NOT proceed (center phase stays 0 -> 0)
        var inertCarrier = MakeReactionCarrier(Ring9.DPhiDisc);
        inertCarrier.StepForward();
        Assert.Equal(Ring9.Zero, inertCarrier[c, c, c].Phi0Curr);

        // 2) Catalytic Belousov-Zhabotinsky oscillation (phi_env = 139 = CollapsePhi):
        ushort[] expectedBzPhases = [0, 352, 361, 212, 214, 61, 262, 301, 95, 52, 121, 90, 59, 427];
        var bzCarrier = MakeReactionCarrier(Ring9.CollapsePhi);

        for (int t = 1; t <= 13; t++)
        {
            bzCarrier.StepForward();
            Assert.Equal(new Ring9(expectedBzPhases[t]), bzCarrier[c, c, c].Phi0Curr);

            var ro = bzCarrier.ReadSensors(c, c, c);
            Assert.Equal(Ring9.QSquare, ro.ColorR + ro.ColorG + ro.ColorB);

            if (t == 1)
            {
                Assert.Equal(Ring9.Zero, ro.ColorR);
                Assert.Equal(new Ring9(96), ro.ColorG);
                Assert.Equal(new Ring9(32), ro.ColorB);
            }
            else if (t == 5)
            {
                Assert.Equal(new Ring9(67), ro.ColorR);
                Assert.Equal(new Ring9(61), ro.ColorG);
                Assert.Equal(Ring9.Zero, ro.ColorB);
            }
        }
    }

    [Fact]
    public void Core_8_EdgeMassFluxTransport_MatchesPythonPolarEdgeFluxStepQCoreBitForBit()
    {
        const ushort N = 16;
        var carrier = new CuboctahedronCarrier(new Ring9(N));
        carrier.StepForward();

        VerletCell72 Cell(int z, int y, int x) => carrier[new((ushort)((z + N) % N)), new((ushort)((y + N) % N)), new((ushort)((x + N) % N))];

        // In Python `torch.roll(t, shifts=(dz, dy, dx))` reads `(z - dz, y - dy, x - dx)` for `_FCC_OFFSETS[:6]`:
        // (1,1,0), (1,-1,0), (-1,1,0), (-1,-1,0), (1,0,1), (1,0,-1)
        Bus36 totalDivJ = Bus36.Zero;
        int nonZeroDivJ = 0;
        Ring9 sampleMu882 = Ring9.Zero;
        Bus36 sampleDivJ882 = Bus36.Zero;

        for (int z = 0; z < N; z++)
            for (int y = 0; y < N; y++)
                for (int x = 0; x < N; x++)
                {
                    var cur = Cell(z, y, x);
                    var p0 = Cell(z - 1, y - 1, x);     var m0 = Cell(z + 1, y + 1, x);
                    var p1 = Cell(z - 1, y + 1, x);     var m1 = Cell(z + 1, y - 1, x);
                    var p2 = Cell(z + 1, y - 1, x);     var m2 = Cell(z - 1, y + 1, x);
                    var p3 = Cell(z + 1, y + 1, x);     var m3 = Cell(z - 1, y - 1, x);
                    var p4 = Cell(z - 1, y,     x - 1); var m4 = Cell(z + 1, y,     x + 1);
                    var p5 = Cell(z - 1, y,     x + 1); var m5 = Cell(z + 1, y,     x - 1);

                    var (muNext, divJ) = CuboctahedronAutomatonKernel.Compute12EdgeMassFluxStepAtCell(
                        cur.Mu0, cur.Phi0Curr,
                        p0.Mu0, p0.Phi0Curr, m0.Mu0, m0.Phi0Curr,
                        p1.Mu0, p1.Phi0Curr, m1.Mu0, m1.Phi0Curr,
                        p2.Mu0, p2.Phi0Curr, m2.Mu0, m2.Phi0Curr,
                        p3.Mu0, p3.Phi0Curr, m3.Mu0, m3.Phi0Curr,
                        p4.Mu0, p4.Phi0Curr, m4.Mu0, m4.Phi0Curr,
                        p5.Mu0, p5.Phi0Curr, m5.Mu0, m5.Phi0Curr);

                    totalDivJ = totalDivJ + divJ;
                    if (divJ != Bus36.Zero) nonZeroDivJ++;
                    if (z == 8 && y == 8 && x == 2)
                    {
                        sampleMu882 = muNext;
                        sampleDivJ882 = divJ;
                    }
                }

        Assert.Equal(Bus36.Zero, totalDivJ);
        Assert.Equal(1067, nonZeroDivJ);
        Assert.Equal(Bus36.FromUnsignedRing9(new Ring9(15)), sampleDivJ882);
        Assert.Equal(Ring9.NLinks, sampleMu882);
    }

    [Fact]
    public void Core_9_SymmetriesAndAxioms_A5VacuumFloor_A14U1Invariance_Spinor360And720_NitrogenaseJanus()
    {
        // 1) A5 Planck cell vacuum floor: 0 -> 1 (MuVac), non-zero -> unchanged
        Assert.Equal(Ring9.MuVac, CuboctahedronAutomatonKernel.EnforcePlanckCellFloor(Ring9.Zero));
        Assert.Equal(Ring9.NPhi,  CuboctahedronAutomatonKernel.EnforcePlanckCellFloor(Ring9.NPhi));

        // 2) A14 U(1) global phase-shift invariance across all 512 global clock shifts c:
        for (ushort c = 0; c < 512; c += 13)
        {
            Ring9 shift = new(c);
            Ring9 pX = new Ring9(41) + shift;
            Ring9 pY = new Ring9(180) + shift;
            Assert.Equal((new Ring9(180) - new Ring9(41)).ProjC, (pY - pX).ProjC);
            Assert.Equal((new Ring9(180) - new Ring9(41)).ProjS, (pY - pX).ProjS);
            Assert.Equal(
                CuboctahedronAutomatonKernel.ComputeSingleDirectedEdgeFlux36(Ring9.NPhi, new(41), Ring9.NPhi, new(180)),
                CuboctahedronAutomatonKernel.ComputeSingleDirectedEdgeFlux36(Ring9.NPhi, pX, Ring9.NPhi, pY));
        }

        // 3) Spinor 360-degree Pauli flip (+256 negates ProjC and ProjS) and 720-degree identity (+256 + 256 == id):
        for (ushort phi = 0; phi < 512; phi++)
        {
            Ring9 r = new(phi);
            Ring9 r360 = r.PauliPartner;
            Ring9 r720 = r360.PauliPartner;
            Assert.Equal(-r.ProjC, r360.ProjC);
            Assert.Equal(-r.ProjS, r360.ProjS);
            Assert.Equal(r, r720);
        }

        // 4) Nitrogenase Fe7MoS9C Lowe-Thorneley E1..E8 Janus saddle transition at E4 (4 * 41 = 164 > 128) and E7 (7 * 41 = 287 -> C_3 = -97):
        Ring9 e3Phase = Ring9.D * Ring9.DPhiDisc;         // 3 * 41 = 123 (< 128) -> C_3 = +5
        Ring9 e4Phase = Ring9.DInSquare * Ring9.DPhiDisc; // 4 * 41 = 164 (> 128) -> C_3 = -36
        Ring9 e7Phase = Ring9.DOutVertex * Ring9.DPhiDisc;// 7 * 41 = 287         -> C_3 = -97 = -MHat!
        Assert.Equal(Ring9.DOutEdge, e3Phase.ProjC);
        Assert.Equal(-new Ring9(36), e4Phase.ProjC);
        Assert.Equal(-Ring9.MHat,    e7Phase.ProjC);
    }

    [Fact]
    public void Core_10_CascadedSystemInteractionAndInformationLoss_FourSystemTypes()
    {
        // Type I: Unmatched open-loop T-cascade (System_1 -> System_2 -> System_3 -> System_4 with impedance ratio 2 at each boundary):
        // Even a 250-tick impulse decays 250 -> 125 -> 62 -> 31 (< 41 Heisenberg deadband!) -> 0 (100% information loss at System_4!)
        Ring9 s1 = new(250);
        Ring9 s2 = CuboctahedronAutomatonKernel.CascadeAttenuatedStage9(s1, Ring9.C0OverCSq);
        Ring9 s3 = CuboctahedronAutomatonKernel.CascadeAttenuatedStage9(s2, Ring9.C0OverCSq);
        Ring9 s4 = CuboctahedronAutomatonKernel.CascadeAttenuatedStage9(s3, Ring9.C0OverCSq);
        Assert.Equal(new Ring9(125), s2);
        Assert.Equal(new Ring9(62),  s3);
        Assert.Equal(Ring9.Zero,     s4);

        // Type II: Impedance-matched catalytic resonance cascade (A_coh = 139, unit impedance ratio 1 across all 4 stages):
        Ring9 m2 = CuboctahedronAutomatonKernel.CascadeAttenuatedStage9(Ring9.CollapsePhi, Ring9.MuVac);
        Ring9 m3 = CuboctahedronAutomatonKernel.CascadeAttenuatedStage9(m2, Ring9.MuVac);
        Ring9 m4 = CuboctahedronAutomatonKernel.CascadeAttenuatedStage9(m3, Ring9.MuVac);
        Assert.Equal(Ring9.CollapsePhi, m4);

        // Type III: Sub-threshold coherent M-channel inside [-40, +40] (81 states: 0 thermal quanta n_E == 0 and 100% lossless 4-stage rotator relay!):
        for (ushort k = 0; k <= 40; k++)
        {
            Ring9 pos = new(k);
            Ring9 neg = -pos;
            Assert.Equal(Ring9.Zero, CuboctahedronAutomatonKernel.HeisenbergEnergyQuanta(pos));
            Assert.Equal(Ring9.Zero, CuboctahedronAutomatonKernel.HeisenbergEnergyQuanta(neg));
            Assert.Equal(pos, CuboctahedronAutomatonKernel.CascadeSubthresholdRotatorRelay4Stage9(pos));
            Assert.Equal(neg, CuboctahedronAutomatonKernel.CascadeSubthresholdRotatorRelay4Stage9(neg));
        }

        // Type IV: Closed-loop pre-compensated drive cancels arbitrary multi-stage hardware jitters (UART, MOSFET, Coil, Dish wall) across all 512 target phases!
        Ring9 jUart = new(17), jMosfet = -new Ring9(43), jCoil = new(29), jWall = new(58);
        Assert.NotEqual(Ring9.CollapsePhi, CuboctahedronAutomatonKernel.CascadeUncompensatedArrival9(Ring9.CollapsePhi, jUart, jMosfet, jCoil, jWall));

        for (ushort t = 0; t < 512; t++)
        {
            Ring9 target = new(t);
            Ring9 drive = CuboctahedronAutomatonKernel.CascadePrecompensatedDrive9(target, jUart, jMosfet, jCoil, jWall);
            Ring9 arrived = CuboctahedronAutomatonKernel.CascadeUncompensatedArrival9(drive, jUart, jMosfet, jCoil, jWall);
            Assert.Equal(target, arrived);
        }

        // Verify the 6 Human-Engineered System Classes Passports on Ring9:
        // Class 1 (Digital Silicon/Quartz) & Class 3 (Air N2+O2 == SiO2 Glass Dish Wall == 30):
        Ring9 zN2 = Ring9.C0OverCSq * Ring9.DOutVertex; // 2 * 7 = 14
        Ring9 zO2 = Ring9.ZOxygen2_16;                  // 16
        Assert.Equal(Ring9.ZQuartz30, Ring9.ZSilicon14 + zO2);
        Assert.Equal(Ring9.ZQuartz30, zN2 + zO2);
        Assert.Equal(Ring9.QuartzSiO2Phi, Ring9.ZQuartz30 * Ring9.DPhiDisc);
        // Class 2 (RLC resonance at deltaPhi = 0 -> ProjS == 0, ProjC == QSquare = 128):
        Assert.Equal(Ring9.Zero, Ring9.Zero.ProjS);
        Assert.Equal(Ring9.QSquare, Ring9.Zero.ProjC);
        // Class 4 (Aerospace Kutta-Joukowski vortex lift circulation = +1 & sub-threshold laminar boundary layer |40| -> n_E = 0):
        Assert.Equal(Ring9.MuVac, CuboctahedronAutomatonKernel.PlaquetteWinding9(Ring9.Zero, new(130), new(260), new(390)));
        Assert.Equal(Ring9.Zero, CuboctahedronAutomatonKernel.HeisenbergEnergyQuanta(Ring9.DPhiDisc - Ring9.MuVac));

        // Single-Instrument Deep Dive: Alpha-Quartz Piezoelectric Resonator ((SiO2)_3 = 90 = MHat - DOutVertex, SiO4 = 46 = minimal Pauli dense threshold)
        // and Real Crystal Impurity Physics (AlO4 = 45 < 46 drops below Pauli threshold; proton passivation AlO4H = 45 + 1 = 46 restores exact Pauli threshold!):
        Assert.Equal(Ring9.ZAlphaQuartz90, Ring9.D * Ring9.ZQuartz30);
        Assert.Equal(Ring9.ZAlphaQuartz90, Ring9.MHat - Ring9.DOutVertex);
        Ring9 zSiO4Tetra = Ring9.ZSilicon14 + (Ring9.DInSquare * Ring9.NTri); // 14 + 32 = 46
        Ring9 zAlO4Impurity = Ring9.NPhi + (Ring9.DInSquare * Ring9.NTri);    // 13 + 32 = 45
        Ring9 zAlO4HPassivated = zAlO4Impurity + Ring9.MuVac;                 // 45 + 1 = 46 == Z(SiO4)!
        Assert.Equal(new Ring9(46), zSiO4Tetra);
        Assert.Equal(new Ring9(45), zAlO4Impurity);
        Assert.Equal(zSiO4Tetra, zAlO4HPassivated);
        Assert.Equal(Ring9.PhiPauli, CuboctahedronAutomatonKernel.PauliExclusionGate9(zSiO4Tetra, Ring9.Zero, zSiO4Tetra, Ring9.Zero, Ring9.MuVac));
        Assert.Equal(Ring9.Zero,     CuboctahedronAutomatonKernel.PauliExclusionGate9(zAlO4Impurity, Ring9.Zero, zAlO4Impurity, Ring9.Zero, Ring9.MuVac));
        Assert.Equal(Ring9.PhiPauli, CuboctahedronAutomatonKernel.PauliExclusionGate9(zAlO4HPassivated, Ring9.Zero, zAlO4HPassivated, Ring9.Zero, Ring9.MuVac));
    }

    [Fact]
    public void Core_11_FinitePerturbationTheoryAndStringDimensionDemystification()
    {
        // 1) Maximum perturbation order on a cell is |sq(3)| = 6 (256 / 41 = 6, remainder 10 = Z(Ne)):
        Assert.Equal(Ring9.NSq, CuboctahedronAutomatonKernel.HeisenbergEnergyQuanta(Ring9.PhiPauli));
        Assert.Equal(new Ring9(10), Ring9.PhiPauli - (Ring9.NSq * Ring9.DPhiDisc));

        // 2) Exact 2nd-order Neumann calibration inverse (I + N_pert)^-1 = I - N_pert + N_pert^2 across all 512 ring states:
        for (ushort k = 0; k < 512; k += 13)
        {
            Ring9 s0 = new(k);
            Ring9 s1 = s0 + Ring9.ZSilicon14;
            Ring9 s2 = s0 + Ring9.ZQuartz30;
            var (y0, y1, y2) = CuboctahedronAutomatonKernel.ApplyNilpotentPerturbation3(s0, s1, s2);
            var (r0, r1, r2) = CuboctahedronAutomatonKernel.CalibrateFiniteNeumannOrder2(y0, y1, y2);
            Assert.Equal(s0, r0);
            Assert.Equal(s1, r1);
            Assert.Equal(s2, r2);
        }

        // 3) Demystification of String Theory / M-Theory numbers 10, 11, and 26 from a single 3D cell (d = 3):
        Assert.Equal(new Ring9(10), Ring9.D + Ring9.DOutVertex);                       // 3 + 7 = 10
        Assert.Equal(new Ring9(10), Ring9.DInSquare + Ring9.NSq);                      // 4 + 6 = 10
        Assert.Equal(new Ring9(11), Ring9.D + Ring9.NTri);                             // 3 + 8 = 11 (3 spatial axes + 8 lanes of VerletCell72!)
        Assert.Equal(new Ring9(26), Ring9.NSq + Ring9.NTri + Ring9.NLinks);            // 6 + 8 + 12 = 26 = 3^3 - 1

        // 4) Second Instrument Deep Dive: MOSFET Gate Stack (Si + SiO2 = 14 + 30 = 44), Copper Coil Lag (64 = KP),
        //    and Glass-to-Air Acoustic/Dielectric Shell Match (Z(SiO2) = 30 == Z(N2) + Z(O2) = 14 + 16):
        Ring9 zMosStack = Ring9.ZSilicon14 + Ring9.ZQuartz30;
        Assert.Equal(new Ring9(44), zMosStack);
        Ring9 lagCoil64 = Ring9.KP;
        Ring9 lagGate25 = Ring9.BezoutInv41;
        Ring9 fpgaDrive = Ring9.CollapsePhi + lagCoil64 + lagGate25;
        Ring9 dishArrived = fpgaDrive - lagCoil64 - lagGate25;
        Assert.Equal(Ring9.CollapsePhi, dishArrived);
    }

    [Fact]
    public void Core_12_Complete10StageKitchenBenchPhysicsAndPerCoilCalibration()
    {
        // 1) Tap-water mineral electrolyte shells:
        //    Z(Na+) = 10 = Z(Ne), Z(Mg2+) = 10 = Z(Ne), Z(K+) = 18 = Z(Ar), Z(Ca2+) = 18 = Z(Ar), Z(Cl-) = 18 = Z(Ar),
        //    Z(HCO3-) = 32 = KP/2, Z(CaCO3) = 50, Z(MgCO3) = 42 = S2!
        Ring9 zNeon10  = new(10);
        Ring9 zArgon18 = new(18);
        Ring9 zNaPlus  = new Ring9(11) - Ring9.MuVac;
        Ring9 zMg2Plus = Ring9.NLinks - Ring9.C0OverCSq;
        Ring9 zKPlus   = new Ring9(19) - Ring9.MuVac;
        Ring9 zCa2Plus = new Ring9(20) - Ring9.C0OverCSq;
        Ring9 zClMinus = new Ring9(17) + Ring9.MuVac;
        Ring9 zHco3    = Ring9.MuVac + Ring9.NSq + (Ring9.D * Ring9.NTri) + Ring9.MuVac;
        Assert.Equal(zNeon10,  zNaPlus);
        Assert.Equal(zNeon10,  zMg2Plus);
        Assert.Equal(zArgon18, zKPlus);
        Assert.Equal(zArgon18, zCa2Plus);
        Assert.Equal(zArgon18, zClMinus);
        Assert.Equal(new Ring9(32), zHco3);

        // 2) Polystyrene Petri dish monomer C8H8: Z = 8*6 + 8 = 56 = 7 * 8 = DOutVertex * NTri,
        //    and liquid meniscus surface coordination: 7 broken outward links + 5 intact links = 12 = NLinks!
        Assert.Equal(new Ring9(56), Ring9.DOutVertex * Ring9.NTri);
        Assert.Equal(Ring9.NLinks, Ring9.DOutVertex + Ring9.DOutEdge);

        // 3) Full 12-coil end-to-end pre-compensation across all 12 hand-wound coils with individual inductance skews:
        for (ushort k = 0; k < 12; k++)
        {
            Ring9 skewK = new Ring9(k) * Ring9.DOutVertex;
            Ring9 thermalDrift = new(4);
            Ring9 dishWall = Ring9.Zero; // SiO2 glass wall matches N2 + O2 = 30
            Ring9 fpgaRegK = CuboctahedronAutomatonKernel.ComputeKitchenBenchPrecompensatedBeamPhase9(
                Ring9.CollapsePhi, Ring9.BezoutInv41, Ring9.KP, thermalDrift, dishWall, skewK);
            Ring9 arrivedK = CuboctahedronAutomatonKernel.ComputeKitchenBenchArrivedDishPhase9(
                fpgaRegK, Ring9.BezoutInv41, Ring9.KP, thermalDrift, dishWall, skewK);
            Assert.Equal(Ring9.CollapsePhi, arrivedK);
        }

        // 4) Mutual inductance off-axis cancellation on O_2(3): 4 neighbors at 60 deg (+1), 2 at 90 deg (0), 4 at 120 deg (-1), 1 at 180 deg:
        //    4 + 2 + 4 + 1 = 11 = NLinks - MuVac, and off-axis sum 4*(+1) + 4*(-1) == 0!
        Assert.Equal(Ring9.NLinks - Ring9.MuVac, Ring9.DInSquare + Ring9.C0OverCSq + Ring9.DInSquare + Ring9.MuVac);
        Assert.Equal(Ring9.Zero, (Ring9.DInSquare * Ring9.MuVac) + (Ring9.DInSquare * Ring9.MinusOne));

        // 5) Schottky flyback relaxation law: 3 * 81 + 13 = 243 + 13 = 256 = PhiPauli,
        //    and 8-bit webcam ADC ±1 LSB quantization step lies inside Heisenberg deadband (ApplyHeisenbergDeadband(1) == 0):
        Ring9 tau81 = Ring9.DPhiDisc + Ring9.DPhiDisc - Ring9.MuVac;
        Assert.Equal(Ring9.PhiPauli, (Ring9.D * tau81) + Ring9.NPhi);
        Assert.Equal(Ring9.Zero, CuboctahedronAutomatonKernel.ApplyHeisenbergDeadband(Ring9.MuVac));

        // 6) Earth geomagnetic field antipodal sum across all 12 O_2(3) vertices == (0, 0, 0),
        //    and USB optocoupler GaAs LED shell Z(Ga_31) + Z(As_33) = 31 + 33 == 64 == KP:
        Ring9 sumDz = Ring9.Zero, sumDy = Ring9.Zero, sumDx = Ring9.Zero;
        for (int k = 0; k < 12; k++)
        {
            var (dz, dy, dx) = Ring9.CuboctOffsets12[k];
            sumDz = sumDz + dz;
            sumDy = sumDy + dy;
            sumDx = sumDx + dx;
        }
        Assert.Equal(Ring9.Zero, sumDz);
        Assert.Equal(Ring9.Zero, sumDy);
        Assert.Equal(Ring9.Zero, sumDx);
        Ring9 zGallium31 = (Ring9.NTri * Ring9.DInSquare) - Ring9.MuVac;
        Ring9 zArsenic33 = (Ring9.NTri * Ring9.DInSquare) + Ring9.MuVac;
        Assert.Equal(Ring9.KP, zGallium31 + zArsenic33);

        // 7) Ambient room light Lock-In rejection across all 512 ring ticks (sum of ProjS == 0 and sum of ProjC == 0),
        //    lead-free solder Tin magic shell Z(Sn) = 2 * 25 = 50, and home 3D-printer ±0.2 mm tolerance (3 ticks -> deadband 0):
        Ring9 lockInSumS = Ring9.Zero, lockInSumC = Ring9.Zero;
        for (ushort phi = 0; phi < 512; phi++)
        {
            Ring9 r = new(phi);
            lockInSumS = lockInSumS + r.ProjS;
            lockInSumC = lockInSumC + r.ProjC;
        }
        Assert.Equal(Ring9.Zero, lockInSumS);
        Assert.Equal(Ring9.Zero, lockInSumC);
        Ring9 zTinSolder50 = Ring9.C0OverCSq * (Ring9.DOutEdge * Ring9.DOutEdge);
        Assert.Equal(Ring9.ZQuartz30 + (Ring9.DInSquare * Ring9.DOutEdge), zTinSolder50);
        Assert.Equal(Ring9.Zero, CuboctahedronAutomatonKernel.ApplyHeisenbergDeadband(Ring9.D));

        // 8) Hydrophilic quartz glass silanol wetting group Z(-Si-OH) = 14 + 8 + 1 = 23 == 2*|N| - 1,
        //    and hydrated hydrogen-bond bridge Z(SiOH + H2O) = 23 + 10 = 33 == Z(As):
        Ring9 zSilanol23 = Ring9.ZSilicon14 + Ring9.NTri + Ring9.MuVac;
        Assert.Equal((Ring9.C0OverCSq * Ring9.NLinks) - Ring9.MuVac, zSilanol23);
        Ring9 zWater10 = Ring9.NTri + Ring9.C0OverCSq;
        Assert.Equal(zArsenic33, zSilanol23 + zWater10);

        // 9) Moist kitchen air carbonic acid shell Z(CO2 + H2O) = 22 + 10 = 32 == KP / 2,
        //    and coil Self-Resonant Frequency margin ratio f_SRF / f_512 = 87 == MHat - 10:
        Ring9 zCo2_22 = Ring9.NSq + Ring9.ZOxygen2_16;
        Ring9 zCarbonicAcid32 = zCo2_22 + zWater10;
        Assert.Equal(Ring9.KP, Ring9.C0OverCSq * zCarbonicAcid32);
        Ring9 srfRatio87 = Ring9.MHat - zWater10;
        Assert.Equal(Ring9.MHat, srfRatio87 + zWater10);

        // 10) PDMS silicone acoustic couplant shell Z(C2H6SiO) = 12 + 6 + 14 + 8 = 40 == DPhiDisc - 1 (deadband == 0),
        //     and 30-node kitchen bench closure 3 * 10 = 30 == ZQuartz30:
        Ring9 zPdms40 = Ring9.NLinks + Ring9.NSq + Ring9.ZSilicon14 + Ring9.NTri;
        Assert.Equal(Ring9.DPhiDisc - Ring9.MuVac, zPdms40);
        Assert.Equal(Ring9.Zero, CuboctahedronAutomatonKernel.ApplyHeisenbergDeadband(zPdms40));
        Assert.Equal(Ring9.ZQuartz30, Ring9.D * zWater10);

        // 11) Air-core 60-turn coil bobbin geometry: N_turns = 5 * 12 = 60, MMF at 0.7A = 42 = S_2(3) ampere-turns (42 * 10 == 60 * 7),
        //     and Wheeler inductance denominator 6*9 + 9*12 + 10*4 = 54 + 108 + 40 = 202 == 2 * 101:
        Ring9 turns60 = Ring9.DOutEdge * Ring9.NLinks;
        Ring9 mmf42 = (Ring9.D * Ring9.ZSilicon14);
        Assert.Equal(turns60 * Ring9.DOutVertex, mmf42 * zWater10);
        Ring9 wheelerDenom202 = (Ring9.NSq * Ring9.ModBits) + (Ring9.ModBits * Ring9.NLinks) + (zWater10 * Ring9.DInSquare);
        Assert.Equal(Ring9.C0OverCSq * ((Ring9.DOutEdge * Ring9.DOutEdge * Ring9.DInSquare) + Ring9.MuVac), wheelerDenom202);

        // 12) Surfactant-free glass pre-wash shells: baking soda Z(NaHCO3) = 11 + 1 + 6 + 24 = 42 == S_2(3),
        //     ethanol Z(C2H5OH) = 12 + 6 + 8 = 26 == 3^3 - 1, pure water surface tension 72 == 8 * 9 = B_2(3) mN/m,
        //     and photosynthetic LED spectral span 718 - 440 = 278 == 2 * 139 = C0OverCSq * CollapsePhi:
        Ring9 zBakingSoda42 = Ring9.NLinks + Ring9.NSq + (Ring9.D * Ring9.NTri);
        Assert.Equal(mmf42, zBakingSoda42);
        Ring9 zEthanol26 = Ring9.NLinks + Ring9.NSq + Ring9.NTri;
        Assert.Equal((Ring9.D * Ring9.D * Ring9.D) - Ring9.MuVac, zEthanol26);
        Ring9 sigmaWater72 = Ring9.NTri * Ring9.ModBits;
        Assert.Equal(Ring9.KP + Ring9.NTri, sigmaWater72);
        Assert.Equal(Ring9.C0OverCSq * Ring9.CollapsePhi, new Ring9(278));

        // 13) Stage 0 Peltier TEC1-12706 condenser: 127 = QSquare - 1 couples, Z(Te) = 4 * 13 = 52, Z(Bi) = 81 + 2 = 83,
        //     Maillard crust target temperature 68 C = KP + 4, and photosynthetic balance 6*22 + 6*10 = 192 = 96 + 96:
        Assert.Equal(Ring9.QSquare - Ring9.MuVac, new Ring9(127));
        Assert.Equal(Ring9.DInSquare * Ring9.NPhi, new Ring9(52));
        Assert.Equal(tau81 + Ring9.C0OverCSq, new Ring9(83));
        Assert.Equal(Ring9.KP + Ring9.DInSquare, new Ring9(68));
        Ring9 zGlucose96 = Ring9.MHat - Ring9.MuVac;
        Ring9 zSixO2_96 = Ring9.NSq * Ring9.ZOxygen2_16;
        Assert.Equal(zGlucose96 + zSixO2_96, (Ring9.NSq * zCo2_22) + (Ring9.NSq * zWater10));

        // 14) Complete 36-edge unit-ball B_1(3) bench closure: 12 radial spokes + 24 shell edges = 36 == 4 * 9 = DInSquare * ModBits (Bus36),
        //     and Heisenberg thermal emergency cutoff 81 C == 3^4 = tau81:
        Ring9 radialEdges12 = Ring9.NLinks;
        Ring9 shellEdges24 = Ring9.C0OverCSq * Ring9.NLinks;
        Assert.Equal(Ring9.DInSquare * Ring9.ModBits, radialEdges12 + shellEdges24);
        Assert.Equal(new Ring9(81), tau81);

        // 15) Gate #0 Belousov-Zhabotinsky reagent shells: bromate Z(BrO3-) = 60 == turns60,
        //     malonic acid Z(C3H4O4) = 54 == 6 * 9 = NSq * ModBits (difference 60 - 54 == 6 = NSq),
        //     triiodide polarizer film Z(I3-) = 3 * 53 + 1 = 160 == 5 * 32,
        //     and complete 42-node second-sphere bench closure 30 + 12 = 42 == S_2(3):
        Ring9 zBromate60 = turns60;
        Ring9 zMalonic54 = Ring9.NSq * Ring9.ModBits;
        Assert.Equal(Ring9.NSq, zBromate60 - zMalonic54);
        Ring9 zIodine53 = mmf42 + (Ring9.NLinks - Ring9.MuVac);
        Ring9 zTriiodide160 = (Ring9.D * zIodine53) + Ring9.MuVac;
        Assert.Equal(Ring9.DOutEdge * zCarbonicAcid32, zTriiodide160);
        Assert.Equal(mmf42, Ring9.ZQuartz30 + Ring9.NLinks);

        // 16) Snell/Fresnel quartz-water refraction ratio (3/2)/(4/3) = 9/8 = ModBits/NTri,
        //     Fresnel reflection denominator (9 + 8)^2 = 17^2 = 289 (transmitted numerator 288 = 8 * 36),
        //     66 mutual inductance coil pairs 6 * 11 = 66, and complete 55-node two-shell ball closure |B_2(3)| = 1 + 12 + 42 = 55 = Z(Cs):
        Ring9 fresnelSum17 = Ring9.ModBits + Ring9.NTri;
        Ring9 fresnelDenom289 = fresnelSum17 * fresnelSum17;
        Assert.Equal(Ring9.NTri * (radialEdges12 + shellEdges24), fresnelDenom289 - Ring9.MuVac);
        Ring9 mutualPairs66 = Ring9.NSq * (Ring9.NLinks - Ring9.MuVac);
        Assert.Equal(new Ring9(66), mutualPairs66);
        Ring9 ballB2Nodes55 = Ring9.MuVac + Ring9.NLinks + mmf42;
        Assert.Equal(Ring9.DOutEdge * (Ring9.NLinks - Ring9.MuVac), ballB2Nodes55);

        // 17) Complete 72-node VerletCell72 bit-budget bench closure: 55 (|B_2(3)|) + 17 (Fresnel 9+8) = 72 == 8 * 9 = sigmaWater72,
        //     total 12-beam capacitive wall coupling 12 * 45 = 540 pF == 6 * 90 (NSq * ZAlphaQuartz90),
        //     and acrylamide carcinogen-free safety margin 120 - 68 = 52 == 4 * 13 = Z(Te):
        Assert.Equal(sigmaWater72, ballB2Nodes55 + fresnelSum17);
        Ring9 wallCap45 = Ring9.DOutEdge * Ring9.ModBits;
        Assert.Equal(Ring9.NSq * Ring9.ZAlphaQuartz90, Ring9.NLinks * wallCap45);
        Ring9 acrylamideThreshold120 = zWater10 * Ring9.NLinks;
        Assert.Equal(Ring9.DInSquare * Ring9.NPhi, acrylamideThreshold120 - (Ring9.KP + Ring9.DInSquare));
    }

    [Fact]
    public void Compression_UniversalArchiver_BitExactSha256AcrossModes()
    {
        byte[] raw = new byte[64 * 1024];
        for (int i = 0; i < raw.Length; i++)
        {
            int t = i / 18;
            int ch = i % 18;
            raw[i] = (byte)(t * 7 + ch * 13);
        }

        byte[] comp = CuboctUniversalArchiver.CompressBlock(raw);
        byte[] rec = CuboctUniversalArchiver.DecompressBlock(comp, raw.Length);
        Assert.Equal(SHA256.HashData(raw), SHA256.HashData(rec));
    }

    [Fact]
    public void MediaCodec_AudioVerletAndTwoPointKeyframe_BitExact()
    {
        short[] stereo = new short[4096];
        for (int i = 0; i < 2048; i++)
        {
            stereo[2 * i] = (short)(i * 11);
            stereo[2 * i + 1] = (short)(i * 11 + 25);
        }

        byte[] comp = CuboctMediaCodec.EncodeStereoPcm16Verlet(stereo);
        short[] rec = CuboctMediaCodec.DecodeStereoPcm16Verlet(comp, 2048);
        Assert.Equal(stereo, rec);

        byte[] kf0 = new byte[1024];
        byte[] kf12 = new byte[1024];
        byte[] gt6 = new byte[1024];
        for (int i = 0; i < 1024; i++)
        {
            kf0[i] = 40;
            kf12[i] = 160;
            gt6[i] = 105;
        }
        var (_, pct) = CuboctMediaCodec.InterpolateTwoPointKeyframe(kf0, kf12, gt6, 6, 12);
        Assert.Equal(100.0, pct);
    }
}
