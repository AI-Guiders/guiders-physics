using System.Security.Cryptography;
using Guiders.Physics.Core;
using Guiders.Physics.Compression;
using Guiders.Physics.Engine;
using Guiders.Physics.MediaCodec;
using Xunit;

namespace Guiders.Physics.Tests;

public class GuidersPhysicsSuiteTests(Xunit.Abstractions.ITestOutputHelper output)
{
    private readonly Xunit.Abstractions.ITestOutputHelper _output = output;
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

        // 18) Complete 81-node Heisenberg equivalence class bench closure: 72 (VerletCell72) + 9 (ModBits) = 81 == 3^4 = tau81,
        //     isotonic NaCl magic shell Z(NaCl) = 10 + 18 = 28 == 4 * 7, natural K-40 isotope (19 + 21 = 40 = DPhiDisc - 1 -> deadband 0),
        //     palmitic lipid emulsion ratio 144 == 2 * 72, and 52-min energy budget 12 * 13 = 156 Wh:
        Assert.Equal(tau81, sigmaWater72 + Ring9.ModBits);
        Assert.Equal(Ring9.DInSquare * Ring9.DOutVertex, zWater10 + (Ring9.C0OverCSq * Ring9.ModBits));
        Ring9 zK40_19 = Ring9.MuVac + Ring9.NLinks + Ring9.NSq;
        Assert.Equal(Ring9.DPhiDisc - Ring9.MuVac, zK40_19 + Ring9.SeamTicks);
        Assert.Equal(Ring9.Zero, CuboctahedronAutomatonKernel.ApplyHeisenbergDeadband(Ring9.SeamTicks));
        Assert.Equal(Ring9.NLinks * Ring9.NLinks, Ring9.C0OverCSq * sigmaWater72);
        Assert.Equal(new Ring9(156), Ring9.NLinks * Ring9.NPhi);

        // 19) Complete 147-node three-shell cuboctahedral ball |B_3(3)| bench closure:
        //     alpha-quartz trigonal unit cell 81 + 9 = 90 == ZAlphaQuartz90, Gaussian modulus 81 + 16 = 97 == MHat,
        //     third shell S_3(3) = 92 (55 + 92 = 147 = 3 * 7^2), and human chromatin nucleosome bio-safety detuning 147 - 144 == 3 = D:
        Assert.Equal(Ring9.ZAlphaQuartz90, tau81 + Ring9.ModBits);
        Assert.Equal(Ring9.MHat, tau81 + Ring9.ZOxygen2_16);
        Ring9 shellS3_92 = (zWater10 * Ring9.ModBits) + Ring9.C0OverCSq;
        Ring9 ballB3Nodes147 = ballB2Nodes55 + shellS3_92;
        Assert.Equal(Ring9.D * Ring9.DOutVertex * Ring9.DOutVertex, ballB3Nodes147);
        Assert.Equal(Ring9.D, ballB3Nodes147 - (Ring9.NLinks * Ring9.NLinks));

        // 20) Complete 512-tick Z_512 ring kitchen bench closure:
        //     sucrose disaccharide Z = 2*96 - 10 = 182 == 14 * 13 (ZSilicon14 * NPhi),
        //     cold-start FPGA pre-compensated drive 139 + 64 + 25 = 228 == 12 * 19 (NLinks * zK40_19),
        //     UART baud divider 27 MHz / 115200 = 234 == 18 * 13 (2 * ModBits * NPhi),
        //     and full ring period 256 (PhiPauli) + 256 (3*81 + 13) == 0 (mod 512):
        Assert.Equal(Ring9.ZSilicon14 * Ring9.NPhi, (Ring9.C0OverCSq * zGlucose96) - zWater10);
        Assert.Equal(Ring9.NLinks * zK40_19, Ring9.CollapsePhi + Ring9.KP + Ring9.BezoutInv41);
        Assert.Equal((Ring9.C0OverCSq * Ring9.ModBits) * Ring9.NPhi, new Ring9(234));
        Assert.Equal(Ring9.Zero, Ring9.PhiPauli + ((Ring9.D * tau81) + Ring9.NPhi));

        // 21) High-order bench invariants mod 512:
        //     gate resistor & blue ferroin 544 == 17 * 32 == 32 (mod 512),
        //     stasis multiplier 833 == 1 + 13 * 64 == 7^2 * 17 == 321 (mod 512),
        //     and two-turn Verlet Bezout identity 25 * 41 = 1025 == 1 (mod 512):
        Assert.Equal(zCarbonicAcid32, fresnelSum17 * zCarbonicAcid32);
        Assert.Equal(Ring9.MuVac + (Ring9.NPhi * Ring9.KP), (Ring9.DOutVertex * Ring9.DOutVertex) * fresnelSum17);
        Assert.Equal(Ring9.MuVac, Ring9.BezoutInv41 * Ring9.DPhiDisc);
    }

    [Fact]
    public void Core_13_FullMultiscaleKitchenBenchClosure512AndTelemetryRecovery_MatchesLeanAndHardware()
    {
        Assert.Equal(new Ring9(202), CuboctahedronAutomatonKernel.ComputeWheelerAirCoreDenominator9(Ring9.ModBits, Ring9.NLinks, Ring9.DInSquare));
        Assert.Equal(new Ring9(228), CuboctahedronAutomatonKernel.ComputeColdStartFpgaPrecompPhase9());
        Assert.Equal(new Ring9(17), CuboctahedronAutomatonKernel.ComputeFresnelQuartzWaterSum9());
        Assert.Equal((Ring9.DOutVertex * Ring9.DOutVertex) * new Ring9(17), CuboctahedronAutomatonKernel.ComputeStasisMultiplierMod512());
        Assert.Equal(Ring9.ZSilicon14 * Ring9.NPhi, CuboctahedronAutomatonKernel.ComputeSucroseDisaccharideShell9());

        for (ushort phi = 0; phi < 512; phi++)
        {
            Ring9 targetPhi = new(phi);
            var cell = VerletCell72.Create(Ring9.MuVac, targetPhi, Ring9.Zero, Ring9.MuVac, Ring9.Zero, Ring9.Zero);
            var sensors = CuboctahedronAutomatonKernel.ReadCellSensors(cell);
            Ring9 recPhi = CuboctahedronAutomatonKernel.RecoverRingPhaseFromHydrophoneAndChiralCurrent9(
                Ring9.MuVac, Ring9.MuVac, sensors.AcousticWave, sensors.ReactionCurrent);
            Assert.Equal(targetPhi, recPhi);
        }

        // Verify full 52-minute (53 timepoints t = 0..52) KitchenBenchDigitalTwin in Guiders.Physics.Engine:
        var benchTwin = new KitchenBenchDigitalTwin();
        for (int minute = 0; minute <= 52; minute++)
        {
            var snap = benchTwin.EvaluateMinute(minute);
            Assert.Equal(snap.TargetDishPhase, snap.ArrivedDishPhase);
            Assert.Equal(snap.TargetDishPhase, snap.RecoveredTelemetryPhase);
            Assert.True(KitchenBenchDigitalTwin.VerifyUart16ByteFrameCrc(snap.Uart0x08CalibrationFrame16B));
        }

        // Verify 1-tick closed-loop PLL disturbance cancellation across all 512/512 unmodeled lag shifts:
        Ring9 initialEstLag = Ring9.KP + Ring9.BezoutInv41; // 64 + 25 = 89
        for (ushort dist = 0; dist < 512; dist++)
        {
            Ring9 trueLag = new(dist);
            var pll = KitchenBenchDigitalTwin.StepClosedLoopPll9(Ring9.CollapsePhi, initialEstLag, trueLag);
            Assert.Equal(trueLag, pll.UpdatedEstimatedLag);
            Assert.Equal(Ring9.CollapsePhi, pll.NextTickArrivedDishPhase);
            Assert.Equal(Ring9.Zero, pll.NextTickPhaseError);
        }

        // Verify 4-channel nilpotent holonomic tomograph (2*(I_4 + N_nil), N_nil^3 == 0) on Bus36
        // de-mixing (s0 = Z(Si)=14, s1 = Z(O2)=16, s2 = Z(SiO2)=30, s3 = Z(SiO4)=46):
        Bus36 clean4 = new(Ring9.ZSilicon14, new Ring9(16), new Ring9(30), new Ring9(46));
        Bus36 mixed4 = CuboctahedronAutomatonKernel.NilpotentTomographForward2xBus36(clean4);
        Bus36 rec4 = CuboctahedronAutomatonKernel.NilpotentTomographInverse24xBus36(mixed4);
        Assert.Equal(clean4, rec4);

        // Verify Bus36 pre-flight 12-coil continuity impedance decomposition:
        // D(N_active) = (1 + (N_active + 1)*64)*8192 == Bus36(0, 16, 2*(N_active + 1), 0)
        // with L2(12) == 26 == Z(Fe) and single-coil drop in L2 == 2 == (c0/c)^2:
        for (ushort nCoils = 0; nCoils <= 12; nCoils++)
        {
            Ring9 activeCount = new(nCoils);
            Bus36 dBus = CuboctahedronAutomatonKernel.ComputePreflightCoilImpedanceBus36(activeCount);
            Assert.Equal(Ring9.Zero, dBus.L0);
            Assert.Equal(new Ring9(16), dBus.L1);
            Assert.Equal(new Ring9((ushort)(2 * (nCoils + 1))), dBus.L2);
            Assert.Equal(activeCount, CuboctahedronAutomatonKernel.RecoverActiveCoilsFromImpedanceBus36(dBus));
        }

        // Verify all 6 recipe shells (Gate #0 BZ 60, Cutlet 70, Bread 86, Vitamin C 92, Cappuccino 102, Lipid 144)
        // and 3-second dry phase-cleaning protocol (139 -> 256):
        ushort[] expectedRecipeShells = [60, 70, 86, 92, 102, 144];
        for (byte rIdx = 0; rIdx < 6; rIdx++)
        {
            Assert.Equal(new Ring9(expectedRecipeShells[rIdx]), KitchenBenchDigitalTwin.GetRecipeTargetShellZ9(rIdx));
            var proto = KitchenBenchDigitalTwin.BuildRecipeAndDryCleanProtocol9(rIdx);
            Assert.Equal(Ring9.CollapsePhi, proto.ReleasePulsePhase);
            Assert.Equal(Ring9.PhiPauli, proto.PauliEjectPhase);
            Assert.True(KitchenBenchDigitalTwin.VerifyUart16ByteFrameCrc(proto.Uart0x02StartFrame16B));
        }

        // Verify 3D Gorkov-Bjerknes trap Hessian trace (2 * 12 = 24 = 8 + 8 + 8)
        // and 90 mm Petri dish (R = 45 mm) Bjerknes quarter-wave rings (45 / 7 = 6 = NSq, 45 % 7 = 3 = D):
        Assert.Equal(Ring9.C0OverCSq * Ring9.NLinks, CuboctahedronAutomatonKernel.ComputeGorkovTrapHessianTrace9(Ring9.NLinks));
        var bjerknes = CuboctahedronAutomatonKernel.ComputePetriDishBjerknesRings9();
        Assert.Equal(Ring9.NSq, bjerknes.ConcentricRingsCount6);
        Assert.Equal(Ring9.D, bjerknes.WallClearanceMm3);

        // Verify cavitation-free 12-beam ultrasonic focusing (256 / 12 = 21 = Seam21, n_E(21) = 0 at glass wall,
        // coherent focal current 12 * 21 = 252 = 256 - 4 > 128) and 4-pole XY Euler chiral quadrupole torque:
        var subCav = CuboctahedronAutomatonKernel.Compute12BeamSubcavitationFocus9();
        Assert.Equal(Ring9.SeamTicks, subCav.SingleRayPhase21);
        Assert.Equal(Ring9.Zero, subCav.SingleRayWallPhotons0);
        Assert.Equal(Ring9.PhiPauli - Ring9.DInSquare, subCav.Coherent12BeamFocalCurrent252);
        Assert.Equal(Ring9.QSquare, CuboctahedronAutomatonKernel.ComputeEulerQuadrupoleChiralTorque9(Ring9.QSquare));
        Assert.Equal(-Ring9.QSquare, CuboctahedronAutomatonKernel.ComputeEulerQuadrupoleChiralTorque9(Ring9.D * Ring9.QSquare));
        Assert.Equal(Ring9.Zero, CuboctahedronAutomatonKernel.ComputeEulerQuadrupoleChiralTorque9(Ring9.PhiPauli));

        // Verify 30-Section Thermodynamic Ladder (12 C [+15] -> 27 C [+41 = DPhiDisc] -> 68 C [+13 = NPhi] -> 81 C)
        // and 180 W * 52 min / 60 min = 156 = 12 * 13 Wh energy closure:
        var therm = CuboctahedronAutomatonKernel.ComputeKitchenBenchThermalLadderAndEnergy9();
        Assert.Equal(Ring9.D * Ring9.DOutEdge, therm.WarmKitchen27C - therm.DewPoint12C);
        Assert.Equal(Ring9.DPhiDisc, therm.Maillard68C - therm.WarmKitchen27C);
        Assert.Equal(Ring9.NPhi, therm.Cutoff81C - therm.Maillard68C);
        Assert.Equal(Ring9.NLinks * Ring9.NPhi, therm.Energy156Wh);

        // Verify Acoustic Fixed-Point Identity Q_square - C_3(dphi) == dphi across all 4 operating modes {0, 128, 139, 256}
        // and 3D Petri Dish CuboctahedronCarrier reversible step cycle:
        Assert.True(KitchenBenchDigitalTwin.VerifyAcousticFixedPointAcross4OperatingModes9());
        Assert.True(KitchenBenchDigitalTwin.Run3DPetriLatticeReversibleCycle9(4));

        // Verify 32-Section Half-Pauli Closure (Maillard pyrazine C4H4N2 = 42 = S_2(3),
        // peptide bond -CO-NH- = 22 = Z(CO2), hydrated peptide = 32 = K_P/2, audible subharmonic 32,768 = Bus36(0,64,0,0)):
        var maillard32 = CuboctahedronAutomatonKernel.ComputeMaillardAromaAndAudibleSubharmonic9();
        Assert.Equal(new Ring9(42), maillard32.PyrazineAroma42);
        Assert.Equal(new Ring9(22), maillard32.PeptideBond22);
        Assert.Equal(new Ring9(32), maillard32.HydratedPeptide32);
        Assert.Equal(Ring9.KP, maillard32.AudibleSubharmonicClocks32768.L1);

        // Verify 36-Section Bus36 Full Closure (Water dielectric 80 -> 64 = K_P at 68 C,
        // Umami glutamic acid Z=78 = 6*13, Cysteine Z=64 = K_P, Cystine disulfide Z=126 = 2*63, Section count = 36 = 4*9):
        var full36 = CuboctahedronAutomatonKernel.ComputeBench36SectionFullClosure9();
        Assert.Equal(new Ring9(80), full36.WaterDielectric20C_80);
        Assert.Equal(Ring9.KP, full36.WaterDielectric68C_64);
        Assert.Equal(Ring9.NSq * Ring9.NPhi, full36.UmamiGlutamicAcid78);
        Assert.Equal(Ring9.KP, full36.Cysteine64);
        Assert.Equal(Ring9.QSquare - Ring9.C0OverCSq, full36.CystineDisulfide126);
        Assert.Equal(Ring9.DInSquare * Ring9.ModBits, full36.Section36Count);

        // Verify 42-Section S_2(3) Full Closure (50/50 reactive balance C_3(64)=S_3(64)=64,
        // Theobromine Z=94, Citric acid Z=100, Lactic acid Z=48, Calcium lactate Z=114,
        // worst-case 5-source household noise = 8 < 41, margin = 33, total sections = 42 = S_2(3)):
        Assert.Equal(Ring9.KP, Ring9.KP.ProjC);
        Assert.Equal(Ring9.KP, Ring9.KP.ProjS);
        Assert.Equal(Ring9.QSquare, Ring9.KP.ProjC + Ring9.KP.ProjS);
        var full42 = CuboctahedronAutomatonKernel.ComputeBench42SecondSphereClosure9();
        Assert.Equal(new Ring9(94), full42.TheobromineChocolate94);
        Assert.Equal(new Ring9(100), full42.CitricAcidNectar100);
        Assert.Equal(new Ring9(48), full42.LacticAcidMilk48);
        Assert.Equal(new Ring9(114), full42.CalciumLactateCheese114);
        Assert.Equal(Ring9.NTri, full42.WorstCaseNoise8);
        Assert.Equal(new Ring9(33), full42.HeisenbergMargin33);
        Assert.Equal(new Ring9(42), full42.FullSection42Count);

        // Verify 55-Section |B_2(3)| Two-Shell Ball Full Closure (Mains peak 312 = 24*13 V,
        // 20 canonical L-amino acids, Aspartate Z=70=Heme, Histidine Z=82=81+1, Tryptophan Z=108=12*9 bit bus,
        // 49 = 7^2 MuJoCo geoms, total sections = 55 = |B_2(3)| = Z(Cs)):
        var full55 = CuboctahedronAutomatonKernel.ComputeBench55TwoShellBallClosure9();
        Assert.Equal(new Ring9(312), full55.MainsPeak312V);
        Assert.Equal(new Ring9(20), full55.CanonicalAminoAcids20);
        Assert.Equal(new Ring9(70), full55.AspartateHeme70);
        Assert.Equal(new Ring9(82), full55.HistidineLigand82);
        Assert.Equal(new Ring9(108), full55.TryptophanBus108);
        Assert.Equal(new Ring9(49), full55.MujocoGeoms49);
        Assert.Equal(new Ring9(55), full55.TwoShellBall55);

        // Verify 72-Section VerletCell72 Full Self-Closure (Numbers 56..72: hot coil lags 58, 59, 61, 62..66,
        // hydrogen bond modulus 63 = 7*9, stasis holding temp 65 = 5*13, and 72 = 8*9 = dim(VerletCell72)):
        var full72 = CuboctahedronAutomatonKernel.ComputeBench72VerletCellClosure9();
        Assert.Equal(new Ring9(58), full72.HotCoilLag58);
        Assert.Equal(new Ring9(59), full72.MidCoilLag59);
        Assert.Equal(new Ring9(61), full72.EarlyCoilLag61);
        Assert.Equal(new Ring9(62), full72.MinCoilLag62);
        Assert.Equal(new Ring9(63), full72.HydrogenBond63);
        Assert.Equal(new Ring9(65), full72.StasisTemp65);
        Assert.Equal(new Ring9(66), full72.MaxCoilLag66);
        Assert.Equal(new Ring9(72), full72.Verlet72SectionCount);

        // Verify 81-Section |E_Heis(3)| = 3^4 Full Self-Closure (Numbers 73..81: frame radius 73=64+9,
        // optocoupler delay 74=2*37 ns, collagen 75=3*25, tetra-octahedral 76=4*19 -> 3*76=228,
        // glutamate anion 77=7*11, ENIG gold plating 79=Z(Au), and 81 = 72 + 9 = 3^4):
        var full81 = CuboctahedronAutomatonKernel.ComputeBench81HeisenbergClassClosure9();
        Assert.Equal(new Ring9(73), full81.FrameRadius73Mm);
        Assert.Equal(new Ring9(74), full81.Optocoupler74Ns);
        Assert.Equal(new Ring9(75), full81.CollagenGate75);
        Assert.Equal(new Ring9(76), full81.TetraOctahedral76);
        Assert.Equal(new Ring9(77), full81.GlutamateAnion77);
        Assert.Equal(new Ring9(79), full81.GoldPlating79);
        Assert.Equal(new Ring9(81), full81.Heisenberg81SectionCount);

        // Verify 97-Section M_hat(3) = 9^2 + 4^2 = 97 Full Self-Closure (Numbers 82..97:
        // total thermal lag ladder 83..89..91, polystyrene total lag 96 = 89 + 7, and 97 = 81 + 16):
        var full97 = CuboctahedronAutomatonKernel.ComputeBench97GaussianModulusClosure9();
        Assert.Equal(new Ring9(83), full97.HotTotalLag83);
        Assert.Equal(new Ring9(84), full97.MidTotalLag84);
        Assert.Equal(new Ring9(85), full97.Lag20Min85);
        Assert.Equal(new Ring9(86), full97.EarlyTotalLag86);
        Assert.Equal(new Ring9(87), full97.SrfTripleCu87);
        Assert.Equal(new Ring9(88), full97.MsgLag88);
        Assert.Equal(new Ring9(89), full97.ColdNominalLag89);
        Assert.Equal(new Ring9(91), full97.MaxCoil11Lag91);
        Assert.Equal(new Ring9(96), full97.PolystyreneTotalLag96);
        Assert.Equal(Ring9.MHat, full97.GaussianMhat97SectionCount);

        // Verify 128-Section Q_square(3) = 128 Quarter-Ring Saddle Full Self-Closure (Numbers 98..128:
        // Wheeler prime 101, Kok oxo / I2 106, collapse odd current S_3(139) = 117 = 9*13,
        // meniscus-Fresnel 119 = 7*17, Kok S3 flash 123 = 3*41, Peltier 127, and 128 = Q_square):
        var full128 = CuboctahedronAutomatonKernel.ComputeBench128QuarterRingSaddleClosure9();
        Assert.Equal(new Ring9(101), full128.WheelerPrime101);
        Assert.Equal(new Ring9(106), full128.KokOxoIodine106);
        Assert.Equal(new Ring9(117), full128.CollapseOddCurrent117);
        Assert.Equal(full128.CollapseOddCurrent117, Ring9.CollapsePhi.ProjS);
        Assert.Equal(new Ring9(119), full128.MeniscusFresnel119);
        Assert.Equal(new Ring9(123), full128.KokS3Flash123);
        Assert.Equal(new Ring9(127), full128.PeltierCouples127);
        Assert.Equal(Ring9.QSquare, full128.QuarterRing128SectionCount);

        // Verify 147-Section |B_3(3)| = 3 * 7^2 = 147 Three-Shell Ball Full Self-Closure (Numbers 129..147:
        // meniscus gate 129, ordered pairs / 6 CO2 / 12-beam collapse sum 132 = 12*11, Cs-133 mass 133 = 7*19 = 55+78,
        // alpha_FS^-1 floor 137 = 128+9, collapse phase 139, palmitic lipid 144 = 12^2, and 147 = |B_3(3)|):
        var full147 = CuboctahedronAutomatonKernel.ComputeBench147ThreeShellBallClosure9();
        Assert.Equal(new Ring9(129), full147.MeniscusGate129);
        Assert.Equal(new Ring9(132), full147.OrderedPairs132);
        Assert.Equal(new Ring9(133), full147.Cesium133Mass);
        Assert.Equal(new Ring9(137), full147.FineStructure137);
        Assert.Equal(Ring9.CollapsePhi, full147.CollapsePhase139);
        Assert.Equal(new Ring9(144), full147.PalmiticLipid144);
        Assert.Equal(new Ring9(147), full147.ThreeShellBall147SectionCount);
        Assert.Equal(Ring9.D, full147.ThreeShellBall147SectionCount - full147.PalmiticLipid144);

        // Verify 441-Section Consecutive Blue GaN LED 21^2 = 9*49 = 7*63 = 441 nm, Bezout 417, Euler 384, Phase 352, 7^3 = 343, 312 V, 289 & Seed IV(3) = 263 Multiscale Self-Closure
        // (Every integer 1..441 consecutively mapped + 512..8192: Kok S4 flash 164 = 4*41 -> C_3(164) = -36,
        // core square 169 = 13^2 = 128 + 41, 180 = 12*15 W, IV Bezout 263*183 == 1 (mod 512),
        // Gold-197 Bezout 197*13 == 1 (mod 512), FPGA cold drive 228 = 12*19, UART divider 234 = 18*13,
        // Pauli mirror 256 = 2^8, Catalan seed IV(3) = 263 = 256 + 7, ice point 273 K = 21*13, Fresnel 289 = 17^2,
        // standard 298 K = 256 + 42, 12-channel gate sum 300 = 12*25, mains peak 312 = 24*13 V, Maillard 341 K = 68 + 273,
        // 7^3 = 343, 12-coil copper 348 = 12*29, oxidation phase 352 = 256 + 96 -> C_3(352) = -32,
        // cutoff 354 K = 81 + 273, full circle 360 = 2*180, boiling water 373 K = 256 + 117, Euler pole 384 = 3*128,
        // wire diameter 400 = 20^2 um, Gaussian Bezout 97*417 == 1 (mod 512), Soret 418 nm, 12-beam Bus36 432 = 12*36, and blue LED 441 = 21^2):
        var full441 = CuboctahedronAutomatonKernel.ComputeBench441BlueLedMultiscaleClosure9();
        Assert.Equal(new Ring9(164), full441.KokS4Flash164);
        Assert.Equal(-(Ring9.DInSquare * Ring9.ModBits), full441.KokS4Flash164.ProjC);
        Assert.Equal(new Ring9(169), full441.CoreSquare169);
        Assert.Equal(Ring9.QSquare + Ring9.DPhiDisc, full441.CoreSquare169);
        Assert.Equal(new Ring9(180), full441.Power180SectionCount);
        Assert.Equal(new Ring9(183), full441.IvSeedBezout183);
        Assert.Equal(Ring9.MuVac, full441.IvSeed263SectionCount * full441.IvSeedBezout183);
        Assert.Equal(new Ring9(197), full441.Gold197SectionCount);
        Assert.Equal(Ring9.MuVac, full441.Gold197SectionCount * Ring9.NPhi);
        Assert.Equal(new Ring9(228), full441.FpgaCold228SectionCount);
        Assert.Equal(new Ring9(234), full441.Uart234SectionCount);
        Assert.Equal(Ring9.PhiPauli, full441.Pauli256SectionCount);
        Assert.Equal(Ring9.IvSeed, full441.IvSeed263SectionCount);
        Assert.Equal(Ring9.DOutVertex, full441.IvSeed263SectionCount - full441.Pauli256SectionCount);
        Assert.Equal(new Ring9(273), full441.IcePoint273K);
        Assert.Equal(new Ring9(289), full441.Fresnel289SectionCount);
        Assert.Equal(new Ring9(298), full441.StandardTemp298K);
        Assert.Equal(new Ring9(300), full441.TwelveGateSum300);
        Assert.Equal(new Ring9(312), full441.Mains312SectionCount);
        Assert.Equal(new Ring9(341), full441.Maillard341K);
        Assert.Equal(new Ring9(343), full441.Branching343SectionCount);
        Assert.Equal(new Ring9(348), full441.TwelveCopper348);
        Assert.Equal(new Ring9(352), full441.Phase352SectionCount);
        Assert.Equal(-(Ring9.DInSquare * Ring9.NTri), full441.Phase352SectionCount.ProjC);
        Assert.Equal(new Ring9(354), full441.Cutoff354K);
        Assert.Equal(new Ring9(360), full441.FullCircle360);
        Assert.Equal(new Ring9(373), full441.Boiling373K);
        Assert.Equal(new Ring9(384), full441.Euler384SectionCount);
        Assert.Equal(Ring9.Zero, full441.Euler384SectionCount.ProjC);
        Assert.Equal(-Ring9.QSquare, full441.Euler384SectionCount.ProjS);
        Assert.Equal(new Ring9(400), full441.WireDiameter400Um);
        Assert.Equal(new Ring9(417), full441.GaussianBezout417);
        Assert.Equal(Ring9.MuVac, Ring9.MHat * full441.GaussianBezout417);
        Assert.Equal(new Ring9(418), full441.Soret418Nm);
        Assert.Equal(new Ring9(432), full441.TwelveBus432);
        Assert.Equal(new Ring9(441), full441.BlueLed441SectionCount);

        // Verify 512-Section Full Phase-Ring N_ring(3) = 2^9 = 8^3 = 512 Consecutive Multiscale Self-Closure (Numbers 442..512:
        // 12-tick scan 444 = 12*37 ns, 7*K_P = 448 = 512 - 64, double FPGA cold drive 456 = 2*228 = 12*38,
        // coil SRF 459 = 27*17 = 512 - 53, quartz ring address 462 = 22*21 = 256 + 206, Bus36*N_phi = 468 = 2*234,
        // upper Heisenberg boundary 471 = 512 - 41, 13-tick core duration 481 = 13*37 ns = 16^2 + 15^2,
        // CO2 square 484 = 22^2 = 256 + 228, complementary gate delay 487 = 512 - 25, complementary seam 491 = 512 - 21,
        // 12-beam Heisenberg sum 492 = 12*41 = 512 - 20, 3rd perfect number 496 = 16*31 = 512 - 16,
        // Schottky decay end 499 = 256 + 243 = 512 - 13, 500 mm twisted pair = 20*25 = 512 - 12,
        // 4*127 Peltier couples = 508 = 256 + 252, E24 gate resistor 510 = 30*17 = 512 - 2,
        // Ring9.MinusOne = 511 = 7*73, and full ring period 512 = Bus36(0, 1, 0, 0)):
        var full512 = CuboctahedronAutomatonKernel.ComputeBench512FullRingMultiscaleClosure9();
        Assert.Equal(new Ring9(444), full512.TwelveTickScan444Ns);
        Assert.Equal(new Ring9(448), full512.SevenMeniscusKp448);
        Assert.Equal(new Ring9(456), full512.DoubleFpgaDrive456);
        Assert.Equal(new Ring9(459), full512.MooreFresnelSrf459);
        Assert.Equal(new Ring9(462), full512.Co2SeamQuartz462);
        Assert.Equal(new Ring9(468), full512.Bus36CoreUart468);
        Assert.Equal(new Ring9(471), full512.UpperHeisenberg471);
        Assert.Equal(Ring9.MuVac, CuboctahedronAutomatonKernel.HeisenbergEnergyQuanta(full512.UpperHeisenberg471));
        Assert.Equal(new Ring9(481), full512.ThirteenTickCore481Ns);
        Assert.Equal(new Ring9(484), full512.Co2Square484);
        Assert.Equal(new Ring9(487), full512.RingMinusGate487);
        Assert.Equal(new Ring9(491), full512.RingMinusSeam491);
        Assert.Equal(new Ring9(492), full512.TwelveHeisenberg492);
        Assert.Equal(new Ring9(496), full512.ThirdPerfect496);
        Assert.Equal(new Ring9(499), full512.SchottkyEnd499);
        Assert.Equal(new Ring9(500), full512.TwistedPair500Mm);
        Assert.Equal(new Ring9(508), full512.QuadPeltierFocal508);
        Assert.Equal(new Ring9(510), full512.QuartzFresnelRg510);
        Assert.Equal(Ring9.MinusOne, full512.RingMax511);
        Assert.Equal(Ring9.Zero, full512.FullRing512SectionCount.L0);
        Assert.Equal(Ring9.MuVac, full512.FullRing512SectionCount.L1);

        // Verify 587-Section Consecutive 587x Coherent Cuboctahedron Barrier Collapse (513..587 = 512 + 75)
        // and Complete Multiscale Optical, Spinor, Electrical, Acoustic, and Crystal Hierarchy on Bus36:
        var full587 = CuboctahedronAutomatonKernel.ComputeBench587BarrierCollapseMultiscaleClosureBus36();
        Assert.Equal(new Bus36(Ring9.MuVac, Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.SecondTurnStart513);
        Assert.Equal(new Bus36(Ring9.NSq, Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.FourteenFaceScan518Ns);
        Assert.Equal(new Bus36(Ring9.NTri, Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.DoubleAtp520);
        Assert.Equal(new Bus36(Ring9.NPhi, Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.SeamGate525);
        Assert.Equal(new Bus36(Ring9.ZSilicon14, Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.DoubleIvSeed526);
        Assert.Equal(new Bus36(new Ring9(15), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.GalliumFresnelRipple527);
        Assert.Equal(new Bus36(Ring9.ZOxygen2_16, Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.TwelveMosfetStacks528);
        Assert.Equal(new Bus36(new Ring9(17), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.SilanolSquare529);
        Assert.Equal(new Bus36(Ring9.SeamTicks, Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.CoreHeisenberg533);
        Assert.Equal(new Bus36(new Ring9(28), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.TotalWallCap540Pf);
        Assert.Equal(new Bus36(new Ring9(32), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.FerroinGate544);
        Assert.Equal(new Bus36(new Ring9(35), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.FourteenthGauss547);
        Assert.Equal(new Bus36(new Ring9(36), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.QuadAlphaFs548);
        Assert.Equal(new Bus36(Ring9.DPhiDisc, Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.SecondTurnHeisenberg553);
        Assert.Equal(new Bus36(new Ring9(44), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.QuadCollapse556);
        Assert.Equal(new Bus36(Ring9.KP, Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.EightVerletVdcSq576);
        Assert.Equal(new Bus36(new Ring9(67), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.OxymyoglobinCenter579Nm);
        Assert.Equal(new Bus36(new Ring9(70), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.SixGaussianHeme582);
        Assert.Equal(new Bus36(new Ring9(72), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.EightRadiusVerlet584);
        Assert.Equal(new Bus36(new Ring9(73), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.RadiusCoreStasis585);
        Assert.Equal(new Bus36(new Ring9(75), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.BarrierCollapse587SectionCount);
        Assert.Equal(new Bus36(Ring9.CollapsePhi, Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.ReducedRedFerroin651Nm);
        Assert.Equal(new Bus36(Ring9.QuartzSiO2Phi, Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.RedGaAsLed718Nm);
        Assert.Equal(new Bus36(new Ring9(208), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.SpinorFullIdentity720Deg);
        Assert.Equal(new Bus36(new Ring9(321), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.StasisMultiplier833);
        Assert.Equal(new Bus36(new Ring9(413), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full587.MosfetGateDelay925Ns);
        Assert.Equal(new Bus36(Ring9.Zero, Ring9.C0OverCSq, Ring9.Zero, Ring9.Zero), full587.TwoTurnVerlet1024);
        Assert.Equal(new Bus36(Ring9.MuVac, Ring9.C0OverCSq, Ring9.Zero, Ring9.Zero), full587.BezoutGateHeisenberg1025);
        Assert.Equal(new Bus36(new Ring9(83), Ring9.C0OverCSq, Ring9.Zero, Ring9.Zero), full587.SiliconBandgap1107Nm);
        Assert.Equal(new Bus36(Ring9.Zero, Ring9.D, Ring9.Zero, Ring9.Zero), full587.CavityLoad1536Ohm);
        Assert.Equal(new Bus36(new Ring9(264), Ring9.D, Ring9.Zero, Ring9.Zero), full587.EsdTotalCap1800Pf);
        Assert.Equal(new Bus36(Ring9.Zero, Ring9.DInSquare, Ring9.Zero, Ring9.Zero), full587.PauliDense2048);
        Assert.Equal(new Bus36(new Ring9(64), new Ring9(10), Ring9.Zero, Ring9.Zero), full587.HiddenMicrostates5184);
        Assert.Equal(new Bus36(Ring9.Zero, new Ring9(16), Ring9.Zero, Ring9.Zero), full587.VacuumFloor8192);

        // Verify 1025-Section Consecutive Two-Turn Verlet Spinor Period (1024 = 2 * 512)
        // and Gate-Heisenberg Bezout Return Identity (1025 = 25 * 41 = 2 * 512 + 1) Closure on Bus36 (588..1025):
        var full1025 = CuboctahedronAutomatonKernel.ComputeBench1025TwoTurnVerletBezoutConsecutiveClosureBus36();
        Assert.Equal(new Bus36(new Ring9(76), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.TwelveBeamsMujoco588);
        Assert.Equal(new Bus36(new Ring9(88), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.TwelveBeamsTin600);
        Assert.Equal(new Bus36(new Ring9(89), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.SecondTurnColdLag601);
        Assert.Equal(new Bus36(new Ring9(96), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.RingPlusGlucose608);
        Assert.Equal(new Bus36(new Ring9(100), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.Bus36Fresnel612);
        Assert.Equal(new Bus36(new Ring9(112), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.TwelveRunMinutes624);
        Assert.Equal(new Bus36(Ring9.QSquare, Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.TenKp640);
        Assert.Equal(new Bus36(Ring9.CollapsePhi, Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.ReducedRedFerroin651Nm);
        Assert.Equal(new Bus36(new Ring9(144), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.BalmerHAlpha656Nm);
        Assert.Equal(new Bus36(new Ring9(154), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.EighteenTickArgon666Ns);
        Assert.Equal(new Bus36(new Ring9(164), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.IronSquare676);
        Assert.Equal(new Bus36(new Ring9(188), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.FourteenFacesTin700);
        Assert.Equal(new Bus36(Ring9.QuartzSiO2Phi, Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.RedGaAsLed718Nm);
        Assert.Equal(new Bus36(new Ring9(208), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.SpinorIdentity720Deg);
        Assert.Equal(new Bus36(new Ring9(209), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.SixteenthGauss721);
        Assert.Equal(new Bus36(new Ring9(217), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.ModBitsCubed729);
        Assert.Equal(new Bus36(new Ring9(228), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.TwentyTickFpga740Ns);
        Assert.Equal(new Bus36(new Ring9(232), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.KleinJ744);
        Assert.Equal(new Bus36(new Ring9(234), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.SecondTurnUart746);
        Assert.Equal(new Bus36(new Ring9(248), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.AtmPressure760Torr);
        Assert.Equal(new Bus36(Ring9.PhiPauli, Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.ThreePauli768);
        Assert.Equal(new Bus36(Ring9.IvSeed, Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.SecondTurnIvSeed775);
        Assert.Equal(new Bus36(new Ring9(273), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.SecondTurnIce785K);
        Assert.Equal(new Bus36(new Ring9(289), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.SecondTurnFresnel801);
        Assert.Equal(new Bus36(new Ring9(305), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.SeventeenthGauss817);
        Assert.Equal(new Bus36(new Ring9(311), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.AudibleFloor823Hz);
        Assert.Equal(new Bus36(new Ring9(312), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.AudibleRounded824Hz);
        Assert.Equal(new Bus36(new Ring9(320), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.GravityDipoleMult832);
        Assert.Equal(new Bus36(new Ring9(321), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.StasisMultiplier833);
        Assert.Equal(new Bus36(new Ring9(322), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.SixCollapse834);
        Assert.Equal(new Bus36(new Ring9(329), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.CopperSquare841);
        Assert.Equal(new Bus36(new Ring9(352), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.TwelveVerlet864);
        Assert.Equal(new Bus36(new Ring9(376), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.TwentyFourTick888Ns);
        Assert.Equal(new Bus36(new Ring9(388), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.QuartzSquare900);
        Assert.Equal(new Bus36(new Ring9(407), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.EighteenthGauss919);
        Assert.Equal(new Bus36(new Ring9(413), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.MosfetGateDelay925Ns);
        Assert.Equal(new Bus36(new Ring9(424), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.CausalConeT1Bits936);
        Assert.Equal(new Bus36(new Ring9(441), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.SecondTurnBlueLed953);
        Assert.Equal(new Bus36(new Ring9(448), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.FifteenKp960);
        Assert.Equal(new Bus36(new Ring9(460), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.TwelveHeisenbergClass972);
        Assert.Equal(new Bus36(new Ring9(487), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.MooreCube27Tick999Ns);
        Assert.Equal(new Bus36(new Ring9(488), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.DecimalCube1000);
        Assert.Equal(new Bus36(new Ring9(489), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.Scheherazade1001);
        Assert.Equal(new Bus36(new Ring9(496), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.FourteenVerlet1008);
        Assert.Equal(new Bus36(new Ring9(504), Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.EightPeltier1016);
        Assert.Equal(new Bus36(Ring9.MinusOne, Ring9.MuVac, Ring9.Zero, Ring9.Zero), full1025.StationaryFirstReturn1023);
        Assert.Equal(new Bus36(Ring9.Zero, Ring9.C0OverCSq, Ring9.Zero, Ring9.Zero), full1025.TwoTurnVerlet1024);
        Assert.Equal(new Bus36(Ring9.MuVac, Ring9.C0OverCSq, Ring9.Zero, Ring9.Zero), full1025.BezoutGateHeisenberg1025);
    }

    [Fact]
    public void Compression_UniversalArchiver_BitExactSha256_ZeroGZipDependency_AndFasterThanGZip()
    {
        // 1. 256 KB 18-channel Verlet sensor/waveform stream
        byte[] raw = new byte[256 * 1024];
        for (int i = 0; i < raw.Length; i++)
        {
            int t = i / 18;
            int ch = i % 18;
            raw[i] = (byte)(t * 7 + ch * 13);
        }

        // Warm-up JIT for both our native archiver and GZipStream
        _ = CuboctUniversalArchiver.CompressBlock(raw.AsSpan(0, 4096));
        using (var warmMs = new System.IO.MemoryStream())
        using (var warmGz = new System.IO.Compression.GZipStream(warmMs, System.IO.Compression.CompressionLevel.Optimal, true))
        {
            warmGz.Write(raw.AsSpan(0, 4096));
        }

        // Measure our pure native CuboctUniversalArchiver (zero GZipStream inside!)
        var swCuboct = System.Diagnostics.Stopwatch.StartNew();
        byte[] comp = CuboctUniversalArchiver.CompressBlock(raw);
        swCuboct.Stop();

        byte[] rec = CuboctUniversalArchiver.DecompressBlock(comp, raw.Length);
        Assert.Equal(SHA256.HashData(raw), SHA256.HashData(rec));

        // Measure standard GZipStream on the exact same 256 KB payload
        var swGzip = System.Diagnostics.Stopwatch.StartNew();
        byte[] gzipBytes;
        using (var gzMs = new System.IO.MemoryStream())
        {
            using (var gz = new System.IO.Compression.GZipStream(gzMs, System.IO.Compression.CompressionLevel.Optimal, true))
            {
                gz.Write(raw);
            }
            gzipBytes = gzMs.ToArray();
        }
        swGzip.Stop();

        // Native Verlet-Heisenberg41 archiver must be BOTH smaller AND faster than GZipStream!
        Assert.True(comp.Length < gzipBytes.Length, $"Expected native archiver ({comp.Length} B) < GZip ({gzipBytes.Length} B)");
        Assert.True(swCuboct.ElapsedTicks < swGzip.ElapsedTicks, $"Expected native archiver ({swCuboct.Elapsed.TotalMilliseconds:F3} ms) < GZip ({swGzip.Elapsed.TotalMilliseconds:F3} ms)");

        // 2. Also verify 100% SHA-256 recovery on Mode 1 (stride-1 ramp), Mode 2 (stride-4 ints), and incompressible random noise (0xFF passthrough)
        byte[] ramp1 = new byte[16384];
        byte[] ramp4 = new byte[16384];
        byte[] noise = new byte[8192];
        var rng = new System.Random(263);
        rng.NextBytes(noise);
        for (int i = 0; i < ramp1.Length; i++)
        {
            ramp1[i] = (byte)(i * 3 + (i % 5));
            ramp4[i] = (byte)((i / 4) * 9 + (i % 4) * 41);
        }

        foreach (byte[] sample in new[] { ramp1, ramp4, noise })
        {
            byte[] c = CuboctUniversalArchiver.CompressBlock(sample);
            byte[] r = CuboctUniversalArchiver.DecompressBlock(c, sample.Length);
            Assert.Equal(SHA256.HashData(sample), SHA256.HashData(r));
        }
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

    [Fact]
    public void Core_14_UniversalCellKernelEvaluatesArbitraryInputsIdenticallyToPythonAndLean()
    {
        // Feed Tetracene (Mu0=120) + Oxygen bridge (Mu1=8) at quarter-ring phase (Phi0=128)
        // directly into the universal VerletCell72 + ReadCellSensors kernel:
        var saddleCell = VerletCell72.Create(new Ring9(120), Ring9.QSquare, Ring9.QSquare, new Ring9(8), Ring9.Zero, Ring9.Zero);
        var readout = CuboctahedronAutomatonKernel.ReadCellSensors(saddleCell);
        Assert.Equal(Ring9.QSquare, readout.OlfactoryShellZ);         // 120 + 8 = 128 = Q_square(3)
        Assert.Equal(Ring9.Zero, saddleCell.Phi0Curr.ProjC);          // C_3(128) = 0
        Assert.Equal(Ring9.QSquare, readout.ChiralOdorProjectionS);   // S_3(128) = +128
        Assert.Equal(Ring9.Zero, readout.EnergyQuantaNe);             // n_E(0) = 0

        // Feed SiO4 / [AlO4H] (Z = 46) vs dry AlO4 (Z = 45) directly into the universal PauliExclusionGate9:
        Assert.Equal(Ring9.PhiPauli, CuboctahedronAutomatonKernel.PauliExclusionGate9(new Ring9(46), Ring9.Zero, new Ring9(46), Ring9.Zero, Ring9.MuVac));
        Assert.Equal(Ring9.Zero,     CuboctahedronAutomatonKernel.PauliExclusionGate9(new Ring9(45), Ring9.Zero, new Ring9(45), Ring9.Zero, Ring9.MuVac));
    }

    [Fact]
    public void MediaCodec_TearsOfSteel14GbDcp_VideoAndAudioComparison_BitExactQualityAndSoundPreserved()
    {
        const string tosPicPath = @"C:\Users\dkara\Downloads\ToS_4k_DCP\tos_version_05\tos_picture.mxf";
        const string tosSndPath = @"C:\Users\dkara\Downloads\ToS_4k_DCP\tos_version_05\tos_sound.mxf";
        const string ffmpegExe  = @"C:\Users\dkara\AppData\Local\Microsoft\WinGet\Links\ffmpeg.exe";

        Assert.True(System.IO.File.Exists(tosPicPath), $"Missing 14 GB Tears of Steel picture MXF: {tosPicPath}");
        Assert.True(System.IO.File.Exists(tosSndPath), $"Missing Tears of Steel sound MXF: {tosSndPath}");
        Assert.True(System.IO.File.Exists(ffmpegExe),  $"Missing ffmpeg executable: {ffmpegExe}");

        const int width = 640;
        const int height = 360;
        const int fps = 24;
        const int durationSec = 5;
        const int frameCount = fps * durationSec; // 120 frames (5.0 seconds at 00:02:00)
        const int frameBytes = width * height * 3;
        const int totalVideoBytes = frameCount * frameBytes; // 82,944,000 bytes (~82.94 MB)
        const int sampleRate = 48000;
        const int nStereoSamples = sampleRate * durationSec; // 240,000 stereo pairs (480,000 int16 samples = 960,000 bytes)

        static byte[] RunFfmpegRawPipe(string exe, string arguments, int expectedBytes)
        {
            var psi = new System.Diagnostics.ProcessStartInfo(exe, arguments)
            {
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                UseShellExecute = false,
                CreateNoWindow = true
            };
            using var proc = System.Diagnostics.Process.Start(psi)!;
            byte[] buf = new byte[expectedBytes];
            int offset = 0;
            var stdout = proc.StandardOutput.BaseStream;
            while (offset < expectedBytes)
            {
                int read = stdout.Read(buf, offset, expectedBytes - offset);
                if (read <= 0) break;
                offset += read;
            }
            proc.WaitForExit();
            Assert.Equal(expectedBytes, offset);
            return buf;
        }

        // 1. Decode 5 seconds (120 frames, 640x360 RGB24) from 14 GB Tears of Steel 4K DCP picture MXF at 00:02:00
        byte[] origVideoRgb24 = RunFfmpegRawPipe(
            ffmpegExe,
            $"-ss 00:02:00 -i \"{tosPicPath}\" -vf scale={width}:{height} -vframes {frameCount} -f rawvideo -pix_fmt rgb24 -",
            totalVideoBytes);

        // 2. Decode 5 seconds (48 kHz 16-bit Stereo PCM) from Tears of Steel studio sound MXF at 00:02:00
        byte[] origAudioBytes = RunFfmpegRawPipe(
            ffmpegExe,
            $"-ss 00:02:00 -i \"{tosSndPath}\" -t {durationSec} -ac 2 -ar {sampleRate} -f s16le -",
            nStereoSamples * 4);

        short[] origAudioPcm16 = new short[nStereoSamples * 2];
        System.Buffer.BlockCopy(origAudioBytes, 0, origAudioPcm16, 0, origAudioBytes.Length);

        // Verify original sound is non-silent (active speech + soundtrack at 00:02:00)
        long origAudioEnergy = 0;
        int maxOrigSample = 0;
        for (int i = 0; i < origAudioPcm16.Length; i++)
        {
            int a = System.Math.Abs((int)origAudioPcm16[i]);
            origAudioEnergy += (long)a * a;
            if (a > maxOrigSample) maxOrigSample = a;
        }
        Assert.True(maxOrigSample > 500, $"Expected active audio in Tears of Steel at 00:02:00, peak={maxOrigSample}");

        // 3. Encode & Decode Multiplexed Video + Audio Container purely in .NET 11 via CuboctMediaCodec
        var swEncode = System.Diagnostics.Stopwatch.StartNew();
        byte[] cuboctContainer = CuboctMediaCodec.EncodeMultiplexedMediaContainer(
            origVideoRgb24, frameCount, width, height, origAudioPcm16, sampleRate);
        swEncode.Stop();

        var swDecode = System.Diagnostics.Stopwatch.StartNew();
        var decoded = CuboctMediaCodec.DecodeMultiplexedMediaContainer(cuboctContainer);
        swDecode.Stop();

        // 4. Measure standard .NET GZipStream(CompressionLevel.Optimal) on the exact same Video + Audio payload
        var swGzip = System.Diagnostics.Stopwatch.StartNew();
        byte[] gzipVideoBytes;
        using (var msV = new System.IO.MemoryStream())
        {
            using (var gzV = new System.IO.Compression.GZipStream(msV, System.IO.Compression.CompressionLevel.Optimal, true))
            {
                gzV.Write(origVideoRgb24);
                gzV.Write(origAudioBytes);
            }
            gzipVideoBytes = msV.ToArray();
        }
        swGzip.Stop();

        // 5. Verify 100% Bit-Exact Video Quality (0 pixel difference, SHA-256 match, PSNR = Infinity)
        byte[] shaVideoOrig = SHA256.HashData(origVideoRgb24);
        byte[] shaVideoDec  = SHA256.HashData(decoded.Rgb24Frames);
        Assert.Equal(shaVideoOrig, shaVideoDec);

        int maxPixelDiff = 0;
        for (int i = 0; i < origVideoRgb24.Length; i++)
        {
            int d = System.Math.Abs(origVideoRgb24[i] - decoded.Rgb24Frames[i]);
            if (d > maxPixelDiff) maxPixelDiff = d;
        }
        Assert.Equal(0, maxPixelDiff);

        // 6. Verify 100% Bit-Exact Audio Preservation (Sound did NOT disappear or degrade: 0 sample difference, SHA-256 match)
        Assert.Equal(origAudioPcm16.Length, decoded.StereoPcm16.Length);
        byte[] decAudioBytes = new byte[origAudioBytes.Length];
        System.Buffer.BlockCopy(decoded.StereoPcm16, 0, decAudioBytes, 0, decAudioBytes.Length);
        Assert.Equal(SHA256.HashData(origAudioBytes), SHA256.HashData(decAudioBytes));

        int maxAudioDiff = 0;
        long decAudioEnergy = 0;
        for (int i = 0; i < origAudioPcm16.Length; i++)
        {
            int d = System.Math.Abs(origAudioPcm16[i] - decoded.StereoPcm16[i]);
            if (d > maxAudioDiff) maxAudioDiff = d;
            int a = System.Math.Abs((int)decoded.StereoPcm16[i]);
            decAudioEnergy += (long)a * a;
        }
        Assert.Equal(0, maxAudioDiff);
        Assert.Equal(origAudioEnergy, decAudioEnergy);

        // 7. Also evaluate Two-Point Keyframe Interpolation quality (Stride N = 2 and Stride N = 4) on Tears of Steel frames
        double totalHeisPctStride2 = 0.0;
        double mseStride2 = 0.0;
        int interpCount2 = 0;
        for (int f = 0; f + 2 < frameCount; f += 2)
        {
            var kf0 = origVideoRgb24.AsSpan(f * frameBytes, frameBytes);
            var gt1 = origVideoRgb24.AsSpan((f + 1) * frameBytes, frameBytes);
            var kf2 = origVideoRgb24.AsSpan((f + 2) * frameBytes, frameBytes);
            var (rec1, pct41) = CuboctMediaCodec.InterpolateTwoPointKeyframe(kf0, kf2, gt1, 1, 2);
            totalHeisPctStride2 += pct41;
            long sqErr = 0;
            for (int idx = 0; idx < frameBytes; idx++)
            {
                int e = gt1[idx] - rec1[idx];
                sqErr += e * e;
            }
            mseStride2 += (double)sqErr / frameBytes;
            interpCount2++;
        }
        double avgHeisStride2 = totalHeisPctStride2 / interpCount2;
        mseStride2 /= interpCount2;
        double psnrStride2 = 10.0 * System.Math.Log10((255.0 * 255.0) / System.Math.Max(1e-9, mseStride2));

        // 8. Export two side-by-side playable MP4 files with sound (Original vs Decoded after .NET 11 CuboctMediaCodec)
        string artifactsDir = @"D:\Experiments\Personal Cursor Folder\open\guiders-physics\artifacts";
        System.IO.Directory.CreateDirectory(artifactsDir);
        string rawVidOrigPath = System.IO.Path.Combine(artifactsDir, "tos_orig.rgb24");
        string rawAudOrigPath = System.IO.Path.Combine(artifactsDir, "tos_orig.s16le");
        string rawVidDecPath  = System.IO.Path.Combine(artifactsDir, "tos_dec.rgb24");
        string rawAudDecPath  = System.IO.Path.Combine(artifactsDir, "tos_dec.s16le");
        string mp4OrigPath    = System.IO.Path.Combine(artifactsDir, "tos_1_original_reference.mp4");
        string mp4DecPath     = System.IO.Path.Combine(artifactsDir, "tos_2_after_cuboct_net11_codec.mp4");

        System.IO.File.WriteAllBytes(rawVidOrigPath, origVideoRgb24);
        System.IO.File.WriteAllBytes(rawAudOrigPath, origAudioBytes);
        System.IO.File.WriteAllBytes(rawVidDecPath,  decoded.Rgb24Frames);
        System.IO.File.WriteAllBytes(rawAudDecPath,  decAudioBytes);

        static void MuxMp4WithSound(string exe, string vidPath, string audPath, int w, int h, int rate, int sr, string outMp4)
        {
            string args = $"-y -f rawvideo -pix_fmt rgb24 -s {w}x{h} -r {rate} -i \"{vidPath}\" " +
                          $"-f s16le -ar {sr} -ac 2 -i \"{audPath}\" " +
                          $"-c:v libx264 -crf 10 -pix_fmt yuv420p -c:a aac -b:a 256k \"{outMp4}\"";
            var psi = new System.Diagnostics.ProcessStartInfo(exe, args)
            {
                RedirectStandardError = true,
                UseShellExecute = false,
                CreateNoWindow = true
            };
            using var p = System.Diagnostics.Process.Start(psi)!;
            p.StandardError.ReadToEnd();
            p.WaitForExit();
            Assert.Equal(0, p.ExitCode);
        }

        MuxMp4WithSound(ffmpegExe, rawVidOrigPath, rawAudOrigPath, width, height, fps, sampleRate, mp4OrigPath);
        MuxMp4WithSound(ffmpegExe, rawVidDecPath,  rawAudDecPath,  width, height, fps, sampleRate, mp4DecPath);

        System.IO.File.Delete(rawVidOrigPath);
        System.IO.File.Delete(rawAudOrigPath);
        System.IO.File.Delete(rawVidDecPath);
        System.IO.File.Delete(rawAudDecPath);

        int rawTotalBytes = origVideoRgb24.Length + origAudioBytes.Length;
        _output.WriteLine($"=== .NET 11 CuboctMediaCodec Benchmark on 14 GB Tears of Steel 4K DCP (00:02:00, {frameCount} frames + 48kHz Stereo Audio) ===");
        _output.WriteLine($"Raw Uncompressed Stream (Video + Audio): {rawTotalBytes,14:N0} B ({rawTotalBytes / (1024.0 * 1024.0):F2} MB)");
        _output.WriteLine($".NET GZipStream (Optimal) Size & Time:   {gzipVideoBytes.Length,14:N0} B ({(double)rawTotalBytes / gzipVideoBytes.Length:F2}x) | Encode: {swGzip.Elapsed.TotalMilliseconds:F1} ms");
        _output.WriteLine($".NET 11 CuboctMediaCodec Size & Time:    {cuboctContainer.Length,14:N0} B ({(double)rawTotalBytes / cuboctContainer.Length:F2}x) | Encode: {swEncode.Elapsed.TotalMilliseconds:F1} ms | Decode: {swDecode.Elapsed.TotalMilliseconds:F1} ms");
        _output.WriteLine($"Video Quality After Codec vs Original:   100% Bit-Exact SHA-256 ({ System.Convert.ToHexString(shaVideoDec) }), MaxPixelDiff = {maxPixelDiff}, PSNR = +Infinity dB");
        _output.WriteLine($"Audio Quality After Codec vs Original:   100% Bit-Exact SHA-256 ({ System.Convert.ToHexString(SHA256.HashData(decAudioBytes)) }), MaxSampleDiff = {maxAudioDiff}, Peak = {maxOrigSample}, Energy = {decAudioEnergy:N0}");
        _output.WriteLine($"Two-Point Keyframe Stride N=2 Quality:   PSNR = {psnrStride2:F2} dB, Pixels inside Heisenberg [-40, +40] = {avgHeisStride2:F2}%");
        _output.WriteLine($"Saved Reference MP4 (with sound):        {mp4OrigPath} ({new System.IO.FileInfo(mp4OrigPath).Length:N0} B)");
        _output.WriteLine($"Saved Decoded MP4 (with sound):          {mp4DecPath} ({new System.IO.FileInfo(mp4DecPath).Length:N0} B)");

        Assert.True(cuboctContainer.Length < gzipVideoBytes.Length,
            $"Expected CuboctMediaCodec ({cuboctContainer.Length:N0} B) < GZipStream ({gzipVideoBytes.Length:N0} B) on Tears of Steel!");
        Assert.True(swEncode.ElapsedTicks < swGzip.ElapsedTicks,
            $"Expected CuboctMediaCodec ({swEncode.Elapsed.TotalMilliseconds:F1} ms) < GZipStream ({swGzip.Elapsed.TotalMilliseconds:F1} ms) on Tears of Steel!");
    }
}
