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
        Ring9 center = Ring9.NLinks; // 12
        Ring9[,,] mu = new Ring9[N, N, N];
        Ring9[,,] p0 = new Ring9[N, N, N];
        Ring9[,,] p1 = new Ring9[N, N, N];
        Ring9[,,] q0 = new Ring9[N, N, N];
        Ring9[,,] q1 = new Ring9[N, N, N];

        for (ushort iz = 0; iz < N; iz++)
            for (ushort iy = 0; iy < N; iy++)
                for (ushort ix = 0; ix < N; ix++)
                {
                    Ring9 dz = new Ring9(iz) - center;
                    Ring9 dy = new Ring9(iy) - center;
                    Ring9 dx = new Ring9(ix) - center;
                    mu[iz, iy, ix] = CuboctahedronAutomatonKernel.Iv263MuAt(dz, dy, dx);
                    p0[iz, iy, ix] = CuboctahedronAutomatonKernel.Iv263PhiAt(dz, dy, dx, Ring9.Zero);
                    p1[iz, iy, ix] = CuboctahedronAutomatonKernel.Iv263PhiAt(dz, dy, dx, Ring9.MuVac);
                    q0[iz, iy, ix] = p0[iz, iy, ix];
                    q1[iz, iy, ix] = p1[iz, iy, ix];
                }

        static (Ring9[,,] Nxt0, Ring9[,,] Nxt1, int ActiveGates) StepGrid3D(
            Ring9[,,] m, Ring9[,,] c0, Ring9[,,] c1, Ring9[,,] pst0, Ring9[,,] pst1)
        {
            Ring9[,,] n0 = new Ring9[N, N, N];
            Ring9[,,] n1 = new Ring9[N, N, N];
            int active = 0;

            Ring9 M(int z, int y, int x) => m[(z + N) % N, (y + N) % N, (x + N) % N];
            Ring9 P0(int z, int y, int x) => c0[(z + N) % N, (y + N) % N, (x + N) % N];
            Ring9 P1(int z, int y, int x) => c1[(z + N) % N, (y + N) % N, (x + N) % N];

            for (int z = 0; z < N; z++)
                for (int y = 0; y < N; y++)
                    for (int x = 0; x < N; x++)
                    {
                        Ring9 gate = CuboctahedronAutomatonKernel.ComputeTwoLaneCellGate(
                            M(z, y, x), P0(z, y, x), M(z, y, x), P1(z, y, x),
                            M(z, y - 1, x - 1), P0(z, y - 1, x - 1), M(z, y - 1, x + 1), P0(z, y - 1, x + 1),
                            M(z, y + 1, x - 1), P0(z, y + 1, x - 1), M(z, y + 1, x + 1), P0(z, y + 1, x + 1),
                            M(z - 1, y, x - 1), P0(z - 1, y, x - 1), M(z - 1, y, x + 1), P0(z - 1, y, x + 1),
                            M(z + 1, y, x - 1), P0(z + 1, y, x - 1), M(z + 1, y, x + 1), P0(z + 1, y, x + 1),
                            M(z - 1, y - 1, x), P0(z - 1, y - 1, x), M(z - 1, y + 1, x), P0(z - 1, y + 1, x),
                            M(z + 1, y - 1, x), P0(z + 1, y - 1, x), M(z + 1, y + 1, x), P0(z + 1, y + 1, x),
                            M(z, y - 1, x - 1), P1(z, y - 1, x - 1), M(z, y - 1, x + 1), P1(z, y - 1, x + 1),
                            M(z, y + 1, x - 1), P1(z, y + 1, x - 1), M(z, y + 1, x + 1), P1(z, y + 1, x + 1),
                            M(z - 1, y, x - 1), P1(z - 1, y, x - 1), M(z - 1, y, x + 1), P1(z - 1, y, x + 1),
                            M(z + 1, y, x - 1), P1(z + 1, y, x - 1), M(z + 1, y, x + 1), P1(z + 1, y, x + 1),
                            M(z - 1, y - 1, x), P1(z - 1, y - 1, x), M(z - 1, y + 1, x), P1(z - 1, y + 1, x),
                            M(z + 1, y - 1, x), P1(z + 1, y - 1, x), M(z + 1, y + 1, x), P1(z + 1, y + 1, x),
                            Ring9.MuVac);

                        if (gate != Ring9.Zero) active++;
                        var cell = VerletCell72.Create(M(z, y, x), c0[z, y, x], pst0[z, y, x], M(z, y, x), c1[z, y, x], pst1[z, y, x]);
                        var stepped = CuboctahedronAutomatonKernel.StepCellForward(cell, gate, gate);
                        n0[z, y, x] = stepped.Phi0Curr;
                        n1[z, y, x] = stepped.Phi1Curr;
                    }
            return (n0, n1, active);
        }

        var (st1_0, st1_1, act0) = StepGrid3D(mu, p0, p1, q0, q1);
        var (st2_0, _,     act1) = StepGrid3D(mu, st1_0, st1_1, p0, p1);

        Assert.Equal(543, act0);
        Assert.Equal(578, act1);

        Ring9 P2(int z, int y, int x) => st2_0[(z + N) % N, (y + N) % N, (x + N) % N];
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
        const int N = 7;
        const int C = 3;

        static (Ring9[,,] Mu, Ring9[,,] P0, Ring9[,,] P1) MakeReactionGrid(Ring9 phiEnv)
        {
            Ring9[,,] m = new Ring9[N, N, N];
            Ring9[,,] p0 = new Ring9[N, N, N];
            Ring9[,,] p1 = new Ring9[N, N, N];
            for (int iz = 0; iz < N; iz++)
            {
                int dz = iz - C;
                for (int iy = 0; iy < N; iy++)
                {
                    int dy = iy - C;
                    for (int ix = 0; ix < N; ix++)
                    {
                        int dx = ix - C;
                        m[iz, iy, ix] = Ring9.NPhi;
                        Ring9 p = (dz == 0 && dy == 0 && dx == 0)
                            ? Ring9.Zero
                            : ((4 * dz + 2 * dy + dx) >= 0 ? phiEnv : -phiEnv);
                        p0[iz, iy, ix] = p;
                        p1[iz, iy, ix] = p;
                    }
                }
            }
            return (m, p0, p1);
        }

        static (Ring9[,,] Nxt0, Ring9[,,] Nxt1) StepReactionGrid(
            Ring9[,,] m, Ring9[,,] c0, Ring9[,,] c1, Ring9[,,] pst0, Ring9[,,] pst1)
        {
            Ring9[,,] n0 = new Ring9[N, N, N];
            Ring9[,,] n1 = new Ring9[N, N, N];
            Ring9 M(int z, int y, int x) => m[(z + N) % N, (y + N) % N, (x + N) % N];
            Ring9 P0(int z, int y, int x) => c0[(z + N) % N, (y + N) % N, (x + N) % N];
            Ring9 P1(int z, int y, int x) => c1[(z + N) % N, (y + N) % N, (x + N) % N];

            for (int z = 0; z < N; z++)
                for (int y = 0; y < N; y++)
                    for (int x = 0; x < N; x++)
                    {
                        Ring9 gate = CuboctahedronAutomatonKernel.ComputeTwoLaneCellGate(
                            M(z, y, x), P0(z, y, x), M(z, y, x), P1(z, y, x),
                            M(z, y - 1, x - 1), P0(z, y - 1, x - 1), M(z, y - 1, x + 1), P0(z, y - 1, x + 1),
                            M(z, y + 1, x - 1), P0(z, y + 1, x - 1), M(z, y + 1, x + 1), P0(z, y + 1, x + 1),
                            M(z - 1, y, x - 1), P0(z - 1, y, x - 1), M(z - 1, y, x + 1), P0(z - 1, y, x + 1),
                            M(z + 1, y, x - 1), P0(z + 1, y, x - 1), M(z + 1, y, x + 1), P0(z + 1, y, x + 1),
                            M(z - 1, y - 1, x), P0(z - 1, y - 1, x), M(z - 1, y + 1, x), P0(z - 1, y + 1, x),
                            M(z + 1, y - 1, x), P0(z + 1, y - 1, x), M(z + 1, y + 1, x), P0(z + 1, y + 1, x),
                            M(z, y - 1, x - 1), P1(z, y - 1, x - 1), M(z, y - 1, x + 1), P1(z, y - 1, x + 1),
                            M(z, y + 1, x - 1), P1(z, y + 1, x - 1), M(z, y + 1, x + 1), P1(z, y + 1, x + 1),
                            M(z - 1, y, x - 1), P1(z - 1, y, x - 1), M(z - 1, y, x + 1), P1(z - 1, y, x + 1),
                            M(z + 1, y, x - 1), P1(z + 1, y, x - 1), M(z + 1, y, x + 1), P1(z + 1, y, x + 1),
                            M(z - 1, y - 1, x), P1(z - 1, y - 1, x), M(z - 1, y + 1, x), P1(z - 1, y + 1, x),
                            M(z + 1, y - 1, x), P1(z + 1, y - 1, x), M(z + 1, y + 1, x), P1(z + 1, y + 1, x),
                            Ring9.MuVac);
                        var cell = VerletCell72.Create(M(z, y, x), c0[z, y, x], pst0[z, y, x], M(z, y, x), c1[z, y, x], pst1[z, y, x]);
                        var stepped = CuboctahedronAutomatonKernel.StepCellForward(cell, gate, gate);
                        n0[z, y, x] = stepped.Phi0Curr;
                        n1[z, y, x] = stepped.Phi1Curr;
                    }
            return (n0, n1);
        }

        // 1) Inert solvent (phi_env = 41): reaction does NOT proceed (center phase stays 0 -> 0)
        var (mInert, p0Inert, p1Inert) = MakeReactionGrid(Ring9.DPhiDisc);
        var (n0Inert, _) = StepReactionGrid(mInert, p0Inert, p1Inert, p0Inert, p1Inert);
        Assert.Equal(Ring9.Zero, n0Inert[C, C, C]);

        // 2) Catalytic Belousov-Zhabotinsky oscillation (phi_env = 139 = CollapsePhi):
        // Exact 13-tick center trajectory from mt_ca/carrier_core.py:
        // [0, 352, 361, 212, 214, 61, 262, 301, 95, 52, 121, 90, 59, 427]
        ushort[] expectedBzPhases = [0, 352, 361, 212, 214, 61, 262, 301, 95, 52, 121, 90, 59, 427];
        var (mCat, cur0, cur1) = MakeReactionGrid(Ring9.CollapsePhi);
        var pst0 = cur0;
        var pst1 = cur1;

        for (int t = 1; t <= 13; t++)
        {
            var (nxt0, nxt1) = StepReactionGrid(mCat, cur0, cur1, pst0, pst1);
            pst0 = cur0;
            pst1 = cur1;
            cur0 = nxt0;
            cur1 = nxt1;

            Assert.Equal(new Ring9(expectedBzPhases[t]), cur0[C, C, C]);
            var centerCell = VerletCell72.Create(Ring9.NPhi, cur0[C, C, C], pst0[C, C, C], Ring9.NPhi, cur1[C, C, C], pst1[C, C, C]);
            var ro = CuboctahedronAutomatonKernel.ReadCellSensors(centerCell);
            Assert.Equal(Ring9.QSquare, ro.ColorR + ro.ColorG + ro.ColorB);

            if (t == 1)
            {
                // At t=1 (phi=352, oxidized blue-green ferroin state): RGB = (0, 96, 32)
                Assert.Equal(Ring9.Zero, ro.ColorR);
                Assert.Equal(new Ring9(96), ro.ColorG);
                Assert.Equal(new Ring9(32), ro.ColorB);
            }
            else if (t == 5)
            {
                // At t=5 (phi=61, reduced red-orange ferroin state): RGB = (67, 61, 0)
                Assert.Equal(new Ring9(67), ro.ColorR);
                Assert.Equal(new Ring9(61), ro.ColorG);
                Assert.Equal(Ring9.Zero, ro.ColorB);
            }
        }
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
