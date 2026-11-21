using System.IO;
using System.Security.Cryptography;
using System.Text.Json;
using Guiders.Physics.Core;
using Guiders.Physics.Compression;
using Guiders.Physics.MediaCodec;
using Xunit;

namespace Guiders.Physics.Tests;

public class GuidersPhysicsSuiteTests
{
    private static int ToSigned512(Ring9 r) =>
        r.Bits9 > 256 ? r.Bits9 - 512 : r.Bits9;

    private static long ToSignedBus36(Bus36 b)
    {
        ulong u = (ulong)b.L0.Bits9 | ((ulong)b.L1.Bits9 << 9) | ((ulong)b.L2.Bits9 << 18) | ((ulong)b.L3.Bits9 << 27);
        return (u & (1UL << 35)) != 0 ? (long)u - (1L << 36) : (long)u;
    }

    [Fact]
    public unsafe void Core_1_PassportAnd72BitLayout_MatchesLean4AndPython()
    {
        Assert.Equal(9, sizeof(VerletCell72));
        Assert.Equal(72, sizeof(VerletCell72) * 8);

        Assert.Equal(3,   Ring9.D.Bits9);
        Assert.Equal(1,   Ring9.MuVac.Bits9);
        Assert.Equal(2,   Ring9.C0OverCSq.Bits9);
        Assert.Equal(12,  Ring9.NLinks.Bits9);
        Assert.Equal(8,   Ring9.NTri.Bits9);
        Assert.Equal(6,   Ring9.NSq.Bits9);
        Assert.Equal(4,   Ring9.DInSquare.Bits9);
        Assert.Equal(5,   Ring9.DOutEdge.Bits9);
        Assert.Equal(7,   Ring9.DOutVertex.Bits9);
        Assert.Equal(9,   Ring9.ModBits.Bits9);
        Assert.Equal(13,  Ring9.NPhi.Bits9);
        Assert.Equal(41,  Ring9.DPhiDisc.Bits9);
        Assert.Equal(21,  Ring9.SeamTicks.Bits9);
        Assert.Equal(64,  Ring9.KP.Bits9);
        Assert.Equal(97,  Ring9.MHat.Bits9);
        Assert.Equal(128, Ring9.QSquare.Bits9);
        Assert.Equal(139, Ring9.CollapsePhi.Bits9);
        Assert.Equal(256, Ring9.PhiPauli.Bits9);
        Assert.Equal(263, Ring9.IvSeed.Bits9);

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
    public void Core_2_Rotator512AndPythagorean_MatchesLean4JsonDumpAll512()
    {
        for (ushort phi = 0; phi < 512; phi++)
        {
            Ring9 r = new(phi);
            int c = ToSigned512(r.ProjC);
            int s = ToSigned512(r.ProjS);

            int wrap = 256 - ((256 - phi) & 511);
            int expectedC = 128 - System.Math.Abs(wrap);
            int wrapS = 256 - ((256 - ((phi - 128) & 511)) & 511);
            int expectedS = 128 - System.Math.Abs(wrapS);

            Assert.Equal(expectedC, c);
            Assert.Equal(expectedS, s);
            Assert.Equal(128, System.Math.Abs(c) + System.Math.Abs(s));

            var (pxBus, pyBus, prBus) = CuboctahedronAutomatonKernel.PythagoreanRotator36(r);
            long px = ToSignedBus36(pxBus);
            long py = ToSignedBus36(pyBus);
            long pr = ToSignedBus36(prBus);

            Assert.Equal(pr * pr, px * px + py * py);
            if (phi == 41)
            {
                Assert.Equal(14703L, px);
                Assert.Equal(10496L, py);
                Assert.Equal(18065L, pr);
            }
        }

        string dumpPath = @"D:\Experiments\Personal Cursor Folder\open\mt-lattice-fluid\lean\mt_carrier_lean_dump.json";
        if (File.Exists(dumpPath))
        {
            using var doc = JsonDocument.Parse(File.ReadAllText(dumpPath));
            var root = doc.RootElement;
            var cLut = root.GetProperty("c_lut");
            var sLut = root.GetProperty("s_lut");
            var pxLut = root.GetProperty("px_lut");
            var pyLut = root.GetProperty("py_lut");
            var prLut = root.GetProperty("pr_lut");

            for (ushort phi = 0; phi < 512; phi++)
            {
                Ring9 r = new(phi);
                Assert.Equal(cLut[phi].GetInt32(), ToSigned512(r.ProjC));
                Assert.Equal(sLut[phi].GetInt32(), ToSigned512(r.ProjS));
                var (pxB, pyB, prB) = CuboctahedronAutomatonKernel.PythagoreanRotator36(r);
                Assert.Equal(pxLut[phi].GetInt64(), ToSignedBus36(pxB));
                Assert.Equal(pyLut[phi].GetInt64(), ToSignedBus36(pyB));
                Assert.Equal(prLut[phi].GetInt64(), ToSignedBus36(prB));
            }
        }
    }

    [Fact]
    public void Core_3_Channel2BranchCarryAndPlaquetteWinding_ExactRuleZ3()
    {
        for (ushort pY = 0; pY < 512; pY += 7)
        {
            for (ushort pX = 0; pX < 512; pX += 11)
            {
                var (_, kCarry) = CuboctahedronAutomatonKernel.EdgeDifferenceWithBranchCarry(new Ring9(pY), new Ring9(pX));
                int dphi = pY - pX;
                int expectedK = (dphi + 255) >= 0 ? (dphi + 255) / 512 : -((-(dphi + 255) + 511) / 512);
                Assert.Equal(expectedK, ToSigned512(kCarry));
            }
        }

        // Plaquette winding around a vortex loop (0 -> 130 -> 260 -> 390 -> 0):
        Ring9 w = CuboctahedronAutomatonKernel.PlaquetteWinding9(new(0), new(130), new(260), new(390));
        Assert.Equal(1, ToSigned512(w));
    }

    [Fact]
    public void Core_4_PauliExclusionAndVerletReversibility_MatchesPythonCheckPauliRepel()
    {
        // Matches check_pauli_repel in verify_checks/axioms.py (line 334):
        // Two dense parallel spinor lanes: mu0 = mu1 = 46 (46^2 = 2116 > 2048), phi0 = phi1 = 0 -> +256 Pauli jump!
        Ring9 pOn = CuboctahedronAutomatonKernel.PauliExclusionGate9(new(46), new(0), new(46), new(0), Ring9.MuVac);
        Ring9 pOff = CuboctahedronAutomatonKernel.PauliExclusionGate9(new(46), new(0), new(46), new(0), Ring9.Zero);
        Ring9 pAnti = CuboctahedronAutomatonKernel.PauliExclusionGate9(new(46), new(0), new(46), Ring9.PhiPauli, Ring9.MuVac);
        Ring9 pSparse = CuboctahedronAutomatonKernel.PauliExclusionGate9(new(13), new(0), new(13), new(0), Ring9.MuVac);

        Assert.Equal(Ring9.PhiPauli, pOn);
        Assert.Equal(Ring9.Zero, pOff);
        Assert.Equal(Ring9.Zero, pAnti);
        Assert.Equal(Ring9.Zero, pSparse);

        // Exhaustive 512/512 history-free Verlet reversibility (polarReverseCell(polarStepCell(c)) == c):
        for (ushort p = 0; p < 512; p++)
        {
            Ring9 curr = new(p);
            Ring9 past = new((ushort)((p * 13 + 97) & 0x1FF));
            Ring9 gate = new((ushort)((p * 41) & 0x1FF));
            var c0 = VerletCell72.Create(new(13), curr, past, new(13), curr + Ring9.DPhiDisc, past + Ring9.DPhiDisc);
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
        // 1. Single-node holonomy checks from TrajectoryIV263.lean:
        Assert.Equal(Ring9.Zero, CuboctahedronAutomatonKernel.Compute12NeighborHolonomyGateAt(new(2), new(2), new(2), Ring9.Zero));
        Ring9 gWallPos = CuboctahedronAutomatonKernel.Compute12NeighborHolonomyGateAt(new(0), new(2), new(2), Ring9.MuVac);
        Ring9 gWallNeg = CuboctahedronAutomatonKernel.Compute12NeighborHolonomyGateAt(new(511), new(2), new(2), Ring9.MuVac);
        Assert.Equal(new Ring9(76), gWallPos);
        Assert.Equal(new Ring9(512 - 77), gWallNeg);
        Assert.Equal(Ring9.MinusOne, gWallPos + gWallNeg);

        // 2. Full 3D 24x24x24 IV(3)=263 t=0..2 evolution using purely Guiders.Physics.Core (matching verify_checks/carrier.py lines 892-919!):
        const int N = 24;
        const int C = N / 2; // 12
        Ring9[,,] mu = new Ring9[N, N, N];
        Ring9[,,] p0 = new Ring9[N, N, N];
        Ring9[,,] p1 = new Ring9[N, N, N];
        Ring9[,,] q0 = new Ring9[N, N, N];
        Ring9[,,] q1 = new Ring9[N, N, N];

        for (int iz = 0; iz < N; iz++)
            for (int iy = 0; iy < N; iy++)
                for (int ix = 0; ix < N; ix++)
                {
                    Ring9 dz = new((ushort)((iz - C) & 0x1FF));
                    Ring9 dy = new((ushort)((iy - C) & 0x1FF));
                    Ring9 dx = new((ushort)((ix - C) & 0x1FF));
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

                        if (gate.Bits9 != 0) active++;
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
        int NXY(int z, int y, int x) => ToSigned512(CuboctahedronAutomatonKernel.PlaquetteWinding9(P2(z, y, x), P2(z, y, x + 1), P2(z, y + 1, x + 1), P2(z, y + 1, x)));
        int NXZ(int z, int y, int x) => ToSigned512(CuboctahedronAutomatonKernel.PlaquetteWinding9(P2(z, y, x), P2(z, y, x + 1), P2(z + 1, y, x + 1), P2(z + 1, y, x)));
        int NYZ(int z, int y, int x) => ToSigned512(CuboctahedronAutomatonKernel.PlaquetteWinding9(P2(z, y, x), P2(z, y + 1, x), P2(z + 1, y + 1, x), P2(z + 1, y, x)));

        int sumXY = 0, sumYZ = 0, sumXZ = 0, maxDiv = 0;
        for (int z = 0; z < N; z++)
            for (int y = 0; y < N; y++)
                for (int x = 0; x < N; x++)
                {
                    sumXY += System.Math.Abs(NXY(z, y, x));
                    sumYZ += System.Math.Abs(NYZ(z, y, x));
                    sumXZ += System.Math.Abs(NXZ(z, y, x));
                    int divCube = (NXY(z + 1, y, x) - NXY(z, y, x))
                                - (NXZ(z, y + 1, x) - NXZ(z, y, x))
                                + (NYZ(z, y, x + 1) - NYZ(z, y, x));
                    maxDiv = System.Math.Max(maxDiv, System.Math.Abs(divCube));
                }

        Assert.Equal(26, sumXY);
        Assert.Equal(26, sumYZ);
        Assert.Equal(34, sumXZ);
        Assert.Equal(86, sumXY + sumYZ + sumXZ);
        Assert.Equal(0, maxDiv);
    }

    [Fact]
    public void Core_6_MultimodalSensorReadout_ExactRgbTriadAndIsqrt4096()
    {
        for (ushort q = 0; q <= 4096; q++)
        {
            Assert.Equal((ushort)System.Math.Sqrt(q), CuboctahedronAutomatonKernel.Isqrt4096(q).Bits9);
        }

        for (ushort phi = 0; phi < 512; phi++)
        {
            var c = VerletCell72.Create(new(13), new(phi), new(phi), new(13), new(phi), new(phi));
            c.Gate0 = new Ring9(139);
            var s = CuboctahedronAutomatonKernel.ReadCellSensors(c);
            Assert.Equal(128, s.ColorR.Bits9 + s.ColorG.Bits9 + s.ColorB.Bits9);
            Assert.Equal(3, s.EnergyQuantaNe.Bits9);
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
            raw[i] = (byte)((t * 7 + ch * 13) & 0xFF);
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
    }
}
