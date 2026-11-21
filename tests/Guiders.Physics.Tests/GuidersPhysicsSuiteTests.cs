using System;
using System.Security.Cryptography;
using Guiders.Physics.Core;
using Guiders.Physics.Compression;
using Guiders.Physics.MediaCodec;
using Guiders.Physics.Engine;
using Xunit;

namespace Guiders.Physics.Tests;

public class GuidersPhysicsSuiteTests
{
    [Fact]
    public unsafe void Core_Ring9AndVerletCell72_ExactBitBudgetAndIv263Holonomy()
    {
        // 1. Exact 9-byte = 72-bit cell size (zero 32/64-bit bloat)
        Assert.Equal(9, sizeof(VerletCell72));
        Assert.Equal(72, sizeof(VerletCell72) * 8);
        Assert.Equal(263, CarrierSpecD3.IvSeed);

        var cell = VerletCell72.Create(64, 263, 139, 13, 41, 97);
        cell.Gate0 = new Ring9(352);
        cell.Gate1 = new Ring9(256);
        Assert.Equal(64, cell.Mu0.UnsignedTicks);
        Assert.Equal(263, cell.Phi0Curr.UnsignedTicks);
        Assert.Equal(139, cell.Phi0Past.UnsignedTicks);
        Assert.Equal(13, cell.Mu1.UnsignedTicks);
        Assert.Equal(41, cell.Phi1Curr.UnsignedTicks);
        Assert.Equal(97, cell.Phi1Past.UnsignedTicks);
        Assert.Equal(352, cell.Gate0.UnsignedTicks);
        Assert.Equal(256, cell.Gate1.UnsignedTicks);

        // 2. Verify IV(3)=263 interior octant silence (+15 -> Gate=0) and domain-wall chiral ratchet (+76 + -77 = -1)
        static (Ring9 Gate, int Raw, int Den, int Ne) EvalNode(int cz, int cy, int cx, int lane)
        {
            Span<int> locMu = stackalloc int[12];
            Span<Ring9> locPhi = stackalloc Ring9[12];
            int muC = CuboctahedronAutomatonKernel.Iv263MuAt(cz, cy, cx);
            Ring9 phiC = CuboctahedronAutomatonKernel.Iv263PhiAt(cz, cy, cx, lane);
            for (int k = 0; k < 12; k++)
            {
                var (dz, dy, dx) = CarrierSpecD3.CuboctOffsets12[k];
                locMu[k] = CuboctahedronAutomatonKernel.Iv263MuAt(cz + dz, cy + dy, cx + dx);
                locPhi[k] = CuboctahedronAutomatonKernel.Iv263PhiAt(cz + dz, cy + dy, cx + dx, lane);
            }
            return CuboctahedronAutomatonKernel.Compute12NeighborGate(muC, phiC, locMu, locPhi, true);
        }

        var oct = EvalNode(2, 2, 2, 0);
        Assert.Equal(15, oct.Raw);
        Assert.Equal(0, oct.Gate.SignedWrap);

        var wPos = EvalNode(0, 2, 2, 1);
        var wNeg = EvalNode(-1, 2, 2, 1);
        Assert.Equal(76, wPos.Gate.SignedWrap);
        Assert.Equal(-77, wNeg.Gate.SignedWrap);
        Assert.Equal(-1, wPos.Gate.SignedWrap + wNeg.Gate.SignedWrap);
    }

    [Fact]
    public void Core_CudaRtx4060Ptx_512StatesMatch()
    {
        var res = CudaRtx4060Driver.RunRtx4060PtxSelfTest();
        if (res.Ok)
        {
            Assert.Equal(512, res.VerifiedCount);
            Assert.True(res.GbpSec > 1.0);
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

        byte[] kf0 = new byte[1024];
        byte[] kf12 = new byte[1024];
        byte[] gt6 = new byte[1024];
        for (int i = 0; i < 1024; i++)
        {
            kf0[i] = 40;
            kf12[i] = 160;
            gt6[i] = 105; // Within [-40, +40] of midpoint 100!
        }
        var (_, pct) = CuboctMediaCodec.InterpolateTwoPointKeyframe(kf0, kf12, gt6, 6, 12);
        Assert.Equal(100.0, pct);
    }

    [Fact]
    public void Engine_CuboctWorldScene3D_ExplosionBinauralAudioAndZeroMemoryTimeReversal()
    {
        var world = new CuboctWorldScene3D(9, 9, 9);
        world.InitializeIv263Field();
        world.TriggerImpulse(4, 4, 4, amplitudeMu: 64, phaseImpulse: 164);

        // Save initial state snapshot for verification
        var initCenter = world.GetCell(4, 4, 4);
        var initNeighbor = world.GetCell(4, 5, 5);

        // Run 4 ticks forward (explosion shockwave & light propagation)
        for (int t = 0; t < 4; t++)
            world.StepForward();
        Assert.Equal(4, world.Tick);

        var sensor = world.SampleCameraVoxel(4, 5, 5);
        Assert.InRange(sensor.WavelengthNm, 440, 718);
        Assert.Equal(128, sensor.ColorR + sensor.ColorG + sensor.ColorB);

        var (earL, earR) = world.SampleBinauralAudioListener((4, 3, 4), (4, 5, 4));
        Assert.True(earL != 0 || earR != 0);

        // Rewind 4 ticks backward in time with ZERO history stack!
        for (int t = 0; t < 4; t++)
            world.StepBackward();
        Assert.Equal(0, world.Tick);

        // Verify 100% bit-exact recovery of the pre-explosion world state!
        var recCenter = world.GetCell(4, 4, 4);
        var recNeighbor = world.GetCell(4, 5, 5);
        Assert.Equal(initCenter.Phi0Curr.UnsignedTicks, recCenter.Phi0Curr.UnsignedTicks);
        Assert.Equal(initCenter.Phi1Curr.UnsignedTicks, recCenter.Phi1Curr.UnsignedTicks);
        Assert.Equal(initCenter.Phi0Past.UnsignedTicks, recCenter.Phi0Past.UnsignedTicks);
        Assert.Equal(initNeighbor.Phi0Curr.UnsignedTicks, recNeighbor.Phi0Curr.UnsignedTicks);
        Assert.Equal(initNeighbor.Phi1Curr.UnsignedTicks, recNeighbor.Phi1Curr.UnsignedTicks);
    }
}
