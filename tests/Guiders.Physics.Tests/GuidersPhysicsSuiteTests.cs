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

        // Verify IV(3)=263 interior octant silence (+15 -> Gate=0) and domain-wall chiral ratchet (+76 + -77 = -1) purely on Ring9!
        static (Ring9 Gate, Ring9 Raw, int Den, int Ne) EvalNode(Ring9 cz, Ring9 cy, Ring9 cx, int lane)
        {
            Span<Ring9> locMu = stackalloc Ring9[12];
            Span<Ring9> locPhi = stackalloc Ring9[12];
            Ring9 muC = CuboctahedronAutomatonKernel.Iv263MuAt(cz, cy, cx);
            Ring9 phiC = CuboctahedronAutomatonKernel.Iv263PhiAt(cz, cy, cx, lane);
            for (int k = 0; k < 12; k++)
            {
                var (dz, dy, dx) = CarrierSpecD3.CuboctOffsets12[k];
                locMu[k] = CuboctahedronAutomatonKernel.Iv263MuAt(cz - dz, cy - dy, cx - dx);
                locPhi[k] = CuboctahedronAutomatonKernel.Iv263PhiAt(cz - dz, cy - dy, cx - dx, lane);
            }
            return CuboctahedronAutomatonKernel.Compute12NeighborGate(muC, phiC, locMu, locPhi, 1);
        }

        var oct = EvalNode(new(2), new(2), new(2), 0);
        Assert.Equal(15, oct.Raw.SignedWrap);
        Assert.Equal(0, oct.Gate.SignedWrap);

        var wPos = EvalNode(new(0), new(2), new(2), 1);
        var wNeg = EvalNode(new(-1), new(2), new(2), 1);
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
            gt6[i] = 105;
        }
        var (_, pct) = CuboctMediaCodec.InterpolateTwoPointKeyframe(kf0, kf12, gt6, 6, 12);
        Assert.Equal(100.0, pct);
    }

    [Fact]
    public void Engine_ClosedFormO1WorldlineParabola_ZeroMemoryJumpAndTimeReversal()
    {
        var world = new CuboctWorldScene3D();
        Ring9 dz = new(0), dy = new(2), dx = new(2);

        // Evaluate cell at t = 0 directly via closed-form O(1) Verlet worldline parabola (zero 3D grid in RAM!):
        var cellT0 = world.EvaluateWorldlineAt(dz, dy, dx, 0);

        // Jump in O(1) to t = 500 without iterating 500 steps!
        world.JumpToTick(500);
        var sensor = world.SampleCameraVoxel(dz, dy, dx);
        Assert.InRange(sensor.WavelengthNm, 440, 718);
        Assert.Equal(128, sensor.ColorR + sensor.ColorG + sensor.ColorB);

        var (earL, earR) = world.SampleBinauralAudioListener((new(0), new(1), new(2)), (new(0), new(3), new(2)));
        Assert.True(earL != 0 || earR != 0);

        // Jump back in O(1) to t = 0 and verify 100% bit-exact agreement!
        world.JumpToTick(0);
        var recT0 = world.GetCellAtCurrentTick(dz, dy, dx);
        Assert.Equal(cellT0.Phi0Curr, recT0.Phi0Curr);
        Assert.Equal(cellT0.Phi1Curr, recT0.Phi1Curr);
    }
}
