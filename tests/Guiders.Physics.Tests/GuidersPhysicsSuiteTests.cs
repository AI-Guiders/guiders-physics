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
    public unsafe void Core_Ring9AndVerletCell72_Pure9BitAlgebraAndBezout()
    {
        Assert.Equal(9, sizeof(VerletCell72));
        Assert.Equal(72, sizeof(VerletCell72) * 8);
        Assert.Equal(new Ring9(263), Ring9.IvSeed);

        var cell = VerletCell72.Create(new(64), new(263), new(139), new(13), new(41), new(97));
        cell.Gate0 = new Ring9(352);
        cell.Gate1 = new Ring9(256);
        Assert.Equal(new Ring9(64), cell.Mu0);
        Assert.Equal(new Ring9(263), cell.Phi0Curr);
        Assert.Equal(new Ring9(139), cell.Phi0Past);
        Assert.Equal(new Ring9(13), cell.Mu1);
        Assert.Equal(new Ring9(41), cell.Phi1Curr);
        Assert.Equal(new Ring9(97), cell.Phi1Past);
        Assert.Equal(new Ring9(352), cell.Gate0);
        Assert.Equal(new Ring9(256), cell.Gate1);

        for (ushort x = 0; x < 512; x++)
        {
            Ring9 rx = new(x);
            Assert.Equal(rx, CuboctahedronAutomatonKernel.DecodeAtomicTimeTick(rx * Ring9.DPhiDisc));
            Assert.Equal(Ring9.Zero, CuboctahedronAutomatonKernel.AtomicTimeServerPfdError(rx, rx));
            Assert.Equal(Ring9.CollapsePhi, CuboctahedronAutomatonKernel.AtomicRelativePhaseAtTarget(rx, Ring9.CollapsePhi));
            Assert.Equal(Ring9.IvSeed, CuboctahedronAutomatonKernel.TwoPointWorldlineInterpolate(Ring9.IvSeed, rx, Ring9.DPhiDisc, Ring9.Zero));
            Assert.Equal(rx, CuboctahedronAutomatonKernel.TwoPointWorldlineInterpolate(Ring9.IvSeed, rx, Ring9.DPhiDisc, Ring9.MuVac));
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
    public void Engine_ClosedFormO1WorldlineParabola_PureRing9JumpAndTimeReversal()
    {
        var world = new CuboctWorldScene3D();
        Ring9 dz = new(0), dy = new(2), dx = new(2);

        var cellT0 = world.EvaluateWorldlineAt(dz, dy, dx, Ring9.Zero);

        world.JumpToTick(new Ring9(500));
        var sensor = world.SampleCameraVoxel(dz, dy, dx);
        Assert.Equal(Ring9.QSquare, sensor.ColorR + sensor.ColorG + sensor.ColorB);
        Assert.Equal(new Ring9(26), sensor.OlfactoryShellZ);

        var (earL, earR) = world.SampleBinauralAudioListener((new(0), new(1), new(2)), (new(0), new(3), new(2)));
        Assert.NotEqual(Ring9.Zero, earL + earR);

        world.JumpToTick(Ring9.Zero);
        var recT0 = world.GetCellAtCurrentTick(dz, dy, dx);
        Assert.Equal(cellT0.Phi0Curr, recT0.Phi0Curr);
        Assert.Equal(cellT0.Phi1Curr, recT0.Phi1Curr);

        for (ushort tau = 0; tau < 512; tau++)
        {
            Ring9 targetTick = new(tau);
            Ring9 dphiVerlet = targetTick * Ring9.DPhiDisc;
            Ring9 atomPhase = Ring9.IvSeed + dphiVerlet;
            world.SyncWithPhysicalAtom(atomPhase, dphiVerlet);
            Assert.Equal(targetTick, world.CurrentTick);
        }

        world.ClockMode = WorldClockMode.FreeRunning;
        world.JumpToTick(new Ring9(10));
        world.AdvanceClock();
        Assert.Equal(new Ring9(11), world.CurrentTick);

        var beforeExplosion = world.GetCellAtCurrentTick(dz, dy, dx);
        world.EmitGameEvent(dz, dy, dx, Ring9.KP, Ring9.CollapsePhi);
        var afterExplosion = world.GetCellAtCurrentTick(dz, dy, dx);
        Assert.NotEqual(beforeExplosion.Phi0Curr, afterExplosion.Phi0Curr);

        world.ClockMode = WorldClockMode.AtomSynced;
        Ring9 extTick = new(77);
        world.AdvanceClock(Ring9.IvSeed + extTick * Ring9.DPhiDisc, extTick * Ring9.DPhiDisc);
        Assert.Equal(extTick, world.CurrentTick);
    }
}
