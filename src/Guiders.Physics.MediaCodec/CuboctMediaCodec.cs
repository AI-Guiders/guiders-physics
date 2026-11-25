using System;
using System.Buffers.Binary;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using Guiders.Physics.Compression;
using Guiders.Physics.Core;

namespace Guiders.Physics.MediaCodec;

/// <summary>
/// 4D Video and Studio Audio Streaming Codec built directly on <see cref="Guiders.Physics.Core"/>
/// (<see cref="CuboctahedronAutomatonKernel"/>) and <see cref="CuboctUniversalArchiver"/>:
///   - Sub-Green reversible color decorrelation (`R - G`, `G`, `B - G`),
///   - Intra-frame (`stride = 3`) and inter-frame (`stride = frameBytes`) SIMD Verlet prediction via
///     <see cref="CuboctahedronAutomatonKernel.ComputeVerletResidualBuffer"/> and
///     <see cref="CuboctahedronAutomatonKernel.ReconstructVerletBuffer"/>,
///   - 2-Point Verlet Stereo PCM16 audio encoding,
///   - Multiplexed `.cuboct` Audio+Video container with 100% bit-exact lossless recovery.
/// </summary>
public static class CuboctMediaCodec
{
    /// <summary>
    /// Encodes a 16-bit Stereo PCM audio buffer losslessly using Mid/Side (L, R-L) decorrelation
    /// and 2nd-order Verlet parabola prediction `s_hat(t) = 2*s(t-1) - s(t-2)` + native Heisenberg-41 compression.
    /// </summary>
    public static byte[] EncodeStereoPcm16Verlet(ReadOnlySpan<short> stereoInterleavedSamples)
    {
        int nStereo = stereoInterleavedSamples.Length / 2;
        short[] verlet = GC.AllocateUninitializedArray<short>(nStereo * 2);
        short prevL1 = 0, prevS1 = 0, prevD1L = 0, prevD1S = 0;

        for (int i = 0; i < nStereo; i++)
        {
            short l = stereoInterleavedSamples[2 * i];
            short r = stereoInterleavedSamples[2 * i + 1];
            short s = (short)(r - l);
            short d1L = (short)(l - prevL1);
            short d1S = (short)(s - prevS1);
            verlet[2 * i] = (short)(d1L - prevD1L);
            verlet[2 * i + 1] = (short)(d1S - prevD1S);
            prevL1 = l; prevS1 = s;
            prevD1L = d1L; prevD1S = d1S;
        }

        byte[] rawBytes = GC.AllocateUninitializedArray<byte>(verlet.Length * 2);
        Buffer.BlockCopy(verlet, 0, rawBytes, 0, rawBytes.Length);
        return CuboctUniversalArchiver.CompressPredecorrelatedResiduals(rawBytes);
    }

    /// <summary>
    /// Decodes a 2-Point Verlet encoded 16-bit Stereo PCM audio buffer 100% bit-exact.
    /// </summary>
    public static short[] DecodeStereoPcm16Verlet(ReadOnlySpan<byte> compressedAudio, int nStereoSamples)
    {
        int byteLen = nStereoSamples * 4;
        byte[] recBytes = CuboctUniversalArchiver.DecompressBlock(compressedAudio, byteLen);

        short[] recVerlet = GC.AllocateUninitializedArray<short>(nStereoSamples * 2);
        Buffer.BlockCopy(recBytes, 0, recVerlet, 0, byteLen);
        short[] pcm = GC.AllocateUninitializedArray<short>(nStereoSamples * 2);
        short rL1 = 0, rS1 = 0, rD1L = 0, rD1S = 0;

        for (int i = 0; i < nStereoSamples; i++)
        {
            rD1L = (short)(rD1L + recVerlet[2 * i]);
            rD1S = (short)(rD1S + recVerlet[2 * i + 1]);
            rL1 = (short)(rL1 + rD1L);
            rS1 = (short)(rS1 + rD1S);
            pcm[2 * i] = rL1;
            pcm[2 * i + 1] = (short)(rL1 + rS1);
        }
        return pcm;
    }

    /// <summary>
    /// Encodes a 4D RGB24 video sequence (`frameCount x height x width x 3`) losslessly using Sub-Green color decorrelation,
    /// kernel intra/inter-frame Verlet prediction (<see cref="CuboctahedronAutomatonKernel.ComputeVerletResidualBuffer"/>),
    /// and single-pass Heisenberg-41 compression.
    /// </summary>
    public static byte[] EncodeRgb24VideoVerlet(ReadOnlySpan<byte> rgb24Frames, int frameCount, int width, int height)
    {
        int frameBytes = width * height * 3;
        int totalBytes = frameCount * frameBytes;
        byte[] subGreen = GC.AllocateUninitializedArray<byte>(totalBytes);

        // 1. Sub-Green color decorrelation: (R - G, G, B - G)
        for (int i = 0; i <= totalBytes - 3; i += 3)
        {
            byte g = rgb24Frames[i + 1];
            subGreen[i]     = (byte)(rgb24Frames[i] - g);
            subGreen[i + 1] = g;
            subGreen[i + 2] = (byte)(rgb24Frames[i + 2] - g);
        }

        byte[] residuals = GC.AllocateUninitializedArray<byte>(totalBytes);

        // 2. Inter-frame 1st-order Verlet residual across all frames via Guiders.Physics.Core kernel
        CuboctahedronAutomatonKernel.ComputeVerletResidualBuffer(
            subGreen, residuals, frameBytes, 0xFF, 0x00);

        // 3. Intra-frame pixel-stride (stride = 3) Verlet residual on Frame 0 via Guiders.Physics.Core kernel
        CuboctahedronAutomatonKernel.ComputeVerletResidualBuffer(
            subGreen.AsSpan(0, frameBytes), residuals.AsSpan(0, frameBytes), 3, 0xFF, 0x00);

        return CuboctUniversalArchiver.CompressPredecorrelatedResiduals(residuals);
    }

    /// <summary>
    /// Decodes a 4D RGB24 video sequence (`frameCount x height x width x 3`) 100% bit-exact via
    /// <see cref="CuboctahedronAutomatonKernel.ReconstructVerletBuffer"/> and inverse Sub-Green color transform.
    /// </summary>
    public static byte[] DecodeRgb24VideoVerlet(ReadOnlySpan<byte> compressedVideo, int frameCount, int width, int height)
    {
        int frameBytes = width * height * 3;
        int totalBytes = frameCount * frameBytes;
        byte[] residuals = CuboctUniversalArchiver.DecompressBlock(compressedVideo, totalBytes);
        byte[] subGreen = GC.AllocateUninitializedArray<byte>(totalBytes);

        // 1. Reconstruct Frame 0 from intra-frame pixel-stride (stride = 3) residuals via Guiders.Physics.Core kernel
        CuboctahedronAutomatonKernel.ReconstructVerletBuffer(
            residuals.AsSpan(0, frameBytes), residuals.AsSpan(0, frameBytes), 3, 0xFF, 0x00);

        // 2. Reconstruct all frames 0..frameCount-1 from inter-frame (stride = frameBytes) residuals via Guiders.Physics.Core kernel
        CuboctahedronAutomatonKernel.ReconstructVerletBuffer(
            residuals, subGreen, frameBytes, 0xFF, 0x00);

        // 3. Inverse Sub-Green color transform: (dR + G, G, dB + G)
        for (int i = 0; i <= totalBytes - 3; i += 3)
        {
            byte g = subGreen[i + 1];
            subGreen[i]     = (byte)(subGreen[i] + g);
            subGreen[i + 2] = (byte)(subGreen[i + 2] + g);
        }

        return subGreen;
    }

    /// <summary>
    /// Multiplexes and compresses a full RGB24 video stream and Stereo PCM16 audio track into a single `.cuboct` container.
    /// </summary>
    public static byte[] EncodeMultiplexedMediaContainer(
        ReadOnlySpan<byte> rgb24Frames, int frameCount, int width, int height,
        ReadOnlySpan<short> stereoInterleavedPcm16, int sampleRate = 48000)
    {
        byte[] compVideo = EncodeRgb24VideoVerlet(rgb24Frames, frameCount, width, height);
        int nStereoSamples = stereoInterleavedPcm16.Length / 2;
        byte[] compAudio = EncodeStereoPcm16Verlet(stereoInterleavedPcm16);

        // Header: 28 bytes (Width, Height, FrameCount, SampleRate, StereoSamples, VideoLen, AudioLen)
        byte[] container = GC.AllocateUninitializedArray<byte>(28 + compVideo.Length + compAudio.Length);
        Span<byte> hdr = container.AsSpan(0, 28);
        BinaryPrimitives.WriteInt32LittleEndian(hdr.Slice(0, 4), width);
        BinaryPrimitives.WriteInt32LittleEndian(hdr.Slice(4, 4), height);
        BinaryPrimitives.WriteInt32LittleEndian(hdr.Slice(8, 4), frameCount);
        BinaryPrimitives.WriteInt32LittleEndian(hdr.Slice(12, 4), sampleRate);
        BinaryPrimitives.WriteInt32LittleEndian(hdr.Slice(16, 4), nStereoSamples);
        BinaryPrimitives.WriteInt32LittleEndian(hdr.Slice(20, 4), compVideo.Length);
        BinaryPrimitives.WriteInt32LittleEndian(hdr.Slice(24, 4), compAudio.Length);

        compVideo.CopyTo(container.AsSpan(28));
        compAudio.CopyTo(container.AsSpan(28 + compVideo.Length));
        return container;
    }

    /// <summary>
    /// Demultiplexes and losslessly decodes both the RGB24 video frames and Stereo PCM16 audio track from a `.cuboct` container.
    /// </summary>
    public static (byte[] Rgb24Frames, short[] StereoPcm16, int Width, int Height, int FrameCount, int SampleRate)
        DecodeMultiplexedMediaContainer(ReadOnlySpan<byte> container)
    {
        int width = BinaryPrimitives.ReadInt32LittleEndian(container.Slice(0, 4));
        int height = BinaryPrimitives.ReadInt32LittleEndian(container.Slice(4, 4));
        int frameCount = BinaryPrimitives.ReadInt32LittleEndian(container.Slice(8, 4));
        int sampleRate = BinaryPrimitives.ReadInt32LittleEndian(container.Slice(12, 4));
        int nStereoSamples = BinaryPrimitives.ReadInt32LittleEndian(container.Slice(16, 4));
        int videoLen = BinaryPrimitives.ReadInt32LittleEndian(container.Slice(20, 4));
        int audioLen = BinaryPrimitives.ReadInt32LittleEndian(container.Slice(24, 4));

        ReadOnlySpan<byte> compVideo = container.Slice(28, videoLen);
        ReadOnlySpan<byte> compAudio = container.Slice(28 + videoLen, audioLen);

        byte[] rgb24 = DecodeRgb24VideoVerlet(compVideo, frameCount, width, height);
        short[] pcm16 = DecodeStereoPcm16Verlet(compAudio, nStereoSamples);
        return (rgb24, pcm16, width, height, frameCount, sampleRate);
    }

    /// <summary>
    /// Reconstructs dropped intermediate RGB24 video frames `t = 1..stride-1` between two keyframes `frame0` and `frameN`
    /// along the straight line through two points: `I_hat(t) = I(0) + (t * (I(N) - I(0))) / stride`,
    /// and measures the fraction of pixels lying inside the Heisenberg deadband `[-40, +40]`.
    /// </summary>
    public static (byte[] ReconstructedFrame, double HeisenbergInsidePercent) InterpolateTwoPointKeyframe(
        ReadOnlySpan<byte> keyframe0, ReadOnlySpan<byte> keyframeN, ReadOnlySpan<byte> groundTruthFrameT, int tOffset, int strideN)
    {
        int len = keyframe0.Length;
        byte[] rec = new byte[len];
        long inside = 0;

        for (int i = 0; i < len; i++)
        {
            int p0 = keyframe0[i];
            int p1 = keyframeN[i];
            int pred = p0 + (tOffset * (p1 - p0)) / strideN;
            rec[i] = (byte)Math.Clamp(pred, 0, 255);
        }

        int gtLen = Math.Min(len, groundTruthFrameT.Length);
        for (int i = 0; i < gtLen; i++)
        {
            Ring9 err9 = new((ushort)(512 + groundTruthFrameT[i] - rec[i]));
            ushort ne = CuboctahedronAutomatonKernel.HeisenbergEnergyQuanta(err9).Bits9;
            inside += 1 / (ne + 1);
        }

        int hasGt = ((0 - gtLen) >> 31) & 1;
        int safeDen = gtLen + (1 - hasGt);
        double pct = (1 - hasGt) * 100.0 + hasGt * (inside * 100.0 / safeDen);
        return (rec, pct);
    }

    /// <summary>
    /// SIMD Vector256 2nd-order Verlet frame residual transform: `res[i] = curr[i] - (2*prev1[i] - prev2[i])`.
    /// </summary>
    public static void ComputeSimdVerletFrameResidual(
        ReadOnlySpan<byte> prev2, ReadOnlySpan<byte> prev1, ReadOnlySpan<byte> curr, Span<byte> residualOut)
    {
        int len = curr.Length;
        int vecCount = len / 32;
        ref byte rP2 = ref MemoryMarshal.GetReference(prev2);
        ref byte rP1 = ref MemoryMarshal.GetReference(prev1);
        ref byte rCur = ref MemoryMarshal.GetReference(curr);
        ref byte rOut = ref MemoryMarshal.GetReference(residualOut);

        for (int i = 0; i < vecCount; i++)
        {
            nuint off = (nuint)(i * 32);
            var vP2 = Vector256.LoadUnsafe(ref rP2, off);
            var vP1 = Vector256.LoadUnsafe(ref rP1, off);
            var vCur = Vector256.LoadUnsafe(ref rCur, off);
            var vPred = (vP1 + vP1) - vP2;
            (vCur - vPred).StoreUnsafe(ref rOut, off);
        }
        for (int i = vecCount * 32; i < len; i++)
        {
            residualOut[i] = (byte)(curr[i] - (2 * prev1[i] - prev2[i]));
        }
    }
}
