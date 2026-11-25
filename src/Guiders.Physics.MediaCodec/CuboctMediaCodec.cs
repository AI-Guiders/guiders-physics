using System;
using System.Buffers;
using System.Buffers.Binary;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using Guiders.Physics.Compression;
using Guiders.Physics.Core;

namespace Guiders.Physics.MediaCodec;

/// <summary>
/// Zero-allocation 4D Video and Studio Audio Streaming Codec built directly on <see cref="Guiders.Physics.Core"/>
/// (<see cref="CuboctahedronAutomatonKernel"/>) and <see cref="CuboctUniversalArchiver"/>:
///   - Zero-allocation <see cref="EncodeStereoPcm16VerletSpan"/>, <see cref="DecodeStereoPcm16VerletSpan"/> (via <see cref="MemoryMarshal.AsBytes{T}(Span{T})"/>),
///   - Zero-allocation <see cref="EncodeRgb24VideoVerletSpan"/> and 100% in-place <see cref="DecodeRgb24VideoVerletSpan"/>
///     (requires zero scratch buffers and allocates `0` bytes on the managed heap!),
///   - Sub-Green reversible color decorrelation (`R - G`, `G`, `B - G`),
///   - Intra-frame (`stride = 3`) and inter-frame (`stride = frameBytes`) SIMD Verlet prediction via
///     <see cref="CuboctUniversalArchiver.ComputeVerletResidualBuffer"/> and
///     <see cref="CuboctUniversalArchiver.ReconstructVerletBuffer"/>.
/// </summary>
public static class CuboctMediaCodec
{
    /// <summary>
    /// Zero-allocation Stereo PCM16 encoder (`0` heap allocations): writes Mid/Side 2-Point Verlet residuals
    /// into <paramref name="verletScratch"/> (length &gt;= <c>stereoInterleavedSamples.Length</c>) and compresses
    /// directly into <paramref name="destination"/> via <see cref="MemoryMarshal.AsBytes{T}(ReadOnlySpan{T})"/>.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int EncodeStereoPcm16VerletSpan(
        ReadOnlySpan<short> stereoInterleavedSamples, Span<short> verletScratch, Span<byte> destination)
    {
        int nStereo = stereoInterleavedSamples.Length / 2;
        Span<short> verlet = verletScratch.Slice(0, nStereo * 2);
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

        ReadOnlySpan<byte> rawBytes = MemoryMarshal.AsBytes((ReadOnlySpan<short>)verlet);
        return CuboctUniversalArchiver.CompressPredecorrelatedSpan(rawBytes, destination);
    }

    /// <summary>
    /// Zero-allocation in-place Stereo PCM16 decoder (`0` heap allocations, zero scratch buffer required):
    /// decodes directly into <paramref name="destinationPcm16"/> (length == <c>nStereoSamples * 2</c>) via
    /// <see cref="MemoryMarshal.AsBytes{T}(Span{T})"/> and reconstructs Mid/Side 2-Point Verlet in-place.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int DecodeStereoPcm16VerletSpan(
        ReadOnlySpan<byte> compressedAudio, int nStereoSamples, Span<short> destinationPcm16)
    {
        Span<short> pcm = destinationPcm16.Slice(0, nStereoSamples * 2);
        Span<byte> pcmBytes = MemoryMarshal.AsBytes(pcm);
        CuboctUniversalArchiver.DecompressSpan(compressedAudio, pcmBytes);

        short rL1 = 0, rS1 = 0, rD1L = 0, rD1S = 0;
        for (int i = 0; i < nStereoSamples; i++)
        {
            rD1L = (short)(rD1L + pcm[2 * i]);
            rD1S = (short)(rD1S + pcm[2 * i + 1]);
            rL1 = (short)(rL1 + rD1L);
            rS1 = (short)(rS1 + rD1S);
            pcm[2 * i] = rL1;
            pcm[2 * i + 1] = (short)(rL1 + rS1);
        }
        return nStereoSamples * 2;
    }

    public static byte[] EncodeStereoPcm16Verlet(ReadOnlySpan<short> stereoInterleavedSamples)
    {
        int len = stereoInterleavedSamples.Length;
        short[] rentedVerlet = ArrayPool<short>.Shared.Rent(Math.Max(16, len));
        byte[] rentedOut = ArrayPool<byte>.Shared.Rent(CuboctUniversalArchiver.GetMaxCompressedLength(len * 2));
        try
        {
            int written = EncodeStereoPcm16VerletSpan(stereoInterleavedSamples, rentedVerlet, rentedOut);
            return rentedOut.AsSpan(0, written).ToArray();
        }
        finally
        {
            ArrayPool<short>.Shared.Return(rentedVerlet);
            ArrayPool<byte>.Shared.Return(rentedOut);
        }
    }

    public static short[] DecodeStereoPcm16Verlet(ReadOnlySpan<byte> compressedAudio, int nStereoSamples)
    {
        short[] pcm = GC.AllocateUninitializedArray<short>(nStereoSamples * 2);
        DecodeStereoPcm16VerletSpan(compressedAudio, nStereoSamples, pcm);
        return pcm;
    }

    /// <summary>
    /// Zero-allocation 4D RGB24 video encoder (`0` heap allocations): uses <paramref name="subGreenScratch"/>
    /// and <paramref name="residualScratch"/> (each length &gt;= <c>frameCount * width * height * 3</c>) and writes
    /// compressed video directly into <paramref name="destination"/>.
    /// </summary>
    public static int EncodeRgb24VideoVerletSpan(
        ReadOnlySpan<byte> rgb24Frames, int frameCount, int width, int height,
        Span<byte> destination, Span<byte> subGreenScratch, Span<byte> residualScratch)
    {
        int frameBytes = width * height * 3;
        int totalBytes = frameCount * frameBytes;
        Span<byte> subGreen = subGreenScratch.Slice(0, totalBytes);
        Span<byte> residuals = residualScratch.Slice(0, totalBytes);

        for (int i = 0; i <= totalBytes - 3; i += 3)
        {
            byte g = rgb24Frames[i + 1];
            subGreen[i]     = (byte)(rgb24Frames[i] - g);
            subGreen[i + 1] = g;
            subGreen[i + 2] = (byte)(rgb24Frames[i + 2] - g);
        }

        CuboctUniversalArchiver.ComputeVerletResidualBuffer(
            subGreen, residuals, frameBytes, 0xFF, 0x00);

        CuboctUniversalArchiver.ComputeVerletResidualBuffer(
            subGreen.Slice(0, frameBytes), residuals.Slice(0, frameBytes), 3, 0xFF, 0x00);

        return CuboctUniversalArchiver.CompressPredecorrelatedSpan(residuals, destination);
    }

    /// <summary>
    /// Zero-allocation 100% in-place 4D RGB24 video decoder (`0` heap allocations, zero scratch buffer required!):
    /// decodes tokens, intra-frame Verlet, inter-frame Verlet, and Sub-Green color directly in-place inside <paramref name="destinationRgb24"/>.
    /// </summary>
    public static int DecodeRgb24VideoVerletSpan(
        ReadOnlySpan<byte> compressedVideo, int frameCount, int width, int height, Span<byte> destinationRgb24)
    {
        int frameBytes = width * height * 3;
        int totalBytes = frameCount * frameBytes;
        Span<byte> dst = destinationRgb24.Slice(0, totalBytes);

        // 1. Decode tokens directly into `dst` (0 allocations)
        CuboctUniversalArchiver.DecompressSpan(compressedVideo, dst);

        // 2. Reconstruct Frame 0 in-place from intra-frame pixel-stride (stride = 3) residuals
        CuboctUniversalArchiver.ReconstructVerletBuffer(
            dst.Slice(0, frameBytes), dst.Slice(0, frameBytes), 3, 0xFF, 0x00);

        // 3. Reconstruct all frames 0..frameCount-1 in-place from inter-frame (stride = frameBytes) residuals
        CuboctUniversalArchiver.ReconstructVerletBuffer(
            dst, dst, frameBytes, 0xFF, 0x00);

        // 4. Inverse Sub-Green color transform in-place: (dR + G, G, dB + G)
        for (int i = 0; i <= totalBytes - 3; i += 3)
        {
            byte g = dst[i + 1];
            dst[i]     = (byte)(dst[i] + g);
            dst[i + 2] = (byte)(dst[i + 2] + g);
        }

        return totalBytes;
    }

    public static byte[] EncodeRgb24VideoVerlet(ReadOnlySpan<byte> rgb24Frames, int frameCount, int width, int height)
    {
        int totalBytes = frameCount * width * height * 3;
        byte[] rentedSubG = ArrayPool<byte>.Shared.Rent(totalBytes);
        byte[] rentedRes  = ArrayPool<byte>.Shared.Rent(totalBytes);
        byte[] rentedOut  = ArrayPool<byte>.Shared.Rent(CuboctUniversalArchiver.GetMaxCompressedLength(totalBytes));
        try
        {
            int written = EncodeRgb24VideoVerletSpan(rgb24Frames, frameCount, width, height, rentedOut, rentedSubG, rentedRes);
            return rentedOut.AsSpan(0, written).ToArray();
        }
        finally
        {
            ArrayPool<byte>.Shared.Return(rentedSubG);
            ArrayPool<byte>.Shared.Return(rentedRes);
            ArrayPool<byte>.Shared.Return(rentedOut);
        }
    }

    public static byte[] DecodeRgb24VideoVerlet(ReadOnlySpan<byte> compressedVideo, int frameCount, int width, int height)
    {
        int totalBytes = frameCount * width * height * 3;
        byte[] dst = GC.AllocateUninitializedArray<byte>(totalBytes);
        DecodeRgb24VideoVerletSpan(compressedVideo, frameCount, width, height, dst);
        return dst;
    }

    public static byte[] EncodeMultiplexedMediaContainer(
        ReadOnlySpan<byte> rgb24Frames, int frameCount, int width, int height,
        ReadOnlySpan<short> stereoInterleavedPcm16, int sampleRate = 48000)
    {
        int totalVideoBytes = frameCount * width * height * 3;
        int nStereoSamples = stereoInterleavedPcm16.Length / 2;
        int maxVidComp = CuboctUniversalArchiver.GetMaxCompressedLength(totalVideoBytes);
        int maxAudComp = CuboctUniversalArchiver.GetMaxCompressedLength(stereoInterleavedPcm16.Length * 2);

        byte[] rentedSubG = ArrayPool<byte>.Shared.Rent(totalVideoBytes);
        byte[] rentedRes  = ArrayPool<byte>.Shared.Rent(totalVideoBytes);
        short[] rentedAud = ArrayPool<short>.Shared.Rent(Math.Max(16, stereoInterleavedPcm16.Length));
        byte[] rentedContainer = ArrayPool<byte>.Shared.Rent(28 + maxVidComp + maxAudComp);
        try
        {
            int vidLen = EncodeRgb24VideoVerletSpan(
                rgb24Frames, frameCount, width, height, rentedContainer.AsSpan(28, maxVidComp), rentedSubG, rentedRes);
            int audLen = EncodeStereoPcm16VerletSpan(
                stereoInterleavedPcm16, rentedAud, rentedContainer.AsSpan(28 + vidLen, maxAudComp));

            Span<byte> hdr = rentedContainer.AsSpan(0, 28);
            BinaryPrimitives.WriteInt32LittleEndian(hdr.Slice(0, 4), width);
            BinaryPrimitives.WriteInt32LittleEndian(hdr.Slice(4, 4), height);
            BinaryPrimitives.WriteInt32LittleEndian(hdr.Slice(8, 4), frameCount);
            BinaryPrimitives.WriteInt32LittleEndian(hdr.Slice(12, 4), sampleRate);
            BinaryPrimitives.WriteInt32LittleEndian(hdr.Slice(16, 4), nStereoSamples);
            BinaryPrimitives.WriteInt32LittleEndian(hdr.Slice(20, 4), vidLen);
            BinaryPrimitives.WriteInt32LittleEndian(hdr.Slice(24, 4), audLen);

            return rentedContainer.AsSpan(0, 28 + vidLen + audLen).ToArray();
        }
        finally
        {
            ArrayPool<byte>.Shared.Return(rentedSubG);
            ArrayPool<byte>.Shared.Return(rentedRes);
            ArrayPool<short>.Shared.Return(rentedAud);
            ArrayPool<byte>.Shared.Return(rentedContainer);
        }
    }

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
