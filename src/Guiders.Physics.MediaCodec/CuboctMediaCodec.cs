using System;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using Guiders.Physics.Compression;
using Guiders.Physics.Core;

namespace Guiders.Physics.MediaCodec;

/// <summary>
/// 4D Video and Studio Audio Streaming Codec based on 2-Point Verlet Parabola Prediction,
/// Sub-Green Color Decorrelation, and Two-Point Keyframe Straight-Line Interpolation in the Heisenberg [-40, +40] Deadband.
/// Uses the pure native CuboctUniversalArchiver (zero System.IO.Compression dependencies).
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
        short[] verlet = new short[nStereo * 2];
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

        byte[] rawBytes = new byte[verlet.Length * 2];
        Buffer.BlockCopy(verlet, 0, rawBytes, 0, rawBytes.Length);
        return CuboctUniversalArchiver.CompressBlock(rawBytes);
    }

    /// <summary>
    /// Decodes a 2-Point Verlet encoded 16-bit Stereo PCM audio buffer 100% bit-exact.
    /// </summary>
    public static short[] DecodeStereoPcm16Verlet(ReadOnlySpan<byte> compressedAudio, int nStereoSamples)
    {
        int byteLen = nStereoSamples * 4;
        byte[] recBytes = CuboctUniversalArchiver.DecompressBlock(compressedAudio, byteLen);

        short[] recVerlet = new short[nStereoSamples * 2];
        Buffer.BlockCopy(recBytes, 0, recVerlet, 0, byteLen);
        short[] pcm = new short[nStereoSamples * 2];
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
        int limit = Ring9.DPhiDisc.Bits9; // 41

        for (int i = 0; i < len; i++)
        {
            int p0 = keyframe0[i];
            int p1 = keyframeN[i];
            int pred = p0 + (tOffset * (p1 - p0)) / strideN;
            rec[i] = (byte)Math.Clamp(pred, 0, 255);
            if (!groundTruthFrameT.IsEmpty && Math.Abs(groundTruthFrameT[i] - pred) < limit)
                inside++;
        }

        double pct = groundTruthFrameT.IsEmpty ? 100.0 : (inside * 100.0 / len);
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
