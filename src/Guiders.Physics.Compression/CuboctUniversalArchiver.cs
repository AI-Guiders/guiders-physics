using System;
using System.Buffers;
using System.Buffers.Binary;
using System.IO;
using System.Numerics;
using System.Runtime.CompilerServices;
using Guiders.Physics.Core;

namespace Guiders.Physics.Compression;

/// <summary>
/// Pure native branchless (zero conditional branching statements) and zero-allocation lossless archiver
/// built directly on <see cref="Guiders.Physics.Core"/> (<see cref="Ring9"/>, <see cref="VerletCell72"/> 8-lane architecture,
/// and <see cref="CuboctahedronAutomatonKernel"/>):
///   - Zero-allocation <see cref="CompressSpan"/>, <see cref="CompressPredecorrelatedSpan"/>, and in-place <see cref="DecompressSpan"/>
///     (`0` bytes allocated on the managed heap per block),
///   - Multi-gigabyte (>14 GB) streaming pipeline <see cref="CompressStream"/> and <see cref="DecompressStream"/> backed by
///     a single reusable <see cref="ArrayPool{T}"/> buffer (constant RAM footprint, zero per-chunk GC allocations),
///   - 64-bit SWAR (8-Lane <see cref="VerletCell72"/> register) quantum bit-width packer (`b in 0..8` bits/lane + O(1) 8-byte ring match `code = 9`).
/// </summary>
public static class CuboctUniversalArchiver
{
    public static readonly int DefaultStride18 = Ring9.D.Bits9 * Ring9.NSq.Bits9;             // 3 * 6 = 18
    public static readonly int HeisenbergMod41 = Ring9.DPhiDisc.Bits9;                        // 41
    public static readonly int HeisenbergLimit = Ring9.DPhiDisc.Bits9 - Ring9.MuVac.Bits9;    // 40

    private const int HashBits = 12;
    private const int HashSize = 1 << HashBits; // 4096-entry stack-allocated 8-byte ring hash table
    private const ulong HashMultiplier64 = 11400714819323198485UL;

    private static readonly int[] ModeStrides = [1, 1, 4, 18, 18];
    private static readonly byte[] ModeMask1  = [0x00, 0xFF, 0xFF, 0xFF, 0xFF];
    private static readonly byte[] ModeMask2  = [0x00, 0x00, 0x00, 0x00, 0xFF];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int GetMaxCompressedLength(int uncompressedLength)
        => uncompressedLength + (uncompressedLength >> 3) + 64;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static byte ZigZagEncodeSignedByte(sbyte delta)
        => (byte)((delta << 1) ^ (delta >> 7));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static sbyte ZigZagDecodeSignedByte(byte u)
        => (sbyte)((u >> 1) ^ -(u & 1));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static ulong ZigZagEncode8LanesSwar(ulong raw8)
    {
        ulong shl1 = (raw8 & 0x7F7F7F7F7F7F7F7FUL) << 1;
        ulong signMask = ((raw8 >> 7) & 0x0101010101010101UL) * 0xFFUL;
        return shl1 ^ signMask;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static ulong ZigZagDecode8LanesSwar(ulong u8)
    {
        ulong half = (u8 >> 1) & 0x7F7F7F7F7F7F7F7FUL;
        ulong signMask = (u8 & 0x0101010101010101UL) * 0xFFUL;
        return half ^ signMask;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int ScoreHeisenbergResidualOnCore(sbyte delta)
        => CuboctahedronAutomatonKernel.ScoreHeisenbergStasis9(new Ring9((ushort)(512 + delta))).Bits9;

    public static int ProbeBestPredictorMode(ReadOnlySpan<byte> src)
    {
        int len = src.Length;
        Span<int> scores = stackalloc int[5];
        scores.Clear();

        int canProbe = ~((len - 128) >> 31);
        int step = Math.Max(1, (len - 36) / 64);
        int limit = (len & canProbe) | (36 & ~canProbe);

        for (int i = 36; i < limit; i += step)
        {
            sbyte v0 = (sbyte)src[i];
            sbyte v1 = (sbyte)(byte)(src[i] - src[i - 1]);
            sbyte v4 = (sbyte)(byte)(src[i] - src[i - 4]);
            sbyte v18 = (sbyte)(byte)(src[i] - src[i - 18]);
            sbyte v18Verlet = (sbyte)(byte)(src[i] - 2 * src[i - 18] + src[i - 36]);

            scores[0] += ScoreHeisenbergResidualOnCore(v0);
            scores[1] += ScoreHeisenbergResidualOnCore(v1);
            scores[2] += ScoreHeisenbergResidualOnCore(v4);
            scores[3] += ScoreHeisenbergResidualOnCore(v18);
            scores[4] += ScoreHeisenbergResidualOnCore(v18Verlet);
        }

        int bestMode = 0;
        int maxS = scores[0];
        for (int m = 1; m <= 4; m++)
        {
            int take = ~((scores[m] - maxS - 9) >> 31);
            maxS = (maxS & ~take) | (scores[m] & take);
            bestMode = (bestMode & ~take) | (m & take);
        }
        return bestMode & canProbe;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void ApplyForwardPredictor(ReadOnlySpan<byte> src, Span<byte> dst, int mode)
    {
        int safeMode = Math.Clamp(mode, 0, 4);
        CuboctahedronAutomatonKernel.ComputeVerletResidualBuffer(
            src, dst, ModeStrides[safeMode], ModeMask1[safeMode], ModeMask2[safeMode]);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void ApplyInversePredictor(ReadOnlySpan<byte> src, Span<byte> dst, int mode)
    {
        int safeMode = Math.Clamp(mode, 0, 4);
        CuboctahedronAutomatonKernel.ReconstructVerletBuffer(
            src, dst, ModeStrides[safeMode], ModeMask1[safeMode], ModeMask2[safeMode]);
    }

    // -------------------------------------------------------------------------
    // Zero-Allocation 8-Lane (64-Bit SWAR) Quantum Bit-Width Encoder
    // -------------------------------------------------------------------------

    private delegate void LaneEmitAction(
        int code, ulong raw8, ulong u8, ushort matchDist, Span<byte> outBuf, ref int outPos);

    private static readonly LaneEmitAction[] LaneEmitters =
    [
        EmitLaneCode0Zero,
        EmitLaneCode1to7BitPacked,
        EmitLaneCode1to7BitPacked,
        EmitLaneCode1to7BitPacked,
        EmitLaneCode1to7BitPacked,
        EmitLaneCode1to7BitPacked,
        EmitLaneCode1to7BitPacked,
        EmitLaneCode1to7BitPacked,
        EmitLaneCode8Raw64,
        EmitLaneCode9RingMatch64
    ];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void EmitLaneCode0Zero(
        int code, ulong raw8, ulong u8, ushort matchDist, Span<byte> outBuf, ref int outPos)
    {
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void EmitLaneCode1to7BitPacked(
        int b, ulong raw8, ulong u8, ushort matchDist, Span<byte> outBuf, ref int outPos)
    {
        ulong packed = (u8 & 0xFFUL)
                     | (((u8 >>  8) & 0xFFUL) << b)
                     | (((u8 >> 16) & 0xFFUL) << (2 * b))
                     | (((u8 >> 24) & 0xFFUL) << (3 * b))
                     | (((u8 >> 32) & 0xFFUL) << (4 * b))
                     | (((u8 >> 40) & 0xFFUL) << (5 * b))
                     | (((u8 >> 48) & 0xFFUL) << (6 * b))
                     | (((u8 >> 56) & 0xFFUL) << (7 * b));
        BinaryPrimitives.WriteUInt64LittleEndian(outBuf.Slice(outPos, 8), packed);
        outPos += b;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void EmitLaneCode8Raw64(
        int code, ulong raw8, ulong u8, ushort matchDist, Span<byte> outBuf, ref int outPos)
    {
        BinaryPrimitives.WriteUInt64LittleEndian(outBuf.Slice(outPos, 8), raw8);
        outPos += 8;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void EmitLaneCode9RingMatch64(
        int code, ulong raw8, ulong u8, ushort matchDist, Span<byte> outBuf, ref int outPos)
    {
        BinaryPrimitives.WriteUInt16LittleEndian(outBuf.Slice(outPos, 2), matchDist);
        outPos += 2;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int Classify8ByteLane(
        ulong raw8, int pos, Span<int> hashPos, Span<ulong> hashWords, out ulong u8, out ushort matchDist)
    {
        int isNonZero = (int)((raw8 | (0UL - raw8)) >> 63);

        u8 = ZigZagEncode8LanesSwar(raw8);
        ulong or4 = u8 | (u8 >> 32);
        ulong or2 = or4 | (or4 >> 16);
        uint orAll = (uint)((or2 | (or2 >> 8)) & 0xFFUL);
        int b = (BitOperations.Log2(orAll | 1u) + 1) * isNonZero;

        int h = (int)((raw8 * HashMultiplier64) >> (64 - HashBits));
        int refPos = hashPos[h];
        ulong refWord = hashWords[h];
        hashPos[h] = pos;
        hashWords[h] = raw8;

        int dist = pos - refPos;
        int validDist = ~((refPos >> 31) | ((7 - dist) >> 31) | ((dist - 65536) >> 31));
        ulong diffWord = raw8 ^ refWord;
        int wordMatch = 1 - (int)((diffWord | (0UL - diffWord)) >> 63);
        int bGe3 = ~((b - 3) >> 31);
        int useMatch = validDist & (-wordMatch) & bGe3;

        matchDist = (ushort)(dist & useMatch);
        return (b & ~useMatch) | (9 & useMatch);
    }

    private delegate void PairEncodeAction(
        ReadOnlySpan<byte> tSpan, int alignedLen, ref int i,
        Span<int> hashPos, Span<ulong> hashWords, Span<byte> outBuf, ref int outPos);

    private static readonly PairEncodeAction[] PairEncoders =
    [
        EncodeAllZero16ByteRun,
        EncodeNonZero16BytePair
    ];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void EncodeAllZero16ByteRun(
        ReadOnlySpan<byte> tSpan, int alignedLen, ref int i,
        Span<int> hashPos, Span<ulong> hashWords, Span<byte> outBuf, ref int outPos)
    {
        int remAfter16 = alignedLen - (i + 16);
        int nonZeroRel = tSpan.Slice(i + 16, remAfter16).IndexOfAnyExcept((byte)0);
        int extraZeroBytes = nonZeroRel + ((remAfter16 - nonZeroRel) & (nonZeroRel >> 31));
        int totalPairs = 1 + (extraZeroBytes >> 4);
        i += totalPairs << 4;

        while (totalPairs > 0)
        {
            int chunk = Math.Min(256, totalPairs);
            outBuf[outPos++] = 0x00;
            outBuf[outPos++] = (byte)(chunk - 1);
            totalPairs -= chunk;
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void EncodeNonZero16BytePair(
        ReadOnlySpan<byte> tSpan, int alignedLen, ref int i,
        Span<int> hashPos, Span<ulong> hashWords, Span<byte> outBuf, ref int outPos)
    {
        ulong raw0 = BinaryPrimitives.ReadUInt64LittleEndian(tSpan.Slice(i, 8));
        ulong raw1 = BinaryPrimitives.ReadUInt64LittleEndian(tSpan.Slice(i + 8, 8));

        int code0 = Classify8ByteLane(raw0, i, hashPos, hashWords, out ulong u0, out ushort dist0);
        int code1 = Classify8ByteLane(raw1, i + 8, hashPos, hashWords, out ulong u1, out ushort dist1);

        outBuf[outPos++] = (byte)(code0 | (code1 << 4));
        LaneEmitters[code0](code0, raw0, u0, dist0, outBuf, ref outPos);
        LaneEmitters[code1](code1, raw1, u1, dist1, outBuf, ref outPos);
        i += 16;
    }

    private delegate void TailEncodeAction(
        ReadOnlySpan<byte> tSpan, int i, int tailLen,
        Span<int> hashPos, Span<ulong> hashWords, Span<byte> outBuf, ref int outPos);

    private static readonly TailEncodeAction[] TailEncoders = [EncodeTailNoOp, EncodeTail1To15Bytes];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void EncodeTailNoOp(
        ReadOnlySpan<byte> tSpan, int i, int tailLen,
        Span<int> hashPos, Span<ulong> hashWords, Span<byte> outBuf, ref int outPos)
    {
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void EncodeTail1To15Bytes(
        ReadOnlySpan<byte> tSpan, int i, int tailLen,
        Span<int> hashPos, Span<ulong> hashWords, Span<byte> outBuf, ref int outPos)
    {
        Span<byte> pad16 = stackalloc byte[16];
        pad16.Clear();
        tSpan.Slice(i, tailLen).CopyTo(pad16);

        ulong raw0 = BinaryPrimitives.ReadUInt64LittleEndian(pad16.Slice(0, 8));
        ulong raw1 = BinaryPrimitives.ReadUInt64LittleEndian(pad16.Slice(8, 8));

        int isNz0 = (int)((raw0 | (0UL - raw0)) >> 63);
        ulong u0 = ZigZagEncode8LanesSwar(raw0);
        ulong or0 = u0 | (u0 >> 32); or0 |= (or0 >> 16);
        int code0 = (BitOperations.Log2((uint)((or0 | (or0 >> 8)) & 0xFFUL) | 1u) + 1) * isNz0;

        int isNz1 = (int)((raw1 | (0UL - raw1)) >> 63);
        ulong u1 = ZigZagEncode8LanesSwar(raw1);
        ulong or1 = u1 | (u1 >> 32); or1 |= (or1 >> 16);
        int code1 = (BitOperations.Log2((uint)((or1 | (or1 >> 8)) & 0xFFUL) | 1u) + 1) * isNz1;

        outBuf[outPos++] = (byte)(code0 | (code1 << 4));
        LaneEmitters[code0](code0, raw0, u0, 0, outBuf, ref outPos);
        LaneEmitters[code1](code1, raw1, u1, 0, outBuf, ref outPos);
    }

    private delegate int SpanFinalizerAction(ReadOnlySpan<byte> rawFallbackSource, Span<byte> destination, int outPos);
    private static readonly SpanFinalizerAction[] SpanFinalizers = [FinalizeSpanCompressed, FinalizeSpanRawPassthrough];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int FinalizeSpanCompressed(ReadOnlySpan<byte> rawFallbackSource, Span<byte> destination, int outPos)
        => outPos;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int FinalizeSpanRawPassthrough(ReadOnlySpan<byte> rawFallbackSource, Span<byte> destination, int outPos)
    {
        destination[0] = 0xFF;
        rawFallbackSource.CopyTo(destination.Slice(1));
        return rawFallbackSource.Length + 1;
    }

    private static int EncodeResidualsToSpanCore(
        ReadOnlySpan<byte> tSpan, byte modeHeader, ReadOnlySpan<byte> rawFallbackSource, Span<byte> destination)
    {
        int len = tSpan.Length;
        destination[0] = modeHeader;
        int outPos = 1;

        Span<int> hashPos = stackalloc int[HashSize];
        Span<ulong> hashWords = stackalloc ulong[HashSize];
        hashPos.Fill(-1);

        int alignedLen = len & ~15;
        int i = 0;

        while (i < alignedLen)
        {
            ulong raw0 = BinaryPrimitives.ReadUInt64LittleEndian(tSpan.Slice(i, 8));
            ulong raw1 = BinaryPrimitives.ReadUInt64LittleEndian(tSpan.Slice(i + 8, 8));
            ulong anyBits = raw0 | raw1;
            int isNonZeroPair = (int)((anyBits | (0UL - anyBits)) >> 63);
            PairEncoders[isNonZeroPair](tSpan, alignedLen, ref i, hashPos, hashWords, destination, ref outPos);
        }

        int tailLen = len - alignedLen;
        int hasTail = ((0 - tailLen) >> 31) & 1;
        TailEncoders[hasTail](tSpan, alignedLen, tailLen, hashPos, hashWords, destination, ref outPos);

        int isPassthrough = ((len - outPos) >> 31) & 1;
        return SpanFinalizers[isPassthrough](rawFallbackSource, destination, outPos);
    }

    private delegate int CompressSpanEntryAction(ReadOnlySpan<byte> src, Span<byte> destination, Span<byte> scratch);
    private static readonly CompressSpanEntryAction[] CompressSpanEntries = [CompressEmptySpan, CompressNonEmptySpan];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int CompressEmptySpan(ReadOnlySpan<byte> src, Span<byte> destination, Span<byte> scratch)
    {
        destination[0] = 0;
        return 1;
    }

    private static int CompressNonEmptySpan(ReadOnlySpan<byte> src, Span<byte> destination, Span<byte> scratch)
    {
        int len = src.Length;
        int mode = ProbeBestPredictorMode(src);
        Span<byte> trans = scratch.Slice(0, len);
        ApplyForwardPredictor(src, trans, mode);
        return EncodeResidualsToSpanCore(trans, (byte)mode, src, destination);
    }

    /// <summary>
    /// Zero-allocation span encoder (`0` bytes allocated on the managed heap).
    /// Writes the compressed block into <paramref name="destination"/> (capacity &gt;= <see cref="GetMaxCompressedLength"/>)
    /// using <paramref name="scratch"/> (capacity &gt;= <c>src.Length</c>) and returns the number of bytes written.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int CompressSpan(ReadOnlySpan<byte> src, Span<byte> destination, Span<byte> scratch)
    {
        int nonEmpty = ((0 - src.Length) >> 31) & 1;
        return CompressSpanEntries[nonEmpty](src, destination, scratch);
    }

    /// <summary>
    /// Zero-allocation span encoder for pre-decorrelated residual buffers (`mode = 0` direct, no scratch buffer needed).
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int CompressPredecorrelatedSpan(ReadOnlySpan<byte> residuals, Span<byte> destination)
    {
        return EncodeResidualsToSpanCore(residuals, 0, residuals, destination);
    }

    /// <summary>
    /// Convenience block wrapper backed by <see cref="ArrayPool{T}.Shared"/> for intermediate buffers.
    /// </summary>
    public static byte[] CompressBlock(ReadOnlySpan<byte> src)
    {
        int len = src.Length;
        byte[] rentedScratch = ArrayPool<byte>.Shared.Rent(Math.Max(16, len));
        byte[] rentedOut = ArrayPool<byte>.Shared.Rent(GetMaxCompressedLength(len));
        try
        {
            int written = CompressSpan(src, rentedOut, rentedScratch);
            return rentedOut.AsSpan(0, written).ToArray();
        }
        finally
        {
            ArrayPool<byte>.Shared.Return(rentedScratch);
            ArrayPool<byte>.Shared.Return(rentedOut);
        }
    }

    /// <summary>
    /// Convenience block wrapper for pre-decorrelated residuals backed by <see cref="ArrayPool{T}.Shared"/>.
    /// </summary>
    public static byte[] CompressPredecorrelatedResiduals(ReadOnlySpan<byte> residuals)
    {
        int len = residuals.Length;
        byte[] rentedOut = ArrayPool<byte>.Shared.Rent(GetMaxCompressedLength(len));
        try
        {
            int written = CompressPredecorrelatedSpan(residuals, rentedOut);
            return rentedOut.AsSpan(0, written).ToArray();
        }
        finally
        {
            ArrayPool<byte>.Shared.Return(rentedOut);
        }
    }

    // -------------------------------------------------------------------------
    // Zero-Allocation In-Place 8-Lane (64-Bit SWAR) Table-Dispatched Decoder
    // -------------------------------------------------------------------------

    private delegate ulong SafeUInt64ReaderAction(ReadOnlySpan<byte> compressed, int inPos);
    private static readonly SafeUInt64ReaderAction[] SafeUInt64Readers =
    [
        ReadUInt64FastAtLeast8,
        ReadUInt64StackPaddedUnder8
    ];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static ulong ReadUInt64FastAtLeast8(ReadOnlySpan<byte> compressed, int inPos)
    {
        int safePos = Math.Min(inPos, compressed.Length - 8);
        int shiftBits = (inPos - safePos) << 3;
        return BinaryPrimitives.ReadUInt64LittleEndian(compressed.Slice(safePos, 8)) >> shiftBits;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static ulong ReadUInt64StackPaddedUnder8(ReadOnlySpan<byte> compressed, int inPos)
    {
        Span<byte> pad16 = stackalloc byte[16];
        pad16.Clear();
        compressed.CopyTo(pad16);
        return BinaryPrimitives.ReadUInt64LittleEndian(pad16.Slice(inPos, 8));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static ulong SafeReadUInt64LittleEndian(ReadOnlySpan<byte> compressed, int inPos)
    {
        int isUnder8 = ((compressed.Length - 8) >> 31) & 1;
        return SafeUInt64Readers[isUnder8](compressed, inPos);
    }

    private delegate void LaneDecodeAction(
        int code, ReadOnlySpan<byte> compressed, ref int inPos, Span<byte> trans, ref int outPos);

    private static readonly LaneDecodeAction[] LaneDecoders =
    [
        DecodeLaneCode0Zero,
        DecodeLaneCode1to7BitPacked,
        DecodeLaneCode1to7BitPacked,
        DecodeLaneCode1to7BitPacked,
        DecodeLaneCode1to7BitPacked,
        DecodeLaneCode1to7BitPacked,
        DecodeLaneCode1to7BitPacked,
        DecodeLaneCode1to7BitPacked,
        DecodeLaneCode8Raw64,
        DecodeLaneCode9RingMatch64
    ];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void DecodeLaneCode0Zero(
        int code, ReadOnlySpan<byte> compressed, ref int inPos, Span<byte> trans, ref int outPos)
    {
        BinaryPrimitives.WriteUInt64LittleEndian(trans.Slice(outPos, 8), 0UL);
        outPos += 8;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void DecodeLaneCode1to7BitPacked(
        int b, ReadOnlySpan<byte> compressed, ref int inPos, Span<byte> trans, ref int outPos)
    {
        ulong packed = SafeReadUInt64LittleEndian(compressed, inPos);
        inPos += b;
        ulong mask = 0xFFUL >> (8 - b);
        ulong u8 = (packed & mask)
                 | (((packed >> b) & mask) << 8)
                 | (((packed >> (2 * b)) & mask) << 16)
                 | (((packed >> (3 * b)) & mask) << 24)
                 | (((packed >> (4 * b)) & mask) << 32)
                 | (((packed >> (5 * b)) & mask) << 40)
                 | (((packed >> (6 * b)) & mask) << 48)
                 | (((packed >> (7 * b)) & mask) << 56);
        ulong decoded8 = ZigZagDecode8LanesSwar(u8);
        BinaryPrimitives.WriteUInt64LittleEndian(trans.Slice(outPos, 8), decoded8);
        outPos += 8;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void DecodeLaneCode8Raw64(
        int code, ReadOnlySpan<byte> compressed, ref int inPos, Span<byte> trans, ref int outPos)
    {
        ulong raw8 = SafeReadUInt64LittleEndian(compressed, inPos);
        inPos += 8;
        BinaryPrimitives.WriteUInt64LittleEndian(trans.Slice(outPos, 8), raw8);
        outPos += 8;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void DecodeLaneCode9RingMatch64(
        int code, ReadOnlySpan<byte> compressed, ref int inPos, Span<byte> trans, ref int outPos)
    {
        int dist = BinaryPrimitives.ReadUInt16LittleEndian(compressed.Slice(inPos, 2));
        inPos += 2;
        ulong match8 = BinaryPrimitives.ReadUInt64LittleEndian(trans.Slice(outPos - dist, 8));
        BinaryPrimitives.WriteUInt64LittleEndian(trans.Slice(outPos, 8), match8);
        outPos += 8;
    }

    private delegate void PairDecodeAction(
        byte ctrl, ReadOnlySpan<byte> compressed, ref int inPos, Span<byte> trans, ref int outPos);

    private static readonly PairDecodeAction[] PairDecoders =
    [
        DecodeAllZero16ByteRun,
        DecodeNonZero16BytePair
    ];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void DecodeAllZero16ByteRun(
        byte ctrl, ReadOnlySpan<byte> compressed, ref int inPos, Span<byte> trans, ref int outPos)
    {
        int extraPairs = compressed[inPos++];
        int zeroBytes = (1 + extraPairs) << 4;
        int avail = Math.Min(zeroBytes, trans.Length - outPos);
        trans.Slice(outPos, avail).Clear();
        outPos += avail;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void DecodeNonZero16BytePair(
        byte ctrl, ReadOnlySpan<byte> compressed, ref int inPos, Span<byte> trans, ref int outPos)
    {
        int code0 = ctrl & 0x0F;
        int code1 = ctrl >> 4;
        LaneDecoders[code0](code0, compressed, ref inPos, trans, ref outPos);
        LaneDecoders[code1](code1, compressed, ref inPos, trans, ref outPos);
    }

    private delegate void TailDecodeAction(
        ReadOnlySpan<byte> compressed, ref int inPos, Span<byte> destination, int alignedOutLen, int tailLen);

    private static readonly TailDecodeAction[] TailDecoders = [DecodeTailNoOp, DecodeTail1To15Bytes];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void DecodeTailNoOp(
        ReadOnlySpan<byte> compressed, ref int inPos, Span<byte> destination, int alignedOutLen, int tailLen)
    {
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void DecodeTail1To15Bytes(
        ReadOnlySpan<byte> compressed, ref int inPos, Span<byte> destination, int alignedOutLen, int tailLen)
    {
        Span<byte> pad16 = stackalloc byte[16];
        pad16.Clear();
        int localOut = 0;
        byte ctrl = compressed[inPos++];
        int isNonZero = ((0 - ctrl) >> 31) & 1;
        PairDecoders[isNonZero](ctrl, compressed, ref inPos, pad16, ref localOut);
        pad16.Slice(0, tailLen).CopyTo(destination.Slice(alignedOutLen, tailLen));
    }

    private static void DecodeTokensDirectlyToSpan(ReadOnlySpan<byte> compressedBlock, Span<byte> destination)
    {
        int uncompressedLength = destination.Length;
        int alignedOutLen = uncompressedLength & ~15;
        int inPos = 1;
        int outPos = 0;

        while (outPos < alignedOutLen)
        {
            byte ctrl = compressedBlock[inPos++];
            int isNonZero = ((0 - ctrl) >> 31) & 1;
            PairDecoders[isNonZero](ctrl, compressedBlock, ref inPos, destination, ref outPos);
        }

        int tailLen = uncompressedLength - alignedOutLen;
        int hasTail = ((0 - (uncompressedLength - outPos)) >> 31) & 1;
        TailDecoders[hasTail](compressedBlock, ref inPos, destination, alignedOutLen, tailLen);
    }

    private delegate int DecompressSpanEntryAction(ReadOnlySpan<byte> compressedBlock, Span<byte> destination);
    private static readonly DecompressSpanEntryAction[] DecompressSpanEntries =
    [
        DecompressEmptySpan,
        DecompressRawPassthroughSpan,
        DecompressTokensMode0Span,
        DecompressTokensWithInPlaceVerletSpan
    ];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int DecompressEmptySpan(ReadOnlySpan<byte> compressedBlock, Span<byte> destination) => 0;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int DecompressRawPassthroughSpan(ReadOnlySpan<byte> compressedBlock, Span<byte> destination)
    {
        compressedBlock.Slice(1, destination.Length).CopyTo(destination);
        return destination.Length;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int DecompressTokensMode0Span(ReadOnlySpan<byte> compressedBlock, Span<byte> destination)
    {
        DecodeTokensDirectlyToSpan(compressedBlock, destination);
        return destination.Length;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int DecompressTokensWithInPlaceVerletSpan(ReadOnlySpan<byte> compressedBlock, Span<byte> destination)
    {
        DecodeTokensDirectlyToSpan(compressedBlock, destination);
        // ReconstructVerletBuffer operates 100% in-place on `destination` with zero scratch buffer!
        ApplyInversePredictor(destination, destination, compressedBlock[0]);
        return destination.Length;
    }

    /// <summary>
    /// Zero-allocation in-place span decoder (`0` bytes allocated on the managed heap, no scratch buffer required).
    /// Decodes directly into <paramref name="destination"/> (whose length must equal the uncompressed block length).
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int DecompressSpan(ReadOnlySpan<byte> compressedBlock, Span<byte> destination)
    {
        int uncompressedLength = destination.Length;
        int nonEmpty = ((0 - uncompressedLength) >> 31) & 1;
        byte mode = compressedBlock[0];
        int isRaw = 1 / ((mode ^ 0xFF) + 1);
        int hasPredictor = ((0 - mode) >> 31) & 1 & (1 - isRaw);
        int entryIdx = nonEmpty * (1 + (1 - isRaw) * (1 + hasPredictor));
        return DecompressSpanEntries[entryIdx](compressedBlock, destination);
    }

    public static byte[] DecompressBlock(ReadOnlySpan<byte> compressedBlock, int uncompressedLength)
    {
        byte[] dst = GC.AllocateUninitializedArray<byte>(uncompressedLength);
        DecompressSpan(compressedBlock, dst);
        return dst;
    }

    // -------------------------------------------------------------------------
    // Multi-Gigabyte (>14 GB) Zero-Allocation Streaming Pipeline
    // -------------------------------------------------------------------------

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int ReadFullChunk(Stream input, Span<byte> buffer)
    {
        int total = 0;
        int lastRead = 1;
        while (total < buffer.Length && lastRead > 0)
        {
            lastRead = input.Read(buffer.Slice(total));
            total += Math.Max(0, lastRead);
        }
        return total;
    }

    /// <summary>
    /// Streams an arbitrary-size input stream (14+ GiB, unlimited 64-bit length) through <see cref="CompressSpan"/>
    /// in fixed-size chunks using a single rented <see cref="ArrayPool{T}.Shared"/> buffer triple (`0` per-chunk GC allocations).
    /// </summary>
    public static (long TotalInputBytes, long TotalCompressedBytes) CompressStream(
        Stream input, Stream output, int chunkSize = 4 * 1024 * 1024)
    {
        byte[] rawBuf = ArrayPool<byte>.Shared.Rent(chunkSize);
        byte[] scratchBuf = ArrayPool<byte>.Shared.Rent(chunkSize);
        byte[] compBuf = ArrayPool<byte>.Shared.Rent(GetMaxCompressedLength(chunkSize) + 16);

        long totalIn = 0;
        long totalOut = 0;
        try
        {
            int readBytes = ReadFullChunk(input, rawBuf.AsSpan(0, chunkSize));
            while (readBytes > 0)
            {
                int compLen = CompressSpan(
                    rawBuf.AsSpan(0, readBytes),
                    compBuf.AsSpan(8),
                    scratchBuf.AsSpan(0, readBytes));
                BinaryPrimitives.WriteInt32LittleEndian(compBuf.AsSpan(0, 4), readBytes);
                BinaryPrimitives.WriteInt32LittleEndian(compBuf.AsSpan(4, 4), compLen);
                output.Write(compBuf.AsSpan(0, 8 + compLen));
                totalIn += readBytes;
                totalOut += 8 + compLen;
                readBytes = ReadFullChunk(input, rawBuf.AsSpan(0, chunkSize));
            }
        }
        finally
        {
            ArrayPool<byte>.Shared.Return(rawBuf);
            ArrayPool<byte>.Shared.Return(scratchBuf);
            ArrayPool<byte>.Shared.Return(compBuf);
        }
        return (totalIn, totalOut);
    }

    /// <summary>
    /// Decompresses an arbitrary-size multi-gigabyte stream produced by <see cref="CompressStream"/>
    /// in-place via <see cref="DecompressSpan"/> using a single rented <see cref="ArrayPool{T}.Shared"/> buffer pair (`0` per-chunk GC allocations).
    /// </summary>
    public static long DecompressStream(
        Stream input, Stream output, int maxChunkSize = 4 * 1024 * 1024)
    {
        byte[] rawBuf = ArrayPool<byte>.Shared.Rent(maxChunkSize);
        byte[] compBuf = ArrayPool<byte>.Shared.Rent(GetMaxCompressedLength(maxChunkSize) + 16);
        Span<byte> hdr8 = stackalloc byte[8];

        long totalOut = 0;
        try
        {
            int hdrRead = ReadFullChunk(input, hdr8);
            while (hdrRead == 8)
            {
                int origLen = BinaryPrimitives.ReadInt32LittleEndian(hdr8.Slice(0, 4));
                int compLen = BinaryPrimitives.ReadInt32LittleEndian(hdr8.Slice(4, 4));
                _ = ReadFullChunk(input, compBuf.AsSpan(0, compLen));
                DecompressSpan(compBuf.AsSpan(0, compLen), rawBuf.AsSpan(0, origLen));
                output.Write(rawBuf.AsSpan(0, origLen));
                totalOut += origLen;
                hdrRead = ReadFullChunk(input, hdr8);
            }
        }
        finally
        {
            ArrayPool<byte>.Shared.Return(rawBuf);
            ArrayPool<byte>.Shared.Return(compBuf);
        }
        return totalOut;
    }
}
