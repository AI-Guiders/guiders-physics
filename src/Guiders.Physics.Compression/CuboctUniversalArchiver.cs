using System;
using System.Buffers.Binary;
using System.Numerics;
using System.Runtime.CompilerServices;
using Guiders.Physics.Core;

namespace Guiders.Physics.Compression;

/// <summary>
/// Pure native branchless (zero conditional branching statements) lossless archiver built directly on <see cref="Guiders.Physics.Core"/>
/// (<see cref="Ring9"/>, <see cref="VerletCell72"/> 8-lane architecture, and <see cref="CuboctahedronAutomatonKernel"/>):
///   - Heisenberg deadband probing via <see cref="CuboctahedronAutomatonKernel.ScoreHeisenbergStasis9"/> (`DPhiDisc = 41`),
///   - Unified branchless SIMD <see cref="CuboctahedronAutomatonKernel.ComputeVerletResidualBuffer"/> and
///     <see cref="CuboctahedronAutomatonKernel.ReconstructVerletBuffer"/> across all 5 predictor modes (`0..4`),
///   - 64-bit SWAR (8-Lane <see cref="VerletCell72"/> register) quantum bit-width packer (`b in 0..8` bits/lane + O(1) 8-byte ring match `code = 9`),
///   - Table-dispatched encoder and decoder (zero `System.IO.Compression` dependencies).
/// </summary>
public static class CuboctUniversalArchiver
{
    public static readonly int DefaultStride18 = Ring9.D.Bits9 * Ring9.NSq.Bits9;             // 3 * 6 = 18
    public static readonly int HeisenbergMod41 = Ring9.DPhiDisc.Bits9;                        // 41
    public static readonly int HeisenbergLimit = Ring9.DPhiDisc.Bits9 - Ring9.MuVac.Bits9;    // 40

    private const int HashBits = 12;
    private const int HashSize = 1 << HashBits; // 4096-entry direct-mapped 8-byte ring hash table
    private const ulong HashMultiplier64 = 11400714819323198485UL;

    private static readonly int[] ModeStrides = [1, 1, 4, 18, 18];
    private static readonly byte[] ModeMask1  = [0x00, 0xFF, 0xFF, 0xFF, 0xFF];
    private static readonly byte[] ModeMask2  = [0x00, 0x00, 0x00, 0x00, 0xFF];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static byte ZigZagEncodeSignedByte(sbyte delta)
        => (byte)((delta << 1) ^ (delta >> 7));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static sbyte ZigZagDecodeSignedByte(byte u)
        => (sbyte)((u >> 1) ^ -(u & 1));

    /// <summary>
    /// 64-bit SWAR simultaneous ZigZag encoding of 8 signed bytes inside a single <see cref="ulong"/> register.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static ulong ZigZagEncode8LanesSwar(ulong raw8)
    {
        ulong shl1 = (raw8 & 0x7F7F7F7F7F7F7F7FUL) << 1;
        ulong signMask = ((raw8 >> 7) & 0x0101010101010101UL) * 0xFFUL;
        return shl1 ^ signMask;
    }

    /// <summary>
    /// 64-bit SWAR simultaneous ZigZag decoding of 8 unsigned bytes inside a single <see cref="ulong"/> register.
    /// </summary>
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

    /// <summary>
    /// Selects the lowest-entropy reversible predictor mode (0..4) in O(1) time (sampling 64 points)
    /// using <see cref="CuboctahedronAutomatonKernel.ScoreHeisenbergStasis9"/> from <see cref="Guiders.Physics.Core"/>.
    /// </summary>
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
    // 8-Lane (64-Bit SWAR) Quantum Bit-Width Encoder & Table-Dispatched Emitters
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
        // 0 payload bytes emitted for an all-zero 8-byte lane!
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
        // 1. Check all-zero 8-byte lane (code = 0)
        int isNonZero = (int)((raw8 | (0UL - raw8)) >> 63); // 1 when raw8 != 0, 0 when raw8 == 0

        // 2. Compute SWAR ZigZag and quantum bit-width b in [1..8]
        u8 = ZigZagEncode8LanesSwar(raw8);
        ulong or4 = u8 | (u8 >> 32);
        ulong or2 = or4 | (or4 >> 16);
        uint orAll = (uint)((or2 | (or2 >> 8)) & 0xFFUL);
        int b = (BitOperations.Log2(orAll | 1u) + 1) * isNonZero; // 0 when raw8 == 0, else 1..8

        // 3. Probe O(1) 8-byte ring match dictionary (only profitable when b >= 3, saving 1..6 bytes)
        int h = (int)((raw8 * HashMultiplier64) >> (64 - HashBits));
        int refPos = hashPos[h];
        ulong refWord = hashWords[h];
        hashPos[h] = pos;
        hashWords[h] = raw8;

        int dist = pos - refPos;
        int validDist = ~((refPos >> 31) | ((7 - dist) >> 31) | ((dist - 65536) >> 31)); // dist in [8..65535]
        ulong diffWord = raw8 ^ refWord;
        int wordMatch = 1 - (int)((diffWord | (0UL - diffWord)) >> 63);
        int bGe3 = ~((b - 3) >> 31);
        int useMatch = validDist & (-wordMatch) & bGe3; // -1 when match should be used, else 0

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
        int extraPairs = Math.Min(255, extraZeroBytes >> 4);

        outBuf[outPos++] = 0x00;
        outBuf[outPos++] = (byte)extraPairs;
        i += (1 + extraPairs) << 4;
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

        // Do not use ring match on padded tail to avoid out-of-bounds reconstruction
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

    private static byte[] EncodeResidualsCore(ReadOnlySpan<byte> tSpan, byte modeHeader, ReadOnlySpan<byte> rawFallbackSource)
    {
        int len = tSpan.Length;
        byte[] outBuf = GC.AllocateUninitializedArray<byte>(len + (len >> 3) + 64);
        outBuf[0] = modeHeader;
        int outPos = 1;

        Span<int> hashPos = stackalloc int[HashSize];
        Span<ulong> hashWords = stackalloc ulong[HashSize];
        hashPos.Fill(-1);

        int alignedLen = len & ~15;
        int i = 0;
        Span<byte> outSpan = outBuf;

        while (i < alignedLen)
        {
            ulong raw0 = BinaryPrimitives.ReadUInt64LittleEndian(tSpan.Slice(i, 8));
            ulong raw1 = BinaryPrimitives.ReadUInt64LittleEndian(tSpan.Slice(i + 8, 8));
            ulong anyBits = raw0 | raw1;
            int isNonZeroPair = (int)((anyBits | (0UL - anyBits)) >> 63); // 0 -> AllZeroRun, 1 -> NonZeroPair
            PairEncoders[isNonZeroPair](tSpan, alignedLen, ref i, hashPos, hashWords, outSpan, ref outPos);
        }

        int tailLen = len - alignedLen;
        int hasTail = ((0 - tailLen) >> 31) & 1;
        TailEncoders[hasTail](tSpan, alignedLen, tailLen, hashPos, hashWords, outSpan, ref outPos);

        int isPassthrough = ((len - outPos) >> 31) & 1;
        return CompressFinalizers[isPassthrough](rawFallbackSource, outBuf, outPos);
    }

    private delegate byte[] CompressEntryAction(ReadOnlySpan<byte> src);
    private static readonly CompressEntryAction[] CompressEntries = [CompressEmptyBlock, CompressNonEmptyBlock];

    private delegate byte[] CompressFinalizerAction(ReadOnlySpan<byte> src, byte[] outBuf, int outPos);
    private static readonly CompressFinalizerAction[] CompressFinalizers = [FinalizeCompressedSlice, FinalizeRawPassthrough];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static byte[] CompressEmptyBlock(ReadOnlySpan<byte> src) => [0];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static byte[] FinalizeCompressedSlice(ReadOnlySpan<byte> src, byte[] outBuf, int outPos)
        => outBuf.AsSpan(0, outPos).ToArray();

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static byte[] FinalizeRawPassthrough(ReadOnlySpan<byte> src, byte[] outBuf, int outPos)
    {
        byte[] rawOut = GC.AllocateUninitializedArray<byte>(src.Length + 1);
        rawOut[0] = 0xFF;
        src.CopyTo(rawOut.AsSpan(1));
        return rawOut;
    }

    private static byte[] CompressNonEmptyBlock(ReadOnlySpan<byte> src)
    {
        int len = src.Length;
        int mode = ProbeBestPredictorMode(src);
        byte[] trans = GC.AllocateUninitializedArray<byte>(len);
        ApplyForwardPredictor(src, trans, mode);
        return EncodeResidualsCore(trans, (byte)mode, src);
    }

    /// <summary>
    /// Single-pass table-dispatched native entropy encoder (zero conditional branching statements).
    /// </summary>
    public static byte[] CompressBlock(ReadOnlySpan<byte> src)
    {
        int nonEmpty = ((0 - src.Length) >> 31) & 1;
        return CompressEntries[nonEmpty](src);
    }

    /// <summary>
    /// Single-pass table-dispatched native entropy encoder for pre-decorrelated residual buffers (`mode = 0` direct).
    /// </summary>
    public static byte[] CompressPredecorrelatedResiduals(ReadOnlySpan<byte> residuals)
    {
        return EncodeResidualsCore(residuals, 0, residuals);
    }

    // -------------------------------------------------------------------------
    // 8-Lane (64-Bit SWAR) Table-Dispatched Decoder
    // -------------------------------------------------------------------------

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
        ulong packed = BinaryPrimitives.ReadUInt64LittleEndian(compressed.Slice(inPos, 8));
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
        ulong raw8 = BinaryPrimitives.ReadUInt64LittleEndian(compressed.Slice(inPos, 8));
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
        trans.Slice(outPos, zeroBytes).Clear();
        outPos += zeroBytes;
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

    private delegate byte[] DecompressEntryAction(ReadOnlySpan<byte> compressedBlock, int uncompressedLength);
    private static readonly DecompressEntryAction[] DecompressEntries =
    [
        DecompressEmptyBlock,
        DecompressRawPassthroughBlock,
        DecompressTokensMode0,
        DecompressTokensWithInverseVerlet
    ];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static byte[] DecompressEmptyBlock(ReadOnlySpan<byte> compressedBlock, int uncompressedLength) => [];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static byte[] DecompressRawPassthroughBlock(ReadOnlySpan<byte> compressedBlock, int uncompressedLength)
        => compressedBlock.Slice(1, uncompressedLength).ToArray();

    private static byte[] DecodeTokenStreamCore(ReadOnlySpan<byte> compressedBlock, int uncompressedLength)
    {
        // Pad compressed buffer by 8 bytes so 64-bit ReadUInt64LittleEndian is always safe at stream end
        byte[] paddedIn = GC.AllocateUninitializedArray<byte>(compressedBlock.Length + 8);
        compressedBlock.CopyTo(paddedIn);
        ReadOnlySpan<byte> inSpan = paddedIn;

        int alignedOutLen = (uncompressedLength + 15) & ~15;
        byte[] transPadded = GC.AllocateUninitializedArray<byte>(alignedOutLen);
        Span<byte> tSpan = transPadded;

        int inPos = 1;
        int outPos = 0;

        while (outPos < uncompressedLength)
        {
            byte ctrl = inSpan[inPos++];
            int isNonZero = ((0 - ctrl) >> 31) & 1;
            PairDecoders[isNonZero](ctrl, inSpan, ref inPos, tSpan, ref outPos);
        }

        return transPadded.AsSpan(0, uncompressedLength).ToArray();
    }

    private static byte[] DecompressTokensMode0(ReadOnlySpan<byte> compressedBlock, int uncompressedLength)
        => DecodeTokenStreamCore(compressedBlock, uncompressedLength);

    private static byte[] DecompressTokensWithInverseVerlet(ReadOnlySpan<byte> compressedBlock, int uncompressedLength)
    {
        byte[] trans = DecodeTokenStreamCore(compressedBlock, uncompressedLength);
        byte[] dst = GC.AllocateUninitializedArray<byte>(uncompressedLength);
        ApplyInversePredictor(trans, dst, compressedBlock[0]);
        return dst;
    }

    /// <summary>
    /// Table-dispatched branchless native entropy decoder (zero conditional branching statements) paired with
    /// <see cref="CuboctahedronAutomatonKernel.ReconstructVerletBuffer"/>.
    /// </summary>
    public static byte[] DecompressBlock(ReadOnlySpan<byte> compressedBlock, int uncompressedLength)
    {
        int nonEmpty = ((0 - uncompressedLength) >> 31) & 1;
        byte mode = compressedBlock[0];
        int isRaw = 1 / ((mode ^ 0xFF) + 1);
        int hasPredictor = ((0 - mode) >> 31) & 1 & (1 - isRaw);
        int entryIdx = nonEmpty * (1 + (1 - isRaw) * (1 + hasPredictor));
        return DecompressEntries[entryIdx](compressedBlock, uncompressedLength);
    }
}
