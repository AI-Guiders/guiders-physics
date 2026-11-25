using System;
using System.Buffers.Binary;
using System.Runtime.CompilerServices;
using Guiders.Physics.Core;

namespace Guiders.Physics.Compression;

/// <summary>
/// Pure native branchless (zero conditional branching statements) lossless archiver built directly on <see cref="Guiders.Physics.Core"/>
/// (<see cref="Ring9"/> and <see cref="CuboctahedronAutomatonKernel"/>):
///   - Heisenberg deadband probing via <see cref="CuboctahedronAutomatonKernel.ScoreHeisenbergStasis9"/> (`DPhiDisc = 41`),
///   - Unified branchless SIMD <see cref="CuboctahedronAutomatonKernel.ComputeVerletResidualBuffer"/> and
///     <see cref="CuboctahedronAutomatonKernel.ReconstructVerletBuffer"/> across all 5 predictor modes (`0..4`),
///   - Hardware SIMD zero-stasis scanner (<see cref="MemoryExtensions.IndexOfAnyExcept{T}(ReadOnlySpan{T}, T)"/>),
///   - Table-dispatched FSM token encoder and decoder (zero `System.IO.Compression` dependencies).
/// </summary>
public static class CuboctUniversalArchiver
{
    public static readonly int DefaultStride18 = Ring9.D.Bits9 * Ring9.NSq.Bits9;             // 3 * 6 = 18
    public static readonly int HeisenbergMod41 = Ring9.DPhiDisc.Bits9;                        // 41
    public static readonly int HeisenbergLimit = Ring9.DPhiDisc.Bits9 - Ring9.MuVac.Bits9;    // 40

    private const int HashBits = 12;
    private const int HashSize = 1 << HashBits; // 4096-entry direct-mapped ring hash table
    private const uint HashMultiplier = 2654435761u;

    // Unified kernel parameters for modes 0..4:
    //   Mode 0: raw (stride=1, m1=0x00, m2=0x00) -> pred = 0
    //   Mode 1: stride-1 delta (stride=1, m1=0xFF, m2=0x00) -> pred = s[i-1]
    //   Mode 2: stride-4 delta (stride=4, m1=0xFF, m2=0x00) -> pred = s[i-4]
    //   Mode 3: stride-18 delta (stride=18, m1=0xFF, m2=0x00) -> pred = s[i-18]
    //   Mode 4: stride-18 2-Point Verlet parabola (stride=18, m1=0xFF, m2=0xFF) -> pred = 2*s[i-18] - s[i-36]
    private static readonly int[] ModeStrides = [1, 1, 4, 18, 18];
    private static readonly byte[] ModeMask1  = [0x00, 0xFF, 0xFF, 0xFF, 0xFF];
    private static readonly byte[] ModeMask2  = [0x00, 0x00, 0x00, 0x00, 0xFF];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static byte ZigZagEncodeSignedByte(sbyte delta)
        => (byte)((delta << 1) ^ (delta >> 7));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static sbyte ZigZagDecodeSignedByte(byte u)
        => (sbyte)((u >> 1) ^ -(u & 1));

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

        int canProbe = ~((len - 128) >> 31); // -1 when len >= 128, 0 when len < 128
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

    /// <summary>
    /// Delegates forward Verlet prediction directly to <see cref="CuboctahedronAutomatonKernel.ComputeVerletResidualBuffer"/>.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void ApplyForwardPredictor(ReadOnlySpan<byte> src, Span<byte> dst, int mode)
    {
        int safeMode = Math.Clamp(mode, 0, 4);
        CuboctahedronAutomatonKernel.ComputeVerletResidualBuffer(
            src, dst, ModeStrides[safeMode], ModeMask1[safeMode], ModeMask2[safeMode]);
    }

    /// <summary>
    /// Delegates inverse Verlet reconstruction directly to <see cref="CuboctahedronAutomatonKernel.ReconstructVerletBuffer"/>.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void ApplyInversePredictor(ReadOnlySpan<byte> src, Span<byte> dst, int mode)
    {
        int safeMode = Math.Clamp(mode, 0, 4);
        CuboctahedronAutomatonKernel.ReconstructVerletBuffer(
            src, dst, ModeStrides[safeMode], ModeMask1[safeMode], ModeMask2[safeMode]);
    }

    // -------------------------------------------------------------------------
    // Table-Dispatched Branchless VarInt & Token Emitters
    // -------------------------------------------------------------------------

    private delegate int VarIntEmitAction(Span<byte> dst, int pos, int value);
    private static readonly VarIntEmitAction[] VarIntEmitters = [EmitVarIntNoOp, EmitVarIntWrite];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int EmitVarIntNoOp(Span<byte> dst, int pos, int value) => pos;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int EmitVarIntWrite(Span<byte> dst, int pos, int value)
    {
        uint v = (uint)value;
        while (v >= 0x80u)
        {
            dst[pos++] = (byte)((v & 0x7Fu) | 0x80u);
            v >>= 7;
        }
        dst[pos++] = (byte)v;
        return pos;
    }

    private delegate void ZeroTokenEmitAction(Span<byte> outBuf, ref int outPos, int zRun);
    private static readonly ZeroTokenEmitAction[] ZeroTokenEmitters = [EmitShortZeroToken, EmitExtendedZeroToken];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void EmitShortZeroToken(Span<byte> outBuf, ref int outPos, int zRun)
    {
        outBuf[outPos++] = (byte)(40 + zRun);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void EmitExtendedZeroToken(Span<byte> outBuf, ref int outPos, int zRun)
    {
        int ext = zRun - 24;
        int clamped6 = Math.Min(ext, 63);
        outBuf[outPos++] = (byte)(0x40 | clamped6);
        int needVar = ((62 - ext) >> 31) & 1;
        outPos = VarIntEmitters[needVar](outBuf, outPos, ext - 63);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void EmitZeroRunToken(Span<byte> outBuf, ref int outPos, int zRun)
    {
        int isLong = ((23 - zRun) >> 31) & 1;
        ZeroTokenEmitters[isLong](outBuf, ref outPos, zRun);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void EmitMatchToken(Span<byte> outBuf, ref int outPos, int mLen, int dist)
    {
        int lenCode = mLen - 4;
        int clamped6 = Math.Min(lenCode, 63);
        outBuf[outPos++] = (byte)(0xC0 | clamped6);
        int needVar = ((62 - lenCode) >> 31) & 1;
        outPos = VarIntEmitters[needVar](outBuf, outPos, lenCode - 63);
        BinaryPrimitives.WriteUInt16LittleEndian(outBuf.Slice(outPos, 2), (ushort)dist);
        outPos += 2;
    }

    // -------------------------------------------------------------------------
    // Table-Dispatched Literal Flushers (0 = Empty, 1 = 1-Byte Heisenberg, 2 = Chunks)
    // -------------------------------------------------------------------------

    private delegate void LiteralFlushAction(
        ReadOnlySpan<byte> data, Span<byte> outBuf, ref int outPos, int litStart, ref int litCount);

    private static readonly LiteralFlushAction[] LiteralFlushers =
    [
        FlushLiteralsNoOp,
        FlushSingleHeisenberg1Byte,
        FlushLiteralChunks
    ];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void FlushLiteralsNoOp(
        ReadOnlySpan<byte> data, Span<byte> outBuf, ref int outPos, int litStart, ref int litCount)
    {
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void FlushSingleHeisenberg1Byte(
        ReadOnlySpan<byte> data, Span<byte> outBuf, ref int outPos, int litStart, ref int litCount)
    {
        outBuf[outPos++] = data[litStart];
        litCount = 0;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void FlushLiteralChunks(
        ReadOnlySpan<byte> data, Span<byte> outBuf, ref int outPos, int litStart, ref int litCount)
    {
        int idx = litStart;
        int rem = litCount;
        while (rem > 0)
        {
            int chunk = Math.Min(rem, 64);
            outBuf[outPos++] = (byte)(0x80 | (chunk - 1));
            data.Slice(idx, chunk).CopyTo(outBuf.Slice(outPos));
            outPos += chunk;
            idx += chunk;
            rem -= chunk;
        }
        litCount = 0;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void FlushLiterals(
        ReadOnlySpan<byte> data, Span<byte> outBuf, ref int outPos, int litStart, ref int litCount)
    {
        int hasLits = ((0 - litCount) >> 31) & 1; // 1 when litCount >= 1, else 0
        int isSingle = (1 / (Math.Abs(litCount - 1) + 1)) & hasLits;
        int safeIdx = Math.Clamp(litStart, 0, Math.Max(0, data.Length - 1));
        int inHeisenberg40 = ((int)data[safeIdx] - 41) >> 31; // -1 when data[safeIdx] <= 40
        int isHeis1 = isSingle & (inHeisenberg40 & 1);
        int flushIdx = hasLits + (hasLits & (1 - isHeis1)); // 0 -> NoOp, 1 -> SingleHeisenberg, 2 -> Chunks
        LiteralFlushers[flushIdx](data, outBuf, ref outPos, litStart, ref litCount);
    }

    // -------------------------------------------------------------------------
    // Table-Dispatched FSM Encoder Step (0 = Literal, 1 = Ring Match, 2 = Zero Stasis)
    // -------------------------------------------------------------------------

    private delegate void EncoderStepAction(
        ReadOnlySpan<byte> tSpan, Span<int> hashTable, Span<byte> outBuf,
        ref int outPos, ref int i, ref int litStart, ref int litCount);

    private static readonly EncoderStepAction[] EncoderSteps =
    [
        Step0AccumulateLiteral,
        Step1EmitRingMatch,
        Step2EmitZeroStasis
    ];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int ProbeZeroStasisRunLength(ReadOnlySpan<byte> tSpan, int i, int litCount)
    {
        int isZeroByte = ((int)tSpan[i] - 1) >> 31; // -1 when tSpan[i] == 0, else 0
        int nonZeroRel = tSpan.Slice(i).IndexOfAnyExcept((byte)0);
        int rawRun = nonZeroRel + ((tSpan.Length - i - nonZeroRel) & (nonZeroRel >> 31));
        int zRun = rawRun & isZeroByte;
        int worthMask = ~((zRun - 2) >> 31) | ((litCount - 1) >> 31); // zRun >= 2 || litCount == 0
        return zRun & worthMask;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int ProbeRingMatchAvailable(ReadOnlySpan<byte> tSpan, Span<int> hashTable, int i)
    {
        int len = tSpan.Length;
        int canRead4 = ~((len - i - 4) >> 31); // -1 when i + 4 <= len, else 0
        int safePos = i & canRead4;
        uint word = BinaryPrimitives.ReadUInt32LittleEndian(tSpan.Slice(safePos, Math.Min(4, len - safePos)));
        int h = (int)((word * HashMultiplier) >> (32 - HashBits));
        int refPos = hashTable[h];

        int dist = i - refPos;
        int validRef = ~((refPos >> 31) | ((0 - dist) >> 31) | ((dist - 65536) >> 31)) & canRead4;
        int safeRef = Math.Max(0, refPos) & validRef;
        uint refWord = BinaryPrimitives.ReadUInt32LittleEndian(tSpan.Slice(safeRef, Math.Min(4, len - safeRef)));
        int wordsEqual = 1 / ((int)((word ^ refWord) & 0x7FFFFFFFu) + (int)((word ^ refWord) >> 31) + 1);
        return (validRef & 1) & wordsEqual;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Step0AccumulateLiteral(
        ReadOnlySpan<byte> tSpan, Span<int> hashTable, Span<byte> outBuf,
        ref int outPos, ref int i, ref int litStart, ref int litCount)
    {
        int len = tSpan.Length;
        int canRead4 = ~((len - i - 4) >> 31);
        int safePos = i & canRead4;
        uint word = BinaryPrimitives.ReadUInt32LittleEndian(tSpan.Slice(safePos, Math.Min(4, len - safePos)));
        int h = (int)((word * HashMultiplier) >> (32 - HashBits));
        hashTable[h] = (hashTable[h] & ~canRead4) | (i & canRead4);

        int isFirstLit = (litCount - 1) >> 31;
        litStart = (litStart & ~isFirstLit) | (i & isFirstLit);
        litCount++;
        i++;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Step1EmitRingMatch(
        ReadOnlySpan<byte> tSpan, Span<int> hashTable, Span<byte> outBuf,
        ref int outPos, ref int i, ref int litStart, ref int litCount)
    {
        uint word = BinaryPrimitives.ReadUInt32LittleEndian(tSpan.Slice(i, 4));
        int h = (int)((word * HashMultiplier) >> (32 - HashBits));
        int refPos = hashTable[h];
        hashTable[h] = i;
        int dist = i - refPos;

        FlushLiterals(tSpan, outBuf, ref outPos, litStart, ref litCount);
        int mLen = 4 + tSpan.Slice(refPos + 4, tSpan.Length - i - 4).CommonPrefixLength(tSpan.Slice(i + 4));
        EmitMatchToken(outBuf, ref outPos, mLen, dist);
        i += mLen;
        litStart = i;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Step2EmitZeroStasis(
        ReadOnlySpan<byte> tSpan, Span<int> hashTable, Span<byte> outBuf,
        ref int outPos, ref int i, ref int litStart, ref int litCount)
    {
        int nonZeroRel = tSpan.Slice(i).IndexOfAnyExcept((byte)0);
        int zRun = nonZeroRel + ((tSpan.Length - i - nonZeroRel) & (nonZeroRel >> 31));
        FlushLiterals(tSpan, outBuf, ref outPos, litStart, ref litCount);
        EmitZeroRunToken(outBuf, ref outPos, zRun);
        i += zRun;
        litStart = i;
    }

    // -------------------------------------------------------------------------
    // Table-Dispatched Compress Entry & Finalizers
    // -------------------------------------------------------------------------

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

        byte[] outBuf = GC.AllocateUninitializedArray<byte>(len + (len >> 4) + 32);
        outBuf[0] = (byte)mode;
        int outPos = 1;

        Span<int> hashTable = stackalloc int[HashSize];
        hashTable.Fill(-1);

        int i = 0;
        int litStart = 0;
        int litCount = 0;
        ReadOnlySpan<byte> tSpan = trans;

        while (i < len)
        {
            int zRun = ProbeZeroStasisRunLength(tSpan, i, litCount);
            int hasZero = ((0 - zRun) >> 31) & 1;
            int hasMatch = (1 - hasZero) & ProbeRingMatchAvailable(tSpan, hashTable, i);
            int stepIdx = (hasZero << 1) | hasMatch; // 2 -> ZeroStasis, 1 -> RingMatch, 0 -> Literal
            EncoderSteps[stepIdx](tSpan, hashTable, outBuf, ref outPos, ref i, ref litStart, ref litCount);
        }

        FlushLiterals(tSpan, outBuf, ref outPos, litStart, ref litCount);

        int isPassthrough = ((len - outPos) >> 31) & 1; // 1 when outPos >= len + 1, else 0
        return CompressFinalizers[isPassthrough](src, outBuf, outPos);
    }

    /// <summary>
    /// Single-pass table-dispatched native entropy encoder (zero conditional branching statements).
    /// </summary>
    public static byte[] CompressBlock(ReadOnlySpan<byte> src)
    {
        int nonEmpty = ((0 - src.Length) >> 31) & 1;
        return CompressEntries[nonEmpty](src);
    }

    // -------------------------------------------------------------------------
    // Table-Dispatched Decoder (0 = Empty, 1 = RawPassthrough, 2 = Mode0, 3 = Verlet 1..4)
    // -------------------------------------------------------------------------

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int ReadVarIntWhen63(ReadOnlySpan<byte> src, ref int pos, int low6)
    {
        uint result = 0;
        int shift = 0;
        int active = ~((low6 - 63) >> 31);
        while (active != 0)
        {
            byte b = src[pos++];
            result |= (uint)(b & 0x7F) << shift;
            shift += 7;
            active = -((b >> 7) & 1);
        }
        return (int)result;
    }

    private delegate void TokenDecoderAction(
        byte token, ReadOnlySpan<byte> compressedBlock, ref int inPos, Span<byte> trans, ref int outPos);

    private static readonly TokenDecoderAction[] TokenDecoders =
    [
        DecodeClass0HeisenbergOrShortZero,
        DecodeClass1ExtendedZero,
        DecodeClass2LiteralChunk,
        DecodeClass3RingMatch
    ];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void DecodeClass0HeisenbergOrShortZero(
        byte token, ReadOnlySpan<byte> compressedBlock, ref int inPos, Span<byte> trans, ref int outPos)
    {
        int isZeroRun = (40 - token) >> 31; // -1 when token > 40, 0 when token <= 40
        int runLen = (1 & ~isZeroRun) | ((token - 40) & isZeroRun);
        byte fillVal = (byte)(token & ~isZeroRun);
        trans.Slice(outPos, runLen).Fill(fillVal);
        outPos += runLen;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void DecodeClass1ExtendedZero(
        byte token, ReadOnlySpan<byte> compressedBlock, ref int inPos, Span<byte> trans, ref int outPos)
    {
        int low6 = token & 0x3F;
        int zRun = 24 + low6 + ReadVarIntWhen63(compressedBlock, ref inPos, low6);
        trans.Slice(outPos, zRun).Clear();
        outPos += zRun;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void DecodeClass2LiteralChunk(
        byte token, ReadOnlySpan<byte> compressedBlock, ref int inPos, Span<byte> trans, ref int outPos)
    {
        int litLen = (token & 0x3F) + 1;
        compressedBlock.Slice(inPos, litLen).CopyTo(trans.Slice(outPos, litLen));
        inPos += litLen;
        outPos += litLen;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void DecodeClass3RingMatch(
        byte token, ReadOnlySpan<byte> compressedBlock, ref int inPos, Span<byte> trans, ref int outPos)
    {
        int low6 = token & 0x3F;
        int rem = 4 + low6 + ReadVarIntWhen63(compressedBlock, ref inPos, low6);
        int dist = BinaryPrimitives.ReadUInt16LittleEndian(compressedBlock.Slice(inPos, 2));
        inPos += 2;

        while (rem > 0)
        {
            int chunk = Math.Min(dist, rem);
            trans.Slice(outPos - dist, chunk).CopyTo(trans.Slice(outPos, chunk));
            outPos += chunk;
            rem -= chunk;
        }
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
        byte[] trans = GC.AllocateUninitializedArray<byte>(uncompressedLength);
        Span<byte> tSpan = trans;
        int inPos = 1;
        int outPos = 0;
        int inLen = compressedBlock.Length;

        while (inPos < inLen && outPos < uncompressedLength)
        {
            byte token = compressedBlock[inPos++];
            TokenDecoders[token >> 6](token, compressedBlock, ref inPos, tSpan, ref outPos);
        }
        return trans;
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
        int isRaw = 1 / ((mode ^ 0xFF) + 1); // 1 when mode == 0xFF, else 0
        int hasPredictor = ((0 - mode) >> 31) & 1 & (1 - isRaw); // 1 when mode in 1..4
        int entryIdx = nonEmpty * (1 + (1 - isRaw) * (1 + hasPredictor));
        return DecompressEntries[entryIdx](compressedBlock, uncompressedLength);
    }
}
