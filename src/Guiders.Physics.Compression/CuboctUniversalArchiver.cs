using System;
using System.IO;
using System.IO.Compression;
using Guiders.Physics.Core;

namespace Guiders.Physics.Compression;

/// <summary>
/// Universal raw-byte block and stream archiver (like GZip) with automatic Heisenberg [-40, +40]
/// deadband probing across reversible 1st-order and 2-point Verlet predictors.
/// Works on ANY arbitrary byte array or file stream without format assumptions.
/// </summary>
public static class CuboctUniversalArchiver
{
    public const int DefaultStride18 = 3 * 6; // 18 bytes (e.g., 6ch x 24-bit or tabular records)
    public const int HeisenbergLimit = 41 - 1; // 40 (< 41)

    /// <summary>
    /// Selects the lowest-entropy reversible predictor mode (0..4) for an arbitrary byte buffer
    /// by probing how many residuals fall inside the Heisenberg deadband [-40, +40].
    /// </summary>
    public static int ProbeBestPredictorMode(ReadOnlySpan<byte> src)
    {
        int len = src.Length;
        if (len < 128) return 0;

        int c0 = 0, c1 = 0, c4 = 0, c18 = 0, c18v = 0;
        int step = Math.Max(1, len / 4096);
        for (int i = 36; i < len; i += step)
        {
            int v0 = (sbyte)src[i];
            if (v0 >= -HeisenbergLimit && v0 <= HeisenbergLimit) c0++;

            int v1 = (sbyte)(byte)(src[i] - src[i - 1]);
            if (v1 >= -HeisenbergLimit && v1 <= HeisenbergLimit) c1++;

            int v4 = (sbyte)(byte)(src[i] - src[i - 4]);
            if (v4 >= -HeisenbergLimit && v4 <= HeisenbergLimit) c4++;

            int v18 = (sbyte)(byte)(src[i] - src[i - DefaultStride18]);
            if (v18 >= -HeisenbergLimit && v18 <= HeisenbergLimit) c18++;

            int v18Verlet = (sbyte)(byte)(src[i] - 2 * src[i - DefaultStride18] + src[i - 2 * DefaultStride18]);
            if (v18Verlet >= -HeisenbergLimit && v18Verlet <= HeisenbergLimit) c18v++;
        }

        int bestMode = 0;
        int maxC = c0;
        if (c1 > maxC + 64) { maxC = c1; bestMode = 1; }
        if (c4 > maxC + 64) { maxC = c4; bestMode = 2; }
        if (c18 > maxC + 64) { maxC = c18; bestMode = 3; }
        if (c18v > maxC + 64) { bestMode = 4; }
        return bestMode;
    }

    public static byte[] CompressBlock(ReadOnlySpan<byte> src, CompressionLevel level = CompressionLevel.Optimal)
    {
        int len = src.Length;
        int mode = ProbeBestPredictorMode(src);
        byte[] transformed = new byte[len];

        if (mode == 0)
        {
            src.CopyTo(transformed);
        }
        else if (mode == 1)
        {
            transformed[0] = src[0];
            for (int i = 1; i < len; i++) transformed[i] = (byte)(src[i] - src[i - 1]);
        }
        else if (mode == 2)
        {
            for (int i = 0; i < 4; i++) transformed[i] = src[i];
            for (int i = 4; i < len; i++) transformed[i] = (byte)(src[i] - src[i - 4]);
        }
        else if (mode == 3)
        {
            for (int i = 0; i < DefaultStride18; i++) transformed[i] = src[i];
            for (int i = DefaultStride18; i < len; i++) transformed[i] = (byte)(src[i] - src[i - DefaultStride18]);
        }
        else
        {
            for (int i = 0; i < DefaultStride18; i++) transformed[i] = src[i];
            for (int i = DefaultStride18; i < 2 * DefaultStride18; i++) transformed[i] = (byte)(src[i] - src[i - DefaultStride18]);
            for (int i = 2 * DefaultStride18; i < len; i++) transformed[i] = (byte)(src[i] - 2 * src[i - DefaultStride18] + src[i - 2 * DefaultStride18]);
        }

        using var ms = new MemoryStream();
        ms.WriteByte((byte)mode);
        using (var gz = new GZipStream(ms, level, true))
            gz.Write(transformed, 0, len);
        byte[] comp = ms.ToArray();

        if (mode != 0)
        {
            using var msRaw = new MemoryStream();
            msRaw.WriteByte(0);
            using (var gzRaw = new GZipStream(msRaw, level, true))
                gzRaw.Write(src);
            if (msRaw.Length < comp.Length)
                return msRaw.ToArray();
        }
        return comp;
    }

    public static byte[] DecompressBlock(ReadOnlySpan<byte> compressedBlock, int uncompressedLength)
    {
        int mode = compressedBlock[0];
        byte[] decTrans = new byte[uncompressedLength];
        using (var msIn = new MemoryStream(compressedBlock[1..].ToArray()))
        using (var gzIn = new GZipStream(msIn, CompressionMode.Decompress))
        {
            int rTot = 0;
            while (rTot < uncompressedLength)
            {
                int n = gzIn.Read(decTrans, rTot, uncompressedLength - rTot);
                if (n <= 0) break;
                rTot += n;
            }
        }

        if (mode == 0) return decTrans;

        byte[] dst = new byte[uncompressedLength];
        if (mode == 1)
        {
            dst[0] = decTrans[0];
            for (int i = 1; i < uncompressedLength; i++) dst[i] = (byte)(dst[i - 1] + decTrans[i]);
        }
        else if (mode == 2)
        {
            for (int i = 0; i < 4; i++) dst[i] = decTrans[i];
            for (int i = 4; i < uncompressedLength; i++) dst[i] = (byte)(dst[i - 4] + decTrans[i]);
        }
        else if (mode == 3)
        {
            for (int i = 0; i < DefaultStride18; i++) dst[i] = decTrans[i];
            for (int i = DefaultStride18; i < uncompressedLength; i++) dst[i] = (byte)(dst[i - DefaultStride18] + decTrans[i]);
        }
        else
        {
            for (int i = 0; i < DefaultStride18; i++) dst[i] = decTrans[i];
            for (int i = DefaultStride18; i < 2 * DefaultStride18; i++) dst[i] = (byte)(dst[i - DefaultStride18] + decTrans[i]);
            for (int i = 2 * DefaultStride18; i < uncompressedLength; i++)
                dst[i] = (byte)(2 * dst[i - DefaultStride18] - dst[i - 2 * DefaultStride18] + decTrans[i]);
        }
        return dst;
    }
}
