using System;
using System.Diagnostics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace Guiders.Physics.Core;

/// <summary>
/// Complete integer passport of the M-layer carrier derived purely from spatial dimension D = 3.
/// Even unity is computed as MuVac = D*D - (1 &lt;&lt; D) = 9 - 8 = 1.
/// </summary>
public static class CarrierSpecD3
{
    public const int D = 3;
    public const int MuVac = (D * D) - (1 << D);              // 1
    public const int C0OverCSq = D - MuVac;                   // 2
    public const int NLinks = 2 * D * (D - MuVac);            // 12
    public const int NTri = 1 << D;                           // 8
    public const int NSq = 2 * D;                             // 6
    public const int DInSquare = NLinks / D;                  // 4
    public const int DOutEdge = NSq - MuVac;                  // 5
    public const int DOutVertex = NTri - MuVac;               // 7
    public const int ModBits = MuVac + NTri;                  // 9
    public const int FracBits = ModBits - D;                  // 6
    public const int NRing = 1 << ModBits;                    // 512
    public const int RingMask = NRing - MuVac;                // 511 = 0x1FF
    public const int KP = 1 << FracBits;                      // 64
    public const int QSquare = NRing / DInSquare;             // 128
    public const int PhiPauli = NRing / C0OverCSq;            // 256
    public const int EpsDenom = KP * QSquare;                 // 8192
    public const int NPhi = MuVac + NLinks;                   // 13
    public const int DPhiDisc = (NRing + NPhi + NTri) / NPhi; // 41
    public const int SeamTicks = NPhi * DPhiDisc - NRing;     // 21
    public const int MHat = MuVac + NLinks * NTri;            // 97
    public const int CollapsePhi = QSquare + NLinks - MuVac;  // 139
    public const int IvSeed = PhiPauli + DOutVertex;          // 263 = 0o407
    public const int B1SingleTickBits = DInSquare * ModBits;  // 36 bits
    public const int B2TwoTickBits = NTri * ModBits;          // 72 bits = 9 bytes!

    public static readonly (int Dz, int Dy, int Dx)[] CuboctOffsets12 =
    [
        ( 0,  1,  1), ( 0,  1, -1), ( 0, -1,  1), ( 0, -1, -1),
        ( 1,  0,  1), ( 1,  0, -1), (-1,  0,  1), (-1,  0, -1),
        ( 1,  1,  0), ( 1, -1,  0), (-1,  1,  0), (-1, -1,  0)
    ];
}

/// <summary>
/// Native 9-Bit Modular Ring Element in Z_512 (bits [8:0]) with 100% branchless
/// Gray-Quadrant Butterfly Rotator (Z_4^Gray x Z_128) and Octal Triad (Z_8^3) representation.
/// </summary>
public readonly record struct Ring9
{
    private readonly ushort _bits9;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public Ring9(int value) => _bits9 = (ushort)(value & CarrierSpecD3.RingMask);

    public int UnsignedTicks
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _bits9;
    }

    public int SignedWrap
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => CarrierSpecD3.PhiPauli - ((CarrierSpecD3.PhiPauli - _bits9) & CarrierSpecD3.RingMask);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public (int ProjC, int ProjS) GrayButterflyRotator()
    {
        int q1 = (_bits9 >> 8) & 1;
        int q0 = (_bits9 >> 7) & 1;
        int m7 = _bits9 & 127;
        int comp7 = CarrierSpecD3.QSquare - m7;
        int magC = (1 - q0) * comp7 + q0 * m7;
        int magS = (1 - q0) * m7 + q0 * comp7;
        int sgnC = 1 - 2 * (q1 ^ q0);
        int sgnS = 1 - 2 * q1;
        return (sgnC * magC, sgnS * magS);
    }

    public int ProjC
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => GrayButterflyRotator().ProjC;
    }

    public int ProjS
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => GrayButterflyRotator().ProjS;
    }

    public (int O2, int O1, int O0) OctalTriad => ((_bits9 >> 6) & 7, (_bits9 >> 3) & 7, _bits9 & 7);

    public Ring9 PauliPartner
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new(_bits9 ^ CarrierSpecD3.PhiPauli);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 operator +(Ring9 a, Ring9 b) => new(a._bits9 + b._bits9);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 operator -(Ring9 a, Ring9 b) => new(a._bits9 - b._bits9);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Ring9 operator *(Ring9 a, int k) => new(a._bits9 * k);

    public override string ToString()
    {
        var (o2, o1, o0) = OctalTriad;
        return $"{_bits9} (0o{o2}{o1}{o0}, wrap={SignedWrap:+#;-#;0})";
    }
}

/// <summary>
/// Exact 72-Bit (9-Byte, Pack = 1) Two-Tick Verlet Polar Cell:
/// packs eight 9-bit ring lanes (mu0, phi0_curr, phi0_prev, mu1, phi1_curr, phi1_prev, gate0, gate1)
/// into literally 9 bytes = 72 bits (`B_2(3) = 8 * 9 = 72 bits`) with ZERO 32/64-bit padding!
/// </summary>
[StructLayout(LayoutKind.Sequential, Pack = 1, Size = 9)]
public struct VerletCell72
{
    private ulong _lo64;
    private byte _hi8;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public readonly Ring9 GetLane9(int laneIdx)
    {
        int bitOffset = laneIdx * 9;
        if (bitOffset <= 54)
            return new Ring9((int)((_lo64 >> bitOffset) & 0x1FFUL));
        int loBit = (int)((_lo64 >> 63) & 1UL);
        return new Ring9(loBit | (_hi8 << 1));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void SetLane9(int laneIdx, Ring9 val)
    {
        ulong v = (ulong)(val.UnsignedTicks & 0x1FF);
        int bitOffset = laneIdx * 9;
        if (bitOffset <= 54)
        {
            ulong mask = ~(0x1FFUL << bitOffset);
            _lo64 = (_lo64 & mask) | (v << bitOffset);
        }
        else
        {
            _lo64 = (_lo64 & 0x7FFF_FFFF_FFFF_FFFFUL) | ((v & 1UL) << 63);
            _hi8 = (byte)(v >> 1);
        }
    }

    public Ring9 Mu0      { readonly get => GetLane9(0); set => SetLane9(0, value); }
    public Ring9 Phi0Curr { readonly get => GetLane9(1); set => SetLane9(1, value); }
    public Ring9 Phi0Past { readonly get => GetLane9(2); set => SetLane9(2, value); }
    public Ring9 Mu1      { readonly get => GetLane9(3); set => SetLane9(3, value); }
    public Ring9 Phi1Curr { readonly get => GetLane9(4); set => SetLane9(4, value); }
    public Ring9 Phi1Past { readonly get => GetLane9(5); set => SetLane9(5, value); }
    public Ring9 Gate0    { readonly get => GetLane9(6); set => SetLane9(6, value); }
    public Ring9 Gate1    { readonly get => GetLane9(7); set => SetLane9(7, value); }

    public static VerletCell72 Create(int mu0, int p0Curr, int p0Past, int mu1, int p1Curr, int p1Past)
    {
        VerletCell72 c = default;
        c.Mu0 = new Ring9(mu0);
        c.Phi0Curr = new Ring9(p0Curr);
        c.Phi0Past = new Ring9(p0Past);
        c.Mu1 = new Ring9(mu1);
        c.Phi1Curr = new Ring9(p1Curr);
        c.Phi1Past = new Ring9(p1Past);
        return c;
    }
}

/// <summary>
/// Closed-form Multimodal Sensor Readout from a 72-bit Verlet cell:
/// Optical wavelength [440..718 nm], RGB triad (R+G+B = 128), acoustic pressure P_sound, concentration, and thermal quanta n_E.
/// </summary>
public readonly record struct MultimodalReadout(
    int ConcentrationNorm,
    int WavelengthNm,
    byte ColorR,
    byte ColorG,
    byte ColorB,
    int AcousticPressure,
    int ReactionCurrent,
    int ThermalQuantaNe);

/// <summary>
/// 100% Branchless 12-Neighbor Cuboctahedron O_2(3) Cellular Automaton Core in .NET 11.
/// </summary>
public static class CuboctahedronAutomatonKernel
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int BranchlessSignNonNeg(int x) => 1 - 2 * ((x >> 31) & 1);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int FloorDivInt(long num, long den)
    {
        long q = num / den;
        long r = num % den;
        if ((r != 0) && ((r < 0) ^ (den < 0))) q--;
        return (int)q;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Ring9 GatePhi, int RawWrapped, int DenomD, int EnergyNe) Compute12NeighborGate(
        int muC, Ring9 phiC, ReadOnlySpan<int> nbMu, ReadOnlySpan<Ring9> nbPhi, bool includeEvenDeficit = true)
    {
        long zrRaw = 0;
        long ziRaw = 0;
        long sumMu = 0;
        for (int k = 0; k < 12; k++)
        {
            Ring9 dphi = nbPhi[k] - phiC;
            var (cProj, sProj) = dphi.GrayButterflyRotator();
            int m = nbMu[k];
            zrRaw += (long)m * cProj;
            ziRaw += (long)m * sProj;
            sumMu += m;
        }
        long zr = (long)muC * zrRaw;
        long zi = (long)muC * ziRaw;
        long deficit = (long)muC * sumMu * CarrierSpecD3.QSquare - zr;
        int sgnZi = BranchlessSignNonNeg((int)zi);
        long zDrive = includeEvenDeficit ? (zi + sgnZi * deficit) : zi;
        long rawDen = zr + (long)muC * muC * CarrierSpecD3.QSquare + CarrierSpecD3.EpsDenom;
        long safeDen = rawDen > 0 ? rawDen : CarrierSpecD3.EpsDenom;

        int qFloor = FloorDivInt(CarrierSpecD3.KP * zDrive, safeDen);
        int rawWrap = new Ring9(qFloor).SignedWrap;
        int absRaw = Math.Abs(rawWrap);
        int ne = absRaw / CarrierSpecD3.DPhiDisc;
        int gateSigned = absRaw >= CarrierSpecD3.DPhiDisc ? rawWrap : 0;
        return (new Ring9(gateSigned), rawWrap, (int)safeDen, ne);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static VerletCell72 StepCellForward(VerletCell72 cell, Ring9 gate0, Ring9 gate1)
    {
        Ring9 p0Next = new(2 * cell.Phi0Curr.UnsignedTicks - cell.Phi0Past.UnsignedTicks + gate0.SignedWrap);
        Ring9 p1Next = new(2 * cell.Phi1Curr.UnsignedTicks - cell.Phi1Past.UnsignedTicks + gate1.SignedWrap);
        var next = VerletCell72.Create(
            cell.Mu0.UnsignedTicks, p0Next.UnsignedTicks, cell.Phi0Curr.UnsignedTicks,
            cell.Mu1.UnsignedTicks, p1Next.UnsignedTicks, cell.Phi1Curr.UnsignedTicks);
        next.Gate0 = gate0;
        next.Gate1 = gate1;
        return next;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static VerletCell72 StepCellReverse(VerletCell72 steppedCell, Ring9 gate0AtT, Ring9 gate1AtT)
    {
        Ring9 p0Curr = steppedCell.Phi0Past;
        Ring9 p1Curr = steppedCell.Phi1Past;
        Ring9 p0PastRec = new(2 * p0Curr.UnsignedTicks - steppedCell.Phi0Curr.UnsignedTicks + gate0AtT.SignedWrap);
        Ring9 p1PastRec = new(2 * p1Curr.UnsignedTicks - steppedCell.Phi1Curr.UnsignedTicks + gate1AtT.SignedWrap);
        return VerletCell72.Create(
            steppedCell.Mu0.UnsignedTicks, p0Curr.UnsignedTicks, p0PastRec.UnsignedTicks,
            steppedCell.Mu1.UnsignedTicks, p1Curr.UnsignedTicks, p1PastRec.UnsignedTicks);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static MultimodalReadout ReadCellSensors(VerletCell72 cell)
    {
        int mu0 = cell.Mu0.UnsignedTicks;
        int mu1 = cell.Mu1.UnsignedTicks;
        Ring9 dphi = cell.Phi1Curr - cell.Phi0Curr;
        var (c3, s3) = dphi.GrayButterflyRotator();

        int conc = mu0 * mu0 + mu1 * mu1;
        int wavelengthNm = 579 + FloorDivInt(CarrierSpecD3.CollapsePhi * c3, CarrierSpecD3.QSquare);
        byte r = (byte)Math.Max(0, c3);
        byte g = (byte)Math.Abs(s3);
        byte b = (byte)Math.Max(0, -c3);
        int soundPressure = mu0 * mu1 * (CarrierSpecD3.QSquare - c3);
        int rxnCurrent = mu0 * mu1 * s3;
        int ne = Math.Abs(cell.Gate0.SignedWrap) / CarrierSpecD3.DPhiDisc;

        return new MultimodalReadout(conc, wavelengthNm, r, g, b, soundPressure, rxnCurrent, ne);
    }

    public static int Iv263MuAt(int dz, int dy, int dx)
    {
        int r2 = dz * dz + dy * dy + dx * dx;
        int inBall = r2 <= CarrierSpecD3.NSq * CarrierSpecD3.NSq ? 1 : 0;
        return CarrierSpecD3.MuVac + CarrierSpecD3.NLinks * inBall;
    }

    public static Ring9 Iv263PhiAt(int dz, int dy, int dx, int lane)
    {
        int r2 = dz * dz + dy * dy + dx * dx;
        int inBall = r2 <= CarrierSpecD3.NSq * CarrierSpecD3.NSq ? 1 : 0;
        int oct = (dz >= 0 ? 4 : 0) + (dy >= 0 ? 2 : 0) + (dx >= 0 ? 1 : 0);
        int shiftOct = lane == 0 ? oct : (oct ^ CarrierSpecD3.D);
        int bit = ((CarrierSpecD3.IvSeed % CarrierSpecD3.PhiPauli) >> shiftOct) & 1;
        int rem13 = ((dx + CarrierSpecD3.NLinks) % CarrierSpecD3.NPhi + CarrierSpecD3.NPhi) % CarrierSpecD3.NPhi;
        int bg = (rem13 * CarrierSpecD3.DPhiDisc) % CarrierSpecD3.NRing;
        int phi = (bg + (CarrierSpecD3.D * CarrierSpecD3.DPhiDisc) * inBall * bit) % CarrierSpecD3.NRing;
        return new Ring9(phi);
    }
}

/// <summary>
/// Direct P/Invoke wrapper for NVIDIA CUDA Driver API (`nvcuda.dll`) on RTX 4060 (`sm_89`).
/// </summary>
public static partial class CudaRtx4060Driver
{
    [LibraryImport("nvcuda.dll")]
    private static partial int cuInit(uint flags);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuDeviceGet(out int device, int ordinal);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuDeviceGetName(Span<byte> name, int len, int dev);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuCtxCreate_v2(out IntPtr pctx, uint flags, int dev);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuCtxDestroy_v2(IntPtr ctx);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuModuleLoadDataEx(out IntPtr module, IntPtr image, uint numOptions, IntPtr options, IntPtr optionValues);

    [LibraryImport("nvcuda.dll", StringMarshalling = StringMarshalling.Utf8)]
    private static partial int cuModuleGetFunction(out IntPtr hfunc, IntPtr hmod, string name);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuMemAlloc_v2(out ulong dptr, nuint bytesize);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuMemFree_v2(ulong dptr);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuMemcpyHtoD_v2(ulong dstDevice, IntPtr srcHost, nuint byteCount);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuMemcpyDtoH_v2(IntPtr dstHost, ulong srcDevice, nuint byteCount);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuLaunchKernel(
        IntPtr f, uint gridDimX, uint gridDimY, uint gridDimZ,
        uint blockDimX, uint blockDimY, uint blockDimZ,
        uint sharedMemBytes, IntPtr hStream, IntPtr kernelParams, IntPtr extra);

    [LibraryImport("nvcuda.dll")]
    private static partial int cuCtxSynchronize();

    private const string RotatorVerletPtx = """
        .version 6.5
        .target sm_70
        .address_size 64

        .visible .entry cuboct_rotator_verlet_ptx(
            .param .u64 p_in,
            .param .u64 p_out,
            .param .u32 n_elems
        )
        {
            .reg .pred %p1;
            .reg .b32 %r_tid, %r_ntid, %r_ctaid, %idx, %n, %phi, %w, %tmp, %c_proj;
            .reg .b64 %base_in, %base_out, %addr_in, %addr_out, %off;

            ld.param.u64 %base_in, [p_in];
            ld.param.u64 %base_out, [p_out];
            ld.param.u32 %n, [n_elems];

            mov.u32 %r_tid, %tid.x;
            mov.u32 %r_ntid, %ntid.x;
            mov.u32 %r_ctaid, %ctaid.x;
            mad.lo.u32 %idx, %r_ctaid, %r_ntid, %r_tid;

            setp.ge.u32 %p1, %idx, %n;
            @%p1 bra DONE;

            mul.wide.u32 %off, %idx, 4;
            add.u64 %addr_in, %base_in, %off;
            add.u64 %addr_out, %base_out, %off;

            ld.global.u32 %phi, [%addr_in];
            and.b32 %phi, %phi, 511;

            mov.s32 %tmp, 256;
            sub.s32 %tmp, %tmp, %phi;
            and.b32 %tmp, %tmp, 511;
            mov.s32 %w, 256;
            sub.s32 %w, %w, %tmp;

            abs.s32 %tmp, %w;
            mov.s32 %c_proj, 128;
            sub.s32 %c_proj, %c_proj, %tmp;

            st.global.s32 [%addr_out], %c_proj;
        DONE:
            ret;
        }
        """;

    public static unsafe (bool Ok, string GpuName, double GbpSec, int VerifiedCount) RunRtx4060PtxSelfTest()
    {
        if (cuInit(0) != 0) return (false, "cuInit failed", 0, 0);
        if (cuDeviceGet(out int dev, 0) != 0) return (false, "cuDeviceGet failed", 0, 0);

        Span<byte> nameBuf = stackalloc byte[128];
        cuDeviceGetName(nameBuf, 128, dev);
        int nullIdx = nameBuf.IndexOf((byte)0);
        string gpuName = System.Text.Encoding.UTF8.GetString(nameBuf[..(nullIdx > 0 ? nullIdx : 128)]);

        if (cuCtxCreate_v2(out IntPtr ctx, 0, dev) != 0) return (false, gpuName, 0, 0);
        try
        {
            byte[] ptxBytes = System.Text.Encoding.ASCII.GetBytes(RotatorVerletPtx.Replace("\r\n", "\n") + "\0");
            byte[] errLog = new byte[4096];
            int loadRc;
            IntPtr mod;
            fixed (byte* pPtx = ptxBytes)
            fixed (byte* pErr = errLog)
            {
                int* opts = stackalloc int[2] { 5, 6 };
                void** vals = stackalloc void*[2] { pErr, (void*)4096 };
                loadRc = cuModuleLoadDataEx(out mod, (IntPtr)pPtx, 2, (IntPtr)opts, (IntPtr)vals);
            }
            if (loadRc != 0)
            {
                int errEnd = Array.IndexOf(errLog, (byte)0);
                string jitMsg = System.Text.Encoding.UTF8.GetString(errLog, 0, errEnd > 0 ? errEnd : 256).Trim();
                return (false, $"{gpuName} (PTX err={loadRc}: {jitMsg})", 0, 0);
            }
            if (cuModuleGetFunction(out IntPtr fn, mod, "cuboct_rotator_verlet_ptx") != 0) return (false, gpuName, 0, 0);

            const int N = 16 * 1024 * 1024;
            int[] hostIn = new int[N];
            int[] hostOut = new int[N];
            for (int i = 0; i < N; i++) hostIn[i] = i & 511;

            nuint byteSize = (nuint)(N * sizeof(int));
            cuMemAlloc_v2(out ulong dIn, byteSize);
            cuMemAlloc_v2(out ulong dOut, byteSize);
            try
            {
                fixed (int* pIn = hostIn)
                    cuMemcpyHtoD_v2(dIn, (IntPtr)pIn, byteSize);

                uint nElems = (uint)N;
                ulong argIn = dIn;
                ulong argOut = dOut;
                void** kParams = stackalloc void*[3];
                kParams[0] = &argIn;
                kParams[1] = &argOut;
                kParams[2] = &nElems;

                cuLaunchKernel(fn, (uint)(N / 256), 1, 1, 256, 1, 1, 0, IntPtr.Zero, (IntPtr)kParams, IntPtr.Zero);
                cuCtxSynchronize();

                var sw = Stopwatch.StartNew();
                const int Iters = 20;
                for (int it = 0; it < Iters; it++)
                {
                    cuLaunchKernel(fn, (uint)(N / 256), 1, 1, 256, 1, 1, 0, IntPtr.Zero, (IntPtr)kParams, IntPtr.Zero);
                }
                cuCtxSynchronize();
                sw.Stop();

                fixed (int* pOut = hostOut)
                    cuMemcpyDtoH_v2((IntPtr)pOut, dOut, byteSize);

                int verified = 0;
                for (int i = 0; i < 512; i++)
                {
                    if (hostOut[i] == new Ring9(i).ProjC) verified++;
                }

                double totalBytes = (double)N * sizeof(int) * 2.0 * Iters;
                double gbpSec = (totalBytes / (1024.0 * 1024.0 * 1024.0)) / sw.Elapsed.TotalSeconds;
                return (verified == 512, gpuName, gbpSec, verified);
            }
            finally
            {
                cuMemFree_v2(dIn);
                cuMemFree_v2(dOut);
            }
        }
        finally
        {
            cuCtxDestroy_v2(ctx);
        }
    }
}
