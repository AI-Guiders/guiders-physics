using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Runtime.InteropServices;
using System.Windows.Forms;
using Guiders.Physics.Core;
using Guiders.Physics.Engine;

namespace Guiders.Physics.ChronovisorApp;

public sealed unsafe class DirectWaveAudioOut : IDisposable
{
    private const int SampleRate = 8000;
    private const int Channels = 2;
    private const int BufferSamples = 256;
    private const int NumBuffers = 4;

    [StructLayout(LayoutKind.Sequential)]
    private struct WaveFormatEx
    {
        public ushort wFormatTag;
        public ushort nChannels;
        public uint nSamplesPerSec;
        public uint nAvgBytesPerSec;
        public ushort nBlockAlign;
        public ushort wBitsPerSample;
        public ushort cbSize;
    }

    [StructLayout(LayoutKind.Sequential)]
    private struct WaveHdr
    {
        public IntPtr lpData;
        public uint dwBufferLength;
        public uint dwBytesRecorded;
        public IntPtr dwUser;
        public uint dwFlags;
        public uint dwLoops;
        public IntPtr lpNext;
        public IntPtr reserved;
    }

    [DllImport("winmm.dll")]
    private static extern int waveOutOpen(out IntPtr hWaveOut, int uDeviceID, ref WaveFormatEx lpFormat, IntPtr dwCallback, IntPtr dwInstance, uint dwFlags);

    [DllImport("winmm.dll")]
    private static extern int waveOutPrepareHeader(IntPtr hWaveOut, IntPtr lpWaveOutHdr, uint uSize);

    [DllImport("winmm.dll")]
    private static extern int waveOutWrite(IntPtr hWaveOut, IntPtr lpWaveOutHdr, uint uSize);

    [DllImport("winmm.dll")]
    private static extern int waveOutReset(IntPtr hWaveOut);

    [DllImport("winmm.dll")]
    private static extern int waveOutClose(IntPtr hWaveOut);

    private IntPtr _hWaveOut;
    private readonly IntPtr[] _hdrPtrs = new IntPtr[NumBuffers];
    private readonly IntPtr[] _bufPtrs = new IntPtr[NumBuffers];
    private int _curBuf;
    public bool IsOpen { get; private set; }

    public DirectWaveAudioOut()
    {
        var fmt = new WaveFormatEx
        {
            wFormatTag = 1,
            nChannels = Channels,
            nSamplesPerSec = SampleRate,
            wBitsPerSample = 16,
            nBlockAlign = (ushort)(Channels * 2),
            nAvgBytesPerSec = SampleRate * Channels * 2,
            cbSize = 0
        };

        if (waveOutOpen(out _hWaveOut, -1, ref fmt, IntPtr.Zero, IntPtr.Zero, 0) == 0)
        {
            IsOpen = true;
            uint byteLen = BufferSamples * Channels * 2;
            for (int i = 0; i < NumBuffers; i++)
            {
                _bufPtrs[i] = Marshal.AllocHGlobal((int)byteLen);
                _hdrPtrs[i] = Marshal.AllocHGlobal(sizeof(WaveHdr));
                var hdr = new WaveHdr
                {
                    lpData = _bufPtrs[i],
                    dwBufferLength = byteLen,
                    dwFlags = 0
                };
                Marshal.StructureToPtr(hdr, _hdrPtrs[i], false);
                waveOutPrepareHeader(_hWaveOut, _hdrPtrs[i], (uint)sizeof(WaveHdr));
            }
        }
    }

    public void PushRawModelAudio(
        CuboctWorldScene3D world,
        (Ring9 Z, Ring9 Y, Ring9 X) leftEar,
        (Ring9 Z, Ring9 Y, Ring9 X) rightEar,
        Ring9 baseTau)
    {
        if (!IsOpen) return;

        IntPtr hdrPtr = _hdrPtrs[_curBuf];
        var hdr = Marshal.PtrToStructure<WaveHdr>(hdrPtr);
        if ((hdr.dwFlags & 0x10) != 0)
        {
            _curBuf = (_curBuf + 1) % NumBuffers;
            hdrPtr = _hdrPtrs[_curBuf];
            hdr = Marshal.PtrToStructure<WaveHdr>(hdrPtr);
            if ((hdr.dwFlags & 0x10) != 0) return;
        }

        short* pcm = (short*)_bufPtrs[_curBuf];
        for (ushort i = 0; i < BufferSamples; i++)
        {
            Ring9 t = baseTau + new Ring9(i);
            var cellL = world.EvaluateWorldlineAt(leftEar.Z, leftEar.Y, leftEar.X, t);
            var cellR = world.EvaluateWorldlineAt(rightEar.Z, rightEar.Y, rightEar.X, t);
            var readL = CuboctahedronAutomatonKernel.ReadCellSensors(cellL);
            var readR = CuboctahedronAutomatonKernel.ReadCellSensors(cellR);

            pcm[2 * i] = (short)((readL.AcousticWave.Bits9 - 256) << 5);
            pcm[2 * i + 1] = (short)((readR.AcousticWave.Bits9 - 256) << 5);
        }

        waveOutWrite(_hWaveOut, hdrPtr, (uint)sizeof(WaveHdr));
        _curBuf = (_curBuf + 1) % NumBuffers;
    }

    public void Dispose()
    {
        if (IsOpen)
        {
            IsOpen = false;
            waveOutReset(_hWaveOut);
            waveOutClose(_hWaveOut);
            for (int i = 0; i < NumBuffers; i++)
            {
                Marshal.FreeHGlobal(_bufPtrs[i]);
                Marshal.FreeHGlobal(_hdrPtrs[i]);
            }
        }
    }
}

public sealed class ChronovisorForm : Form
{
    private const int ViewW = 96;
    private const int ViewH = 72;

    private readonly CuboctWorldScene3D _world = new(WorldClockMode.AtomSynced);
    private readonly Bitmap _bmp = new(ViewW, ViewH, PixelFormat.Format32bppRgb);
    private readonly int[] _pixels = new int[ViewW * ViewH];
    private readonly System.Windows.Forms.Timer _timer = new() { Interval = 33 };
    private readonly DirectWaveAudioOut _audio = new();

    private Ring9 _camX = Ring9.Zero;
    private Ring9 _camY = Ring9.Zero;
    private Ring9 _camZ = Ring9.Zero;
    private Ring9 _tauOffset = Ring9.Zero;
    private bool _soundOn = true;

    public ChronovisorForm()
    {
        Text = "Guiders.Physics — Pure 3D ReadCellSensors Output (Zero Hand-Drawn Boosts / Zero Synthetic Audio)";
        ClientSize = new Size(960, 720);
        StartPosition = FormStartPosition.CenterScreen;
        DoubleBuffered = true;
        KeyPreview = true;
        BackColor = Color.Black;

        _world.SyncWithSiliconThermalTrng();
        Render3DAndPushRawAudio();

        _timer.Tick += (_, _) =>
        {
            _world.AdvanceClock();
            Render3DAndPushRawAudio();
            Invalidate();
        };
        _timer.Start();
    }

    protected override void OnFormClosed(FormClosedEventArgs e)
    {
        _timer.Stop();
        _audio.Dispose();
        base.OnFormClosed(e);
    }

    protected override void OnKeyDown(KeyEventArgs e)
    {
        base.OnKeyDown(e);
        switch (e.KeyCode)
        {
            case Keys.Space:
                _camX = Ring9.Zero;
                _camY = Ring9.Zero;
                _camZ = Ring9.Zero;
                _tauOffset = Ring9.Zero;
                _world.SyncWithSiliconThermalTrng();
                break;
            case Keys.Left:
            case Keys.A:
                _camX = _camX - Ring9.MuVac;
                break;
            case Keys.Right:
            case Keys.D:
                _camX = _camX + Ring9.MuVac;
                break;
            case Keys.Up:
            case Keys.W:
                _camZ = _camZ + Ring9.MuVac;
                break;
            case Keys.Down:
            case Keys.S:
                _camZ = _camZ - Ring9.MuVac;
                break;
            case Keys.Q:
                _camY = _camY - Ring9.MuVac;
                break;
            case Keys.E:
                _camY = _camY + Ring9.MuVac;
                break;
            case Keys.OemOpenBrackets:
            case Keys.PageDown:
                _tauOffset = _tauOffset - Ring9.MuVac;
                break;
            case Keys.OemCloseBrackets:
            case Keys.PageUp:
                _tauOffset = _tauOffset + Ring9.MuVac;
                break;
            case Keys.M:
                _world.ClockMode = _world.ClockMode == WorldClockMode.AtomSynced
                    ? WorldClockMode.FreeRunning
                    : WorldClockMode.AtomSynced;
                break;
            case Keys.V:
                _soundOn = !_soundOn;
                break;
        }
        Render3DAndPushRawAudio();
        Invalidate();
    }

    private void Render3DAndPushRawAudio()
    {
        Ring9 targetTau = _world.CurrentTick + _tauOffset;
        int halfW = ViewW / 2;
        int halfH = ViewH / 2;

        // 3D depth integration across all 3 dimensions (Z, Y, X) using purely ReadCellSensors(cell) output:
        for (int py = 0; py < ViewH; py++)
        {
            Ring9 dy = _camY + new Ring9((ushort)((halfH - py) & 0x1FF));
            int rowOff = py * ViewW;
            for (int px = 0; px < ViewW; px++)
            {
                Ring9 dx = _camX + new Ring9((ushort)((px - halfW) & 0x1FF));
                var cell = _world.EvaluateWorldlineAt(_camZ, dy, dx, targetTau);
                var readout = CuboctahedronAutomatonKernel.ReadCellSensors(cell);

                // Direct unmodified readout.ColorR, ColorG, ColorB from ReadCellSensors:
                int r = readout.ColorR.Bits9;
                int g = readout.ColorG.Bits9;
                int b = readout.ColorB.Bits9;
                _pixels[rowOff + px] = (255 << 24) | (r << 16) | (g << 8) | b;
            }
        }

        var rect = new Rectangle(0, 0, ViewW, ViewH);
        var bmpData = _bmp.LockBits(rect, ImageLockMode.WriteOnly, PixelFormat.Format32bppRgb);
        Marshal.Copy(_pixels, 0, bmpData.Scan0, _pixels.Length);
        _bmp.UnlockBits(bmpData);

        if (_soundOn)
        {
            _audio.PushRawModelAudio(
                _world,
                (_camZ, _camY, _camX - Ring9.MuVac),
                (_camZ, _camY, _camX + Ring9.MuVac),
                targetTau);
        }
    }

    protected override void OnPaint(PaintEventArgs e)
    {
        var g = e.Graphics;
        g.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.NearestNeighbor;
        g.PixelOffsetMode = System.Drawing.Drawing2D.PixelOffsetMode.Half;
        g.DrawImage(_bmp, new Rectangle(0, 0, ClientSize.Width, ClientSize.Height));
    }

    [STAThread]
    public static int Main(string[] args)
    {
        // Universal domain-agnostic 12-neighbor cuboctahedron cell evaluator:
        // Usage: --cell <mu0> <phi0> <mu1> <phi1> <nbMu> <nbPhi> <t>
        if (args.Length >= 8 && args[0] == "--cell")
        {
            Ring9 mu0 = new(ushort.Parse(args[1]));
            Ring9 phi0 = new(ushort.Parse(args[2]));
            Ring9 mu1 = new(ushort.Parse(args[3]));
            Ring9 phi1 = new(ushort.Parse(args[4]));
            Ring9 nbMu = new(ushort.Parse(args[5]));
            Ring9 nbPhi = new(ushort.Parse(args[6]));
            Ring9 t = new(ushort.Parse(args[7]));

            Ring9 gate0 = CuboctahedronAutomatonKernel.Compute12NeighborGateFromCells(
                mu0, phi0,
                nbMu, nbPhi, nbMu, nbPhi, nbMu, nbPhi, nbMu, nbPhi,
                nbMu, nbPhi, nbMu, nbPhi, nbMu, nbPhi, nbMu, nbPhi,
                nbMu, nbPhi, nbMu, nbPhi, nbMu, nbPhi, nbMu, nbPhi);
            Ring9 pauli = CuboctahedronAutomatonKernel.PauliExclusionGate9(mu0, phi0, mu1, phi1, Ring9.MuVac);

            var cell = VerletCell72.Create(mu0, phi0, phi0, mu1, phi1, phi1);
            cell.Gate0 = gate0;
            cell.Gate1 = gate0;
            var stepped = cell.EvaluateClosedFormAtTick(t);
            var sensors = CuboctahedronAutomatonKernel.ReadCellSensors(stepped);

            Console.WriteLine(
                $"{{\"Gate0\":{gate0.Bits9}," +
                $"\"PauliGate\":{pauli.Bits9}," +
                $"\"Phi0AtT\":{stepped.Phi0Curr.Bits9}," +
                $"\"ProjC\":{stepped.Phi0Curr.ProjC.Bits9}," +
                $"\"ProjS\":{stepped.Phi0Curr.ProjS.Bits9}," +
                $"\"ShellZ\":{sensors.OlfactoryShellZ.Bits9}," +
                $"\"AcousticWave\":{sensors.AcousticWave.Bits9}," +
                $"\"ReactionCurrent\":{sensors.ReactionCurrent.Bits9}," +
                $"\"EnergyQuantaNe\":{sensors.EnergyQuantaNe.Bits9}}}");
            return 0;
        }

        ApplicationConfiguration.Initialize();
        Application.Run(new ChronovisorForm());
        return 0;
    }
}
