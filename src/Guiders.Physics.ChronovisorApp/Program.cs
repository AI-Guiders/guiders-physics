using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Runtime.InteropServices;
using System.Windows.Forms;
using Guiders.Physics.Core;
using Guiders.Physics.Engine;

namespace Guiders.Physics.ChronovisorApp;

public sealed unsafe class LiveWinMmStereoAudioOut : IDisposable
{
    private const int SampleRate = 22050;
    private const int Channels = 2;
    private const int BufferSamples = 735;
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
    private Ring9 _carrierPhaseL = Ring9.Zero;
    private Ring9 _carrierPhaseR = Ring9.Zero;
    public bool IsOpen { get; private set; }

    public LiveWinMmStereoAudioOut()
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

    public void PushDirectReadoutAudio(MultimodalReadout9 leftReadout, MultimodalReadout9 rightReadout)
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

        Ring9 stepL = Ring9.DOutVertex + new Ring9((ushort)(leftReadout.AcousticWave.Bits9 >> 6));
        Ring9 stepR = Ring9.DOutVertex + new Ring9((ushort)(rightReadout.AcousticWave.Bits9 >> 6));

        short* pcm = (short*)_bufPtrs[_curBuf];
        for (int i = 0; i < BufferSamples; i++)
        {
            _carrierPhaseL = _carrierPhaseL + stepL;
            _carrierPhaseR = _carrierPhaseR + stepR;

            short waveL = (short)((_carrierPhaseL.ProjS.Bits9 - 256) * 18);
            short waveR = (short)((_carrierPhaseR.ProjS.Bits9 - 256) * 18);
            pcm[2 * i] = waveL;
            pcm[2 * i + 1] = waveR;
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
    private const int ViewW = 48;
    private const int ViewH = 36;

    private readonly CuboctWorldScene3D _world = new(WorldClockMode.AtomSynced);
    private readonly Bitmap _bmp = new(ViewW, ViewH, PixelFormat.Format32bppRgb);
    private readonly int[] _pixels = new int[ViewW * ViewH];
    private readonly System.Windows.Forms.Timer _timer = new() { Interval = 33 };
    private readonly LiveWinMmStereoAudioOut _audio = new();

    private Ring9 _camX = Ring9.Zero;
    private Ring9 _camY = Ring9.Zero;
    private Ring9 _camZ = Ring9.Zero;
    private Ring9 _tauOffset = Ring9.Zero;
    private bool _soundOn = true;

    public ChronovisorForm()
    {
        Text = "Guiders.Physics — Direct ReadCellSensors (MultimodalReadout9 -> Screen RGB & WinMM Stereo Audio)";
        ClientSize = new Size(960, 720);
        StartPosition = FormStartPosition.CenterScreen;
        DoubleBuffered = true;
        KeyPreview = true;
        BackColor = Color.Black;

        _world.SyncWithSiliconThermalTrng();
        RenderViewportAndPushAudio();

        _timer.Tick += (_, _) =>
        {
            _world.AdvanceClock();
            RenderViewportAndPushAudio();
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
        RenderViewportAndPushAudio();
        Invalidate();
    }

    private void RenderViewportAndPushAudio()
    {
        Ring9 targetTau = _world.CurrentTick + _tauOffset;
        int halfW = ViewW / 2;
        int halfH = ViewH / 2;

        for (int py = 0; py < ViewH; py++)
        {
            Ring9 dy = _camY + new Ring9((ushort)((halfH - py) & 0x1FF));
            int rowOff = py * ViewW;
            for (int px = 0; px < ViewW; px++)
            {
                Ring9 dx = _camX + new Ring9((ushort)((px - halfW) & 0x1FF));
                var cell = _world.EvaluateWorldlineAt(_camZ, dy, dx, targetTau);
                var readout = CuboctahedronAutomatonKernel.ReadCellSensors(cell);

                int boost = cell.Mu0.Bits9 > 1 ? 2 : 1;
                int r = Math.Clamp(readout.ColorR.Bits9 * boost, 0, 255);
                int g = Math.Clamp(readout.ColorG.Bits9 * boost, 0, 255);
                int b = Math.Clamp(readout.ColorB.Bits9 * boost, 0, 255);
                _pixels[rowOff + px] = (255 << 24) | (r << 16) | (g << 8) | b;
            }
        }

        var rect = new Rectangle(0, 0, ViewW, ViewH);
        var bmpData = _bmp.LockBits(rect, ImageLockMode.WriteOnly, PixelFormat.Format32bppRgb);
        Marshal.Copy(_pixels, 0, bmpData.Scan0, _pixels.Length);
        _bmp.UnlockBits(bmpData);

        if (_soundOn)
        {
            var leftReadout = CuboctahedronAutomatonKernel.ReadCellSensors(
                _world.EvaluateWorldlineAt(_camZ, _camY, _camX - Ring9.C0OverCSq, targetTau));
            var rightReadout = CuboctahedronAutomatonKernel.ReadCellSensors(
                _world.EvaluateWorldlineAt(_camZ, _camY, _camX + Ring9.C0OverCSq, targetTau));
            _audio.PushDirectReadoutAudio(leftReadout, rightReadout);
        }
    }

    protected override void OnPaint(PaintEventArgs e)
    {
        var g = e.Graphics;
        int viewH = ClientSize.Height - 115;
        g.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.NearestNeighbor;
        g.PixelOffsetMode = System.Drawing.Drawing2D.PixelOffsetMode.Half;
        g.DrawImage(_bmp, new Rectangle(0, 0, ClientSize.Width, viewH));

        Ring9 targetTau = _world.CurrentTick + _tauOffset;
        var centerCell = _world.EvaluateWorldlineAt(_camZ, _camY, _camX, targetTau);
        var centerReadout = CuboctahedronAutomatonKernel.ReadCellSensors(centerCell);

        using var hudBrush = new SolidBrush(Color.FromArgb(230, 12, 16, 24));
        g.FillRectangle(hudBrush, 0, viewH, ClientSize.Width, 115);

        using var font = new Font("Consolas", 10.5f, FontStyle.Bold);
        string l1 = $"[CLOCK MODE: {_world.ClockMode} (M)]  CurrentTick={_world.CurrentTick}  tauOffset={_tauOffset}  ->  TargetTau={targetTau}";
        string l2 = $"[CAMERA POS]  Z={_camZ} (W/S)  Y={_camY} (Q/E)  X={_camX} (A/D)  |  Cell Mu={centerCell.Mu0}  Gate={centerCell.Gate0}";
        string l3 = $"[DIRECT READOUT]  RGB=({centerReadout.ColorR.Bits9},{centerReadout.ColorG.Bits9},{centerReadout.ColorB.Bits9}) sum=128 | Acoustic={centerReadout.AcousticWave.Bits9} | OdorZ={centerReadout.OlfactoryShellZ.Bits9} | Quanta n_E={centerReadout.EnergyQuantaNe.Bits9} | Audio={(_soundOn ? "ON (V)" : "OFF (V)")}";
        string l4 = "CONTROLS: [Space]=HERE & NOW | [WASD/Arrows]=Move X/Z | [Q/E]=Move Y | [ [ / ] ]=Time Knob dTau | [M]=AtomSynced(RDSEED)/FreeRunning | [V]=Mute";

        g.DrawString(l1, font, Brushes.Cyan, 12, viewH + 8);
        g.DrawString(l2, font, Brushes.Yellow, 12, viewH + 34);
        g.DrawString(l3, font, Brushes.LightGreen, 12, viewH + 60);
        g.DrawString(l4, font, Brushes.White, 12, viewH + 86);
    }

    [STAThread]
    public static void Main()
    {
        ApplicationConfiguration.Initialize();
        Application.Run(new ChronovisorForm());
    }
}
