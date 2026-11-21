using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Runtime.InteropServices;
using System.Windows.Forms;
using Guiders.Physics.Core;
using Guiders.Physics.Engine;

namespace Guiders.Physics.ChronovisorApp;

public sealed class ChronovisorForm : Form
{
    private const int ViewW = 256;
    private const int ViewH = 192;

    private readonly CuboctWorldScene3D _world = new(WorldClockMode.AtomSynced);
    private readonly Bitmap _bmp = new(ViewW, ViewH, PixelFormat.Format32bppRgb);
    private readonly int[] _pixels = new int[ViewW * ViewH];
    private readonly System.Windows.Forms.Timer _timer = new() { Interval = 33 };

    private int _offsetX;
    private int _offsetY;
    private int _offsetZ;
    private int _offsetTau;
    private bool _freezeAnchor;
    private bool _audioEnabled;

    private Ring9 _lastPhiSi = Ring9.IvSeed;
    private Ring9 _lastDphiSi = Ring9.Zero;
    private Ring9 _anchorXNow = Ring9.Zero;
    private Ring9 _anchorTauNow = Ring9.Zero;

    [DllImport("winmm.dll")]
    private static extern bool Beep(uint dwFreq, uint dwDuration);

    public ChronovisorForm()
    {
        Text = "Guiders.Physics — Cuboctahedron Chronovisor (RDSEED 48 0F C7 F8 Silicon TRNG -> Ring9 -> RGB/Audio)";
        ClientSize = new Size(960, 720);
        StartPosition = FormStartPosition.CenterScreen;
        DoubleBuffered = true;
        KeyPreview = true;
        BackColor = Color.Black;

        SampleFreshSiliconAnchor();

        _timer.Tick += (_, _) =>
        {
            if (!_freezeAnchor)
                SampleFreshSiliconAnchor();
            RenderChronovisorViewport();
            Invalidate();
        };
        _timer.Start();
    }

    private void SampleFreshSiliconAnchor()
    {
        var (phiNow, dphiVerlet) = CuboctahedronAutomatonKernel.SampleSiliconThermalTwoTickPair9();
        _lastPhiSi = phiNow;
        _lastDphiSi = dphiVerlet;
        _world.SyncWithPhysicalAtom(phiNow, dphiVerlet);
        _anchorXNow = _world.SpatialAnchorX;
        _anchorTauNow = _world.CurrentTick;
    }

    protected override void OnKeyDown(KeyEventArgs e)
    {
        base.OnKeyDown(e);
        switch (e.KeyCode)
        {
            case Keys.Space:
                _offsetX = 0;
                _offsetY = 0;
                _offsetZ = 0;
                _offsetTau = 0;
                SampleFreshSiliconAnchor();
                break;
            case Keys.Left:
            case Keys.A:
                _offsetX = (_offsetX - 1) & 0x1FF;
                break;
            case Keys.Right:
            case Keys.D:
                _offsetX = (_offsetX + 1) & 0x1FF;
                break;
            case Keys.Up:
            case Keys.W:
                _offsetZ = (_offsetZ + 1) & 0x1FF;
                break;
            case Keys.Down:
            case Keys.S:
                _offsetZ = (_offsetZ - 1) & 0x1FF;
                break;
            case Keys.Q:
                _offsetY = (_offsetY - 1) & 0x1FF;
                break;
            case Keys.E:
                _offsetY = (_offsetY + 1) & 0x1FF;
                break;
            case Keys.OemOpenBrackets:
            case Keys.PageDown:
                _offsetTau = (_offsetTau - 1) & 0x1FF;
                break;
            case Keys.OemCloseBrackets:
            case Keys.PageUp:
                _offsetTau = (_offsetTau + 1) & 0x1FF;
                break;
            case Keys.M:
                _freezeAnchor = !_freezeAnchor;
                break;
            case Keys.F:
                _world.EmitGameEvent(
                    new Ring9((ushort)_offsetZ),
                    new Ring9((ushort)_offsetY),
                    new Ring9((ushort)_offsetX),
                    Ring9.KP,
                    Ring9.CollapsePhi);
                break;
            case Keys.C:
                _world.ClearGameEvents();
                break;
            case Keys.V:
                _audioEnabled = !_audioEnabled;
                break;
        }
        RenderChronovisorViewport();
        Invalidate();
    }

    private void RenderChronovisorViewport()
    {
        Ring9 targetTau = new((ushort)((_anchorTauNow.Bits9 + _offsetTau) & 0x1FF));
        Ring9 baseZ = new((ushort)(_offsetZ & 0x1FF));
        int halfW = ViewW / 2;
        int halfH = ViewH / 2;

        for (int py = 0; py < ViewH; py++)
        {
            Ring9 dy = new((ushort)((_offsetY + (halfH - py)) & 0x1FF));
            int rowOff = py * ViewW;
            for (int px = 0; px < ViewW; px++)
            {
                Ring9 dx = new((ushort)((_anchorXNow.Bits9 + _offsetX + (px - halfW)) & 0x1FF));
                Ring9 phase = CuboctahedronAutomatonKernel.ChronovisorArbitraryTargetPhase(
                    _lastPhiSi, _lastDphiSi, dx + dy + baseZ, targetTau);

                var cell = _world.EvaluateWorldlineAt(baseZ, dy, dx, targetTau);
                cell.Phi1Curr = cell.Phi0Curr + phase;
                var sensor = CuboctahedronAutomatonKernel.ReadCellSensors(cell);

                int r = Math.Clamp(sensor.ColorR.Bits9 * 2, 0, 255);
                int g = Math.Clamp(sensor.ColorG.Bits9 * 2, 0, 255);
                int b = Math.Clamp(sensor.ColorB.Bits9 * 2, 0, 255);
                _pixels[rowOff + px] = (255 << 24) | (r << 16) | (g << 8) | b;
            }
        }

        var rect = new Rectangle(0, 0, ViewW, ViewH);
        var bmpData = _bmp.LockBits(rect, ImageLockMode.WriteOnly, PixelFormat.Format32bppRgb);
        Marshal.Copy(_pixels, 0, bmpData.Scan0, _pixels.Length);
        _bmp.UnlockBits(bmpData);
    }

    protected override void OnPaint(PaintEventArgs e)
    {
        var g = e.Graphics;
        g.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.NearestNeighbor;
        g.PixelOffsetMode = System.Drawing.Drawing2D.PixelOffsetMode.Half;
        g.DrawImage(_bmp, new Rectangle(0, 0, ClientSize.Width, ClientSize.Height - 110));

        Ring9 targetTau = new((ushort)((_anchorTauNow.Bits9 + _offsetTau) & 0x1FF));
        var centerSensor = _world.SampleCameraVoxel(
            new Ring9((ushort)_offsetZ),
            new Ring9((ushort)_offsetY),
            new Ring9((ushort)(_anchorXNow.Bits9 + _offsetX)));

        using var hudBrush = new SolidBrush(Color.FromArgb(220, 12, 16, 24));
        g.FillRectangle(hudBrush, 0, ClientSize.Height - 110, ClientSize.Width, 110);

        using var font = new Font("Consolas", 10.5f, FontStyle.Bold);
        string modeStr = _freezeAnchor ? "LOCKED-SNAPSHOT (Press M to live-stream RDSEED)" : "LIVE RDSEED SILICON TRNG (48 0F C7 F8)";
        string line1 = $"[MODE: {modeStr}]  Silicon Anchor: phi_Si={_lastPhiSi}  dphi_Si={_lastDphiSi}  ->  x_now={_anchorXNow.Bits9}, tau_now={_anchorTauNow.Bits9}";
        string line2 = $"[VIEWPORT TARGET]  dX={_offsetX}  dY={_offsetY}  dZ(Depth)={_offsetZ}  |  dTau(Time Knob)={_offsetTau}  ->  Target Tick tau={targetTau}";
        string line3 = $"[CENTER SENSOR]    RGB=({centerSensor.ColorR.Bits9},{centerSensor.ColorG.Bits9},{centerSensor.ColorB.Bits9}) sum=128 | Acoustic={centerSensor.AcousticWave.Bits9} | OdorZ={centerSensor.OlfactoryShellZ.Bits9} | Quanta n_E={centerSensor.EnergyQuantaNe.Bits9}";
        string line4 = "CONTROLS: [Space]=HERE & NOW (0,0,0,0) | [Arrows/WASD]=Move X/Z | [Q/E]=Move Y | [ [ / ] ]=Rewind/Forward Time dTau | [M]=Freeze/Live RDSEED | [F]=Emit Impulse";

        g.DrawString(line1, font, Brushes.Cyan, 12, ClientSize.Height - 102);
        g.DrawString(line2, font, Brushes.Yellow, 12, ClientSize.Height - 78);
        g.DrawString(line3, font, Brushes.LightGreen, 12, ClientSize.Height - 54);
        g.DrawString(line4, font, Brushes.White, 12, ClientSize.Height - 28);
    }

    [STAThread]
    public static void Main()
    {
        ApplicationConfiguration.Initialize();
        Application.Run(new ChronovisorForm());
    }
}
