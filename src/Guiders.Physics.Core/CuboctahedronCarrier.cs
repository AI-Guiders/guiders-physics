using System.Runtime.CompilerServices;

namespace Guiders.Physics.Core;

public sealed class CuboctahedronCarrier
{
    private VerletCell72[,,] _grid;
    private VerletCell72[,,] _scratch;
    private readonly ushort _n;

    public Ring9 Size => new(_n);
    public Ring9 Tick { get; private set; }

    public CuboctahedronCarrier(Ring9 size)
    {
        _n = size.Bits9;
        _grid = new VerletCell72[_n, _n, _n];
        _scratch = new VerletCell72[_n, _n, _n];
        Tick = Ring9.Zero;
        ushort c = (ushort)(_n >> 1);
        Ring9 center = new(c);
        for (ushort iz = 0; iz < _n; iz++)
            for (ushort iy = 0; iy < _n; iy++)
                for (ushort ix = 0; ix < _n; ix++)
                {
                    Ring9 dz = new Ring9(iz) - center;
                    Ring9 dy = new Ring9(iy) - center;
                    Ring9 dx = new Ring9(ix) - center;
                    Ring9 mu = CuboctahedronAutomatonKernel.Iv263MuAt(dz, dy, dx);
                    Ring9 p0 = CuboctahedronAutomatonKernel.Iv263PhiAt(dz, dy, dx, Ring9.Zero, center);
                    Ring9 p1 = CuboctahedronAutomatonKernel.Iv263PhiAt(dz, dy, dx, Ring9.MuVac, center);
                    _grid[iz, iy, ix] = VerletCell72.Create(mu, p0, p0, mu, p1, p1);
                }
    }

    public VerletCell72 this[Ring9 z, Ring9 y, Ring9 x]
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _grid[(z.Bits9 + _n) % _n, (y.Bits9 + _n) % _n, (x.Bits9 + _n) % _n];
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set => _grid[(z.Bits9 + _n) % _n, (y.Bits9 + _n) % _n, (x.Bits9 + _n) % _n] = value;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public MultimodalReadout9 ReadSensors(Ring9 z, Ring9 y, Ring9 x) =>
        this[z, y, x].ReadSensors();

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private VerletCell72 CellAt(ushort z, ushort y, ushort x) =>
        _grid[z % _n, y % _n, x % _n];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private Ring9 ComputeGateAtNode(ushort z, ushort y, ushort x, Ring9 usePastFlag)
    {
        ushort n = _n;
        ushort zm = (ushort)(z + n - 1), zp = (ushort)(z + 1);
        ushort ym = (ushort)(y + n - 1), yp = (ushort)(y + 1);
        ushort xm = (ushort)(x + n - 1), xp = (ushort)(x + 1);

        var c   = CellAt(z,  y,  x);
        var n0  = CellAt(z,  ym, xm);
        var n1  = CellAt(z,  ym, xp);
        var n2  = CellAt(z,  yp, xm);
        var n3  = CellAt(z,  yp, xp);
        var n4  = CellAt(zm, y,  xm);
        var n5  = CellAt(zm, y,  xp);
        var n6  = CellAt(zp, y,  xm);
        var n7  = CellAt(zp, y,  xp);
        var n8  = CellAt(zm, ym, x);
        var n9  = CellAt(zm, yp, x);
        var n10 = CellAt(zp, ym, x);
        var n11 = CellAt(zp, yp, x);

        Ring9 selPast = usePastFlag;
        Ring9 selCurr = Ring9.MuVac - selPast;

        Ring9 P0(VerletCell72 v) => selCurr * v.Phi0Curr + selPast * v.Phi0Past;
        Ring9 P1(VerletCell72 v) => selCurr * v.Phi1Curr + selPast * v.Phi1Past;

        return CuboctahedronAutomatonKernel.ComputeTwoLaneCellGate(
            c.Mu0, P0(c), c.Mu1, P1(c),
            n0.Mu0,  P0(n0),  n1.Mu0,  P0(n1),  n2.Mu0,  P0(n2),  n3.Mu0,  P0(n3),
            n4.Mu0,  P0(n4),  n5.Mu0,  P0(n5),  n6.Mu0,  P0(n6),  n7.Mu0,  P0(n7),
            n8.Mu0,  P0(n8),  n9.Mu0,  P0(n9),  n10.Mu0, P0(n10), n11.Mu0, P0(n11),
            n0.Mu1,  P1(n0),  n1.Mu1,  P1(n1),  n2.Mu1,  P1(n2),  n3.Mu1,  P1(n3),
            n4.Mu1,  P1(n4),  n5.Mu1,  P1(n5),  n6.Mu1,  P1(n6),  n7.Mu1,  P1(n7),
            n8.Mu1,  P1(n8),  n9.Mu1,  P1(n9),  n10.Mu1, P1(n10), n11.Mu1, P1(n11),
            Ring9.MuVac);
    }

    public void StepForward()
    {
        ushort n = _n;
        for (ushort z = 0; z < n; z++)
            for (ushort y = 0; y < n; y++)
                for (ushort x = 0; x < n; x++)
                {
                    Ring9 gate = ComputeGateAtNode(z, y, x, Ring9.Zero);
                    _scratch[z, y, x] = _grid[z, y, x].StepForward(gate);
                }
        (_grid, _scratch) = (_scratch, _grid);
        Tick = Tick + Ring9.MuVac;
    }

    public void StepBackward()
    {
        ushort n = _n;
        for (ushort z = 0; z < n; z++)
            for (ushort y = 0; y < n; y++)
                for (ushort x = 0; x < n; x++)
                {
                    Ring9 gateAtT = ComputeGateAtNode(z, y, x, Ring9.MuVac);
                    _scratch[z, y, x] = _grid[z, y, x].StepReverse(gateAtT);
                }
        (_grid, _scratch) = (_scratch, _grid);
        Tick = Tick - Ring9.MuVac;
    }
}
