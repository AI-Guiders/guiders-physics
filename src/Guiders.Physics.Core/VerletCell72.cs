using System;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace Guiders.Physics.Core;

[StructLayout(LayoutKind.Sequential, Pack = 1, Size = 9)]
public struct VerletCell72
{
    private ulong _lo64;
    private byte _hi8;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public readonly Ring9 GetLane9(Ring9 laneIdx)
    {
        UInt128 w = (UInt128)_lo64 | ((UInt128)_hi8 << 64);
        ushort bits = (ushort)((w >> ((laneIdx.Bits9 & 7) * 9)) & 0x1FFUL);
        return new Ring9(bits);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void SetLane9(Ring9 laneIdx, Ring9 val)
    {
        ushort shift = (ushort)((laneIdx.Bits9 & 7) * 9);
        UInt128 w = (UInt128)_lo64 | ((UInt128)_hi8 << 64);
        UInt128 mask = ~((UInt128)0x1FFUL << shift);
        w = (w & mask) | ((UInt128)val.Bits9 << shift);
        _lo64 = (ulong)w;
        _hi8 = (byte)(w >> 64);
    }

    public Ring9 Mu0      { readonly get => GetLane9(new(0)); set => SetLane9(new(0), value); }
    public Ring9 Phi0Curr { readonly get => GetLane9(new(1)); set => SetLane9(new(1), value); }
    public Ring9 Phi0Past { readonly get => GetLane9(new(2)); set => SetLane9(new(2), value); }
    public Ring9 Mu1      { readonly get => GetLane9(new(3)); set => SetLane9(new(3), value); }
    public Ring9 Phi1Curr { readonly get => GetLane9(new(4)); set => SetLane9(new(4), value); }
    public Ring9 Phi1Past { readonly get => GetLane9(new(5)); set => SetLane9(new(5), value); }
    public Ring9 Gate0    { readonly get => GetLane9(new(6)); set => SetLane9(new(6), value); }
    public Ring9 Gate1    { readonly get => GetLane9(new(7)); set => SetLane9(new(7), value); }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static VerletCell72 Create(Ring9 mu0, Ring9 p0Curr, Ring9 p0Past, Ring9 mu1, Ring9 p1Curr, Ring9 p1Past)
    {
        VerletCell72 c = default;
        c.Mu0 = mu0;
        c.Phi0Curr = p0Curr;
        c.Phi0Past = p0Past;
        c.Mu1 = mu1;
        c.Phi1Curr = p1Curr;
        c.Phi1Past = p1Past;
        return c;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public readonly VerletCell72 EvaluateClosedFormAtTick(Ring9 t)
    {
        Ring9 v0_0 = Phi0Curr - Phi0Past;
        Ring9 v0_1 = Phi1Curr - Phi1Past;
        Ring9 tPrev = t - Ring9.MuVac;
        Ring9 p0AtT = Phi0Curr.JumpVerletParabola(v0_0, Gate0, t);
        Ring9 p0AtTMinus1 = Phi0Curr.JumpVerletParabola(v0_0, Gate0, tPrev);
        Ring9 p1AtT = Phi1Curr.JumpVerletParabola(v0_1, Gate1, t);
        Ring9 p1AtTMinus1 = Phi1Curr.JumpVerletParabola(v0_1, Gate1, tPrev);
        var c = Create(Mu0, p0AtT, p0AtTMinus1, Mu1, p1AtT, p1AtTMinus1);
        c.Gate0 = Gate0;
        c.Gate1 = Gate1;
        return c;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public readonly VerletCell72 StepForward(Ring9 gate) =>
        CuboctahedronAutomatonKernel.StepCellForward(this, gate, gate);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public readonly VerletCell72 StepReverse(Ring9 gateAtT) =>
        CuboctahedronAutomatonKernel.StepCellReverse(this, gateAtT, gateAtT);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public readonly MultimodalReadout9 ReadSensors() =>
        CuboctahedronAutomatonKernel.ReadCellSensors(this);
}

public readonly record struct MultimodalReadout9(
    Ring9 ColorR,
    Ring9 ColorG,
    Ring9 ColorB,
    Ring9 AcousticWave,
    Ring9 ReactionCurrent,
    Ring9 EnergyQuantaNe,
    Ring9 OlfactoryShellZ,
    Ring9 ChiralOdorProjectionS);
