# AIGuiders.Physics.Engine

Unified 3D world simulation engine (`CuboctWorldScene3D`) on 72-bit `VerletCell72` lattices:
- Solid, liquid, gas, light (`440..718 nm`), heat ($\hat{n}_E$), and shockwaves in a single branchless Verlet step;
- Live binaural stereo audio readout `(LeftEar, RightEar)` directly from lattice pressure;
- Zero-memory time reversal (`StepBackward`) for gameplay rewind and deterministic rollback netcode.
