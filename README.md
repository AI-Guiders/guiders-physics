# AIGuiders.Physics (`guiders-physics`)

Native `.NET 11` + `AVX2/AVX-512` + `NVIDIA CUDA Driver API (PTX)` packages based on the discrete 12-neighbor Cuboctahedron $O_2(3)$ carrier calculus (`Ring9` $\in \mathbb{Z}_{512}$, 72-bit two-tick `VerletCell72`).

## NuGet Packages
- **`AIGuiders.Physics.Core`** (`Guiders.Physics.Core`): Native 9-bit ring `Ring9`, 72-bit cell `VerletCell72` (`9 bytes = 72 bits`, zero padding), branchless Gray-quadrant butterfly rotator, 12-neighbor cuboctahedron holonomy gate, multimodal sensor readout (optical wavelength `440..718 nm`, RGB triad `R+G+B=128`, acoustic pressure, thermal quanta), and direct `nvcuda.dll` PTX runner.
- **`AIGuiders.Physics.Compression`** (`Guiders.Physics.Compression`): Universal raw-byte stream & file archiver (like `GZipStream`) with automatic Heisenberg `[-40, +40]` deadband probing and reversible 1st-order / 2-point Verlet predictors.
- **`AIGuiders.Physics.MediaCodec`** (`Guiders.Physics.MediaCodec`): 4D lossless & two-point keyframe-drop (`Stride N = 2, 4, 8, 12`) video/audio streaming codec using SIMD `Vector256<byte>` and 2nd-order Verlet parabola prediction.
- **`AIGuiders.Physics.Engine`** (`Guiders.Physics.Engine`): Unified 3D physical, optical, thermal, fracture, and binaural acoustic world simulation engine on 72-bit `VerletCell72` lattices with zero-memory time reversal.
