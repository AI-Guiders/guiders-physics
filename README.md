# AIGuiders.Physics (`guiders-physics`)

Native `.NET 11` + `AVX2/AVX-512` + `NVIDIA CUDA Driver API (PTX)` packages based on the discrete 12-neighbor Cuboctahedron $O_2(3)$ carrier calculus (`Ring9` $\in \mathbb{Z}_{512}$, 72-bit two-tick `VerletCell72`).

## NuGet Packages
- **`AIGuiders.Physics.Core`** (`Guiders.Physics.Core`): Native 9-bit ring `Ring9`, 72-bit cell `VerletCell72` (`9 bytes = 72 bits`, zero padding), branchless Gray-quadrant butterfly rotator, 12-neighbor cuboctahedron holonomy gate, multimodal sensor readout (optical wavelength `440..718 nm`, RGB triad `R+G+B=128`, acoustic pressure, thermal quanta), and direct `nvcuda.dll` PTX runner.
- **`AIGuiders.Physics.Compression`** (`Guiders.Physics.Compression`): Universal raw-byte stream & file archiver (like `GZipStream`) with automatic Heisenberg `[-40, +40]` deadband probing and reversible 1st-order / 2-point Verlet predictors.
- **`AIGuiders.Physics.MediaCodec`** (`Guiders.Physics.MediaCodec`): 4D lossless & two-point keyframe-drop (`Stride N = 2, 4, 8, 12`) video/audio streaming codec using SIMD `Vector256<byte>` and 2nd-order Verlet parabola prediction.
- **`AIGuiders.Physics.Engine`** (`Guiders.Physics.Engine`): Unified 3D physical, optical, thermal, fracture, and binaural acoustic world simulation engine on 72-bit `VerletCell72` lattices with zero-memory time reversal.

---

## Open Science Invitation, Community Verification & Honest Disclaimer

**Science and engineering do not belong to isolated individuals — they belong to humanity and to all thinking actors on this planet.**

- **Working Pure CA Core (`Guiders.Physics.Core`):** The 5-file pure cellular automaton core (`Ring9.cs`, `Bus36.cs`, `VerletCell72.cs`, `CuboctahedronCarrier.cs`, `CuboctahedronAutomatonKernel.cs`) is branchless, bit-exact, and verified by unit tests (`17/17`).
- **Honest Disclaimer:** Higher-level engine adapters, SIMD/PTX kernels, and documentation across the broader ecosystem **may contain bugs, typos, or edge cases that require independent testing, benchmarking, and refinement**. One person working with AI agents cannot exhaustively verify every hardware target and domain layer alone.
- **Open Like Wikipedia:** All code is open-source. We invite `.NET`, GPU/CUDA, and physics simulation engineers, researchers, and AI actors worldwide to inspect the code, run benchmarks, report issues, and contribute Pull Requests.

