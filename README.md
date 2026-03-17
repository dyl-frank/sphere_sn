# 1D Spherical Sn Transport Code

## Dependencies

- **Eigen 3.3+** - linear algebra
- **Axom** (Inlet + Sidre) - YAML input parsing
- **MFEM** - Gauss-Legendre quadrature

All three are managed via [Spack](https://spack.io).

## Building

The easiest way to build is with the provided script, which handles Spack setup, dependency installation, and CMake configuration automatically:

```bash
./build.sh
```

Common options:

| Flag | Effect |
|------|--------|
| `--clean` | Wipe build directory and rebuild from scratch |
| `--debug` | Build in Debug mode |
| `--no-deps` | Skip dependency checks (assume already loaded) |
| `--jobs N` | Override parallel job count |

To build manually:

```bash
spack load eigen axom mfem
mkdir build && cd build
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DEigen3_DIR=$(spack location -i eigen)/share/eigen3/cmake \
  -Daxom_DIR=$(spack location -i axom)/lib/cmake/axom \
  -DMFEM_DIR=$(spack location -i mfem)
make -j$(nproc)
```

## Running

```bash
./build/sn_transport examples/spherical_sn.yaml
```

## Visualization

Output `.ult` files can be visualized with [pydv](https://github.com/LLNL/pydv):

```bash
pydv scalar_flux.ult
```
