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

## Input File Reference

The solver is driven by a YAML input deck. The full structure is shown below; all fields with a default are optional.

```yaml
name: "my_problem"   # string label prepended to output file names

geometry:
  shells:            # list of concentric shells, innermost first
    - id: 1          # integer identifier (used for ordering)
      outer_radius:  # outer radius of this shell [cm]  (required)
      num_cells:     # number of uniform spatial cells  (required)
      material:      # material id string               (required)
      source:        # distributed source id, or "none" (default: "none")

angular_quadrature:
  sn_order:          # Sn order N — must be even        (required)

materials:
  - id:              # material identifier              (required)
    sigma_t:         # total cross section [1/cm]       (required)
    sigma_a:         # absorption cross section [1/cm]  (required)
    scattering_order: 0          # Pn scattering order K; 0 = isotropic (default: 0)
    sigma_s_moments: [s1, s2, …] # scattering moments σ_sk for k=1..K
                                 # (required when scattering_order > 0)

sources:

  # --- Volumetric distributed sources ---
  distributed:
    - id:            # source identifier                (required)
      strength:      # uniform source strength [p/cm³/s]
      per_cell_strength: [q1, q2, …]  # per-cell values; overrides 'strength'
      normalize: false  # if true, scale so ∫q dV = 1   (default: false)

  # --- Outer-boundary incident flux ---
  # Exactly one of scalar_flux / isotropic_flux / per_direction should be set.
  boundary:
    scalar_flux:      # isotropic scalar flux φ₀ [p/cm²/s]; solver sets ψ = φ₀/2
                      # (default: 0.0)
    isotropic_flux:   # isotropic angular flux ψ applied to all incoming directions
                      # (default: 0.0; use scalar_flux instead when φ is known)
    per_direction: [ψ₁, ψ₂, …]  # explicit ψ_m for each μ_m < 0, length = N/2
    normalize: false  # if true, scale so incoming half-range current J⁺ = 1
                      # (default: false)

solver:
  acceleration: none        # "none" or "dsa"                 (default: none)
  convergence_tolerance: 1.0e-4  # L∞ relative change in φ   (default: 1e-4)
  max_iterations: 500       # iteration cap                   (default: 500)

output:
  directory: ./results/          # output directory            (default: ./results/)
  scalar_flux_csv: true          # write φᵢ to CSV             (default: true)
  scalar_flux_pdv: true          # write φᵢ to pydv .ult file  (default: true)
  balance_table: true            # write global balance table  (default: true)
  angular_flux_boundary: true    # write ψ_m at outer boundary (default: true)
  starting_direction_origin: true # write ψ₁/₂ at r = 0       (default: true)
```

### Boundary condition field priority

When `normalize: false` the incoming angular flux for direction $m$ is resolved in this order:

1. `per_direction[m]` — if the list is provided
2. `isotropic_flux` — if non-zero (sets ψ = value for every incoming direction)
3. `scalar_flux / 2` — converts an isotropic scalar flux to an angular flux

When `normalize: true` the raw ψ values above are scaled so that the incoming half-range current equals 1:

$$J^+ = \sum_{\mu_m < 0} (-\mu_m)\, w_m\, \psi_m = 1$$

### Multi-shell example

```yaml
geometry:
  shells:
    - id: 1
      outer_radius: 5.0
      num_cells: 50
      material: fuel
      source: src_fuel
    - id: 2
      outer_radius: 6.0
      num_cells: 10
      material: reflector
      source: none
```

Shells are sorted by `id` internally, so order in the file does not matter.

## Visualization

Output `.ult` files can be visualized with [pydv](https://github.com/LLNL/pydv):

```bash
pydv scalar_flux.ult
```
