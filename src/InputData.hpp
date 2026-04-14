// ===========================================================================
//  InputData.hpp — Plain data structures for the 1D Spherical Sn input deck
// ===========================================================================
#pragma once

#include "Functions.hpp"

// ---------------------------------------------------------------------------
//  Geometry
// ---------------------------------------------------------------------------
struct ShellInput
{
    int         id           = 0;
    double      outer_radius = 0.0;   // [cm]
    int         num_cells    = 0;
    std::string material;             // references a MaterialInput::id
    std::string source;               // references a DistributedSourceInput::id, or "none"
};

// ---------------------------------------------------------------------------
//  Materials
// ---------------------------------------------------------------------------
struct MaterialInput
{
    std::string          id;
    double               sigma_t          = 0.0;
    double               sigma_a          = 0.0;
    int                  scattering_order = 0;
    std::vector<double>  sigma_s_moments;   // length = scattering_order (k=1..K)
};

// ---------------------------------------------------------------------------
//  Sources
// ---------------------------------------------------------------------------
struct DistributedSourceInput
{
    std::string          id;
    double               strength = 0.0;     // uniform q0 [p/cm³/s]
    std::vector<double>  per_cell_strength;  // per-cell override; takes precedence if non-empty
    bool                 normalize = false;
};

struct BoundaryInput
{
    double               isotropic_flux   = 0.0;  // angular flux ψ, applied isotropically
    double               scalar_flux      = 0.0;  // isotropic scalar flux φ; solver sets ψ = φ/2
    std::vector<double>  per_direction;            // explicit ψ_m for µm < 0, length = N/2
    bool                 normalize = false;
};

// ---------------------------------------------------------------------------
//  Solver
// ---------------------------------------------------------------------------
struct SolverInput
{
    std::string acceleration          = "none";   // "none" | "dsa"
    double      convergence_tolerance = 1.0e-4;
    int         max_iterations        = 500;
};

// ---------------------------------------------------------------------------
//  Angular Quadrature
// ---------------------------------------------------------------------------
struct AngularQuadrature
{
    int                 order = 4;
    std::vector<double> mu;   // direction cosines μ_m, sorted ascending ∈ [-1,1]
    std::vector<double> w;    // quadrature weights, Σ w_m = 2
};

// ---------------------------------------------------------------------------
//  Output
// ---------------------------------------------------------------------------
struct OutputInput
{
    std::string directory                  = "./results/";
    bool        scalar_flux_csv            = true;
    bool        scalar_flux_pdv            = true;
    bool        balance_table              = true;
    bool        angular_flux_boundary      = true;
    bool        starting_direction_origin  = true;
};

// ---------------------------------------------------------------------------
//  Top-level container
// ---------------------------------------------------------------------------
struct InputData
{
    // name
    std::string                         name;
    // geometry
    std::vector<ShellInput>             shells;
    // angular quadrature
    AngularQuadrature                   quadrature;
    // materials
    std::vector<MaterialInput>          materials;
    // sources
    std::vector<DistributedSourceInput> distributed_sources;
    BoundaryInput                       boundary;
    // solver
    SolverInput                         solver;
    // output
    OutputInput                         output;
};
