// ===========================================================================
//  1D Spherical Sn Transport Code
// ===========================================================================
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include "Functions.hpp"
#include "InputParser.hpp"
#include "Solver.hpp"

namespace fs = std::filesystem;


// ---------------------------------------------------------------------------
//  Usage message
// ---------------------------------------------------------------------------
static void PrintUsage(const char* prog)
{
    std::cerr << "Usage: " << prog << " <input.yaml>\n";
}

// ---------------------------------------------------------------------------
//  Print a summary of the parsed InputData to stdout
// ---------------------------------------------------------------------------
static void PrintInputSummary(const InputData& d)
{
    std::cout << "\n== Input Summary =====================================\n";

    // Geometry
    std::cout << "Geometry: " << d.shells.size() << " shell(s)\n";
    for (const auto& s : d.shells)
    {
        std::cout << "  Shell " << s.id
                  << "  r_out=" << s.outer_radius << " cm"
                  << "  cells=" << s.num_cells
                  << "  mat=" << s.material
                  << "  src=" << s.source << "\n";
    }

    // Angular Quadrature table
    const auto& q = d.quadrature;
    std::cout << "Angular Quadrature: S" << q.order
              << "  (" << q.mu.size() << " directions)\n";
    std::cout << "  " << std::setw(3) << "m"
              << "  " << std::setw(20) << "mu_m"
              << "  " << std::setw(18) << "w_m" << "\n";
    std::cout << "  " << std::string(3+2+20+2+18, '-') << "\n";
    for (int m = 0; m < static_cast<int>(q.mu.size()); ++m)
    {
        std::cout << "  " << std::setw(3) << (m + 1)
                  << "  " << std::setw(20) << std::setprecision(12) << std::fixed << q.mu[m]
                  << "  " << std::setw(18) << std::setprecision(12) << std::fixed << q.w[m]
                  << "\n";
    }
    std::cout << std::defaultfloat;

    // Materials
    std::cout << "Materials: " << d.materials.size() << "\n";
    for (const auto& m : d.materials)
    {
        std::cout << "  [" << m.id << "]"
                  << "  sigma_t=" << m.sigma_t
                  << "  sigma_a=" << m.sigma_a
                  << "  P" << m.scattering_order;
        if (!m.sigma_s_moments.empty())
        {
            std::cout << "  sigma_s1=" << m.sigma_s_moments[0];
        }
        std::cout << "\n";
    }

    // Distributed sources
    std::cout << "Distributed sources: " << d.distributed_sources.size() << "\n";
    for (const auto& src : d.distributed_sources)
    {
        std::cout << "  [" << src.id << "]";
        if (!src.per_cell_strength.empty())
            std::cout << "  per-cell (" << src.per_cell_strength.size() << " values)";
        else
            std::cout << "  strength=" << src.strength;
        std::cout << "  normalize=" << (src.normalize ? "true" : "false") << "\n";
    }

    // Boundary
    std::cout << "Boundary: isotropic_flux=" << d.boundary.isotropic_flux
              << "  normalize=" << (d.boundary.normalize ? "true" : "false") << "\n";

    // Solver
    std::cout << "Solver: acceleration=" << d.solver.acceleration
              << "  tol=" << d.solver.convergence_tolerance
              << "  max_iter=" << d.solver.max_iterations << "\n";

    // Output
    std::cout << "Output: dir=" << d.output.directory << "\n";

    std::cout << "───────────────────────────────────────────────────────\n\n";
}

// ===========================================================================
int main(int argc, char* argv[])
{
    //  Command-line argument: path to YAML input file
    if (argc < 2)
    {
        PrintUsage(argv[0]);
        return EXIT_FAILURE;
    }
    const std::string input_path = argv[1];

    if (!fs::exists(input_path))
    {
        std::cerr << "Error: input file not found: " << input_path << "\n";
        return EXIT_FAILURE;
    }

    //  Parse the YAML input deck
    InputData input;
    try
    {
        InputParser parser(input_path);
        input = parser.parse();
    }
    catch (const std::exception& ex)
    {
        std::cerr << "Error: " << ex.what() << "\n";
        return EXIT_FAILURE;
    }

    PrintInputSummary(input);

    //  Run the transport solver
    try {
        SphericalSnSolver solver(input);
        solver.solve();
    }
    catch (const std::exception& ex) {
        std::cerr << "Solver error: " << ex.what() << "\n";
        return EXIT_FAILURE;
    }

    std::cout << "\nDone.\n";
    return EXIT_SUCCESS;
}
