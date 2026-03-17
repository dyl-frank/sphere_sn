// ===========================================================================
//  Solver.hpp — 1D Spherical-Geometry Sn Transport Solver
//
//  Primary reference: AllNotes.pdf, Lecture 25
//    "Discretization and Solution of the Spherical-Geometry Sn Equations"
//
//  Project assignment sections addressed (project_2026.pdf):
//    Task 1 – Build the spatial mesh from shell definitions
//    Task 2 – Compute the angular quadrature grid (\alpha-coefficients, \beta weights)
//    Task 3 – Transport sweep with diamond-difference (DD) / weighted-diamond
//              (WD) and starting-direction treatment
//    Task 4 – Source iteration loop (convergence in L_\inf scalar-flux norm)
//    Task 5 – DSA acceleration (optional; activated by solver.acceleration=dsa)
//    Task 6 – Output: scalar flux, balance table, boundary angular flux,
//              starting-direction flux at origin
//
//  Equation numbers below refer to Lecture 25 unless noted otherwise.
//  Linear algebra uses mfem::Vector and mfem::DenseMatrix (no FEM required).
// ===========================================================================
#pragma once

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

#include <mfem.hpp>

#include "Functions.hpp"
#include "InputData.hpp"

namespace fs = std::filesystem;

// ---------------------------------------------------------------------------
//  Helper: Thomas algorithm for tridiagonal system  a*x_{i-1} + b*x_i + c*x_{i+1} = d
//  Used by apply_dsa() to solve the 1D spherical diffusion correction equation.
// ---------------------------------------------------------------------------
static mfem::Vector thomas(const mfem::Vector& a, const mfem::Vector& b,
                            const mfem::Vector& c, const mfem::Vector& d)
{
    int n = b.Size();
    mfem::Vector cp(n), dp(n), x(n);
    cp(0) = c(0) / b(0);
    dp(0) = d(0) / b(0);
    for (int i = 1; i < n; ++i) {
        double denom = b(i) - a(i) * cp(i-1);
        cp(i) = (i < n-1) ? c(i) / denom : 0.0;
        dp(i) = (d(i) - a(i) * dp(i-1)) / denom;
    }
    x(n-1) = dp(n-1);
    for (int i = n-2; i >= 0; --i)
        x(i) = dp(i) - cp(i) * x(i+1);
    return x;
}

// ===========================================================================
class SphericalSnSolver
{
public:
    explicit SphericalSnSolver(const InputData& inp) : msh_inp(inp) {}

    void solve()
    {
        using Clock = std::chrono::steady_clock;
        using Ms    = std::chrono::duration<double, std::milli>;

        auto t0 = Clock::now();
        build_mesh();
        auto t1 = Clock::now();

        build_material_arrays();
        auto t2 = Clock::now();

        build_source_arrays();
        auto t3 = Clock::now();

        compute_angular_grid();
        compute_bc_scale();
        auto t4 = Clock::now();

        source_iteration();
        auto t5 = Clock::now();

        write_output();
        auto t6 = Clock::now();

        double total = Ms(t6 - t0).count();
        std::cout << "\n--- Runtime breakdown ---\n"
                  << std::fixed << std::setprecision(3)
                  << "  build_mesh          : " << std::setw(10) << Ms(t1-t0).count() << " ms\n"
                  << "  build_material_arrays: " << std::setw(9) << Ms(t2-t1).count() << " ms\n"
                  << "  build_source_arrays : " << std::setw(10) << Ms(t3-t2).count() << " ms\n"
                  << "  compute_angular_grid: " << std::setw(10) << Ms(t4-t3).count() << " ms\n"
                  << "  source_iteration    : " << std::setw(10) << Ms(t5-t4).count() << " ms\n"
                  << "  write_output        : " << std::setw(10) << Ms(t6-t5).count() << " ms\n"
                  << "  TOTAL               : " << std::setw(10) << total              << " ms\n";
    }

private:
    const InputData& msh_inp;

    //  mesh 
    int              msh_I = 0;
    mfem::Vector     msh_r_edge;      // r_{i+1/2}, size I+1;  [0] = 0
    mfem::Vector     msh_r_ctr;       // cell centers,   size I
    mfem::Vector     msh_A;           // size I+1
    mfem::Vector     msh_V;           // size I
    mfem::Vector     msh_dA;          // A[i+1]−A[i],    size I
    mfem::Array<int> msh_shell_idx;   // which shell each cell belongs to, size I

    //  quadrature 
    // Lec 25, Eqs. 2, 9, 12
    int          msh_N_angle = 0;
    mfem::Vector msh_mu;        // direction cosines sorted ascending,  size N
    mfem::Vector msh_w_angle;   // weights (sum = 2),                   size N
    mfem::Vector msh_mu_half;   // angular cell edges \mu_{m+1/2},        size N+1
    mfem::Vector msh_alpha;     // \alpha_{m+1/2},                            size N+1
    mfem::Vector msh_beta;      // \beta_m,                                  size N

    //  per-cell cross sections 
    mfem::Vector      msh_sigma_t;    // \Sigma_t,      size I
    mfem::Vector      msh_sigma_a;    // \Sigma_a,      size I
    mfem::Vector      msh_sigma_s0;   // \Sigma_{s,0},  size I
    int               msh_K_max = 0;  // max scattering order across all materials
    mfem::DenseMatrix msh_sigma_sk;   // \Sigma_{s,k}(k,i), k=1..K_max; size (K_max x I)

    //  external source 
    mfem::Vector msh_q_dist;    // isotropic volumetric q_0 [p/cm^3/s], size I

    //  solution arrays 
    mfem::DenseMatrix msh_psi_cell;          // \psi_{i,m},  size (I x N)
    mfem::DenseMatrix msh_phi_l;             // \phi_k(i),   size (K_max+1 x I)
    mfem::Vector      msh_psi_start_cell;    // \psi_s cell averages, size I
    double            msh_psi_start_origin = 0.0;  // \psi_s at r = 0

    //  output bookkeeping
    mfem::Vector msh_psi_bnd;   // \psi_{I+1/2,m} at outer boundary, size N

    //  boundary normalization scale (1.0 unless bc.normalize = true)
    double msh_bc_scale = 1.0;

    // =======================================================================
    //  TASK 1 — BUILD SPATIAL MESH
    // =======================================================================
    void build_mesh()
    {
        // Pre-compute total number of cells
        msh_I = 0;
        for (const auto& shell : msh_inp.shells)
            msh_I += shell.num_cells;

        // Allocate all mesh arrays
        msh_r_edge.SetSize(msh_I + 1);
        msh_r_ctr.SetSize(msh_I);
        msh_A.SetSize(msh_I + 1);
        msh_V.SetSize(msh_I);
        msh_dA.SetSize(msh_I);
        msh_shell_idx.SetSize(msh_I);

        // Fill radial edges
        msh_r_edge(0) = 0.0;
        int cell = 0;
        double r_prev = 0.0;
        for (int s = 0; s < static_cast<int>(msh_inp.shells.size()); ++s) {
            const auto& shell = msh_inp.shells[s];
            double r_out = shell.outer_radius;
            int    n_cells = shell.num_cells;
            double h       = (r_out - r_prev) / n_cells;
            for (int j = 0; j < n_cells; ++j) {
                msh_r_edge(cell + j + 1) = r_prev + (j + 1) * h;
                msh_shell_idx[cell + j]  = s;
            }
            cell  += n_cells;
            r_prev = r_out;
        }

        for (int i = 0; i <= msh_I; ++i) {
            double r = msh_r_edge(i);
            msh_A(i) = 4.0 * M_PI * r * r;
        }

        // Cell volumes, centers, and dA (Lec 25, Eq. 19)
        for (int i = 0; i < msh_I; ++i) {
            double r_in  = msh_r_edge(i);
            double r_out = msh_r_edge(i + 1);
            msh_r_ctr(i) = 0.5 * (r_in + r_out);
            msh_V(i)     = (4.0 / 3.0) * M_PI * (r_out*r_out*r_out - r_in*r_in*r_in);
            msh_dA(i)    = msh_A(i + 1) - msh_A(i);
        }
    }

    // =======================================================================
    //  MATERIAL ARRAYS
    // =======================================================================
    void build_material_arrays()
    {
        std::map<std::string, int> mat_idx;
        for (int j = 0; j < static_cast<int>(msh_inp.materials.size()); ++j)
            mat_idx[msh_inp.materials[j].id] = j;

        for (const auto& mat : msh_inp.materials)
            msh_K_max = std::max(msh_K_max, mat.scattering_order);

        msh_sigma_t.SetSize(msh_I);
        msh_sigma_a.SetSize(msh_I);
        msh_sigma_s0.SetSize(msh_I);
        // Size (K_max x I); use at least 1 row to avoid zero-size DenseMatrix
        msh_sigma_sk.SetSize(msh_K_max > 0 ? msh_K_max : 1, msh_I);
        msh_sigma_sk = 0.0;

        for (int i = 0; i < msh_I; ++i) {
            int s = msh_shell_idx[i];
            const auto& mat = msh_inp.materials[mat_idx.at(msh_inp.shells[s].material)];
            msh_sigma_t(i)  = mat.sigma_t;
            msh_sigma_a(i)  = mat.sigma_a;
            msh_sigma_s0(i) = mat.sigma_t - mat.sigma_a;
            for (int k = 0; k < mat.scattering_order; ++k)
                msh_sigma_sk(k, i) = mat.sigma_s_moments[k];
        }
    }

    // =======================================================================
    //  SOURCE ARRAYS
    // =======================================================================
    void build_source_arrays()
    {
        std::map<std::string, int> src_idx;
        for (int j = 0; j < static_cast<int>(msh_inp.distributed_sources.size()); ++j)
            src_idx[msh_inp.distributed_sources[j].id] = j;

        msh_q_dist.SetSize(msh_I);
        msh_q_dist = 0.0;
        mfem::Array<int> shell_local_ctr(static_cast<int>(msh_inp.shells.size()));
        shell_local_ctr = 0;

        for (int i = 0; i < msh_I; ++i) {
            int s = msh_shell_idx[i];
            const std::string& src_id = msh_inp.shells[s].source;
            if (src_id == "none" || src_id.empty()) { shell_local_ctr[s]++; continue; }
            auto it = src_idx.find(src_id);
            if (it == src_idx.end())
                throw std::runtime_error("Source id not found: " + src_id);
            const auto& src = msh_inp.distributed_sources[it->second];
            int local = shell_local_ctr[s]++;
            msh_q_dist(i) = !src.per_cell_strength.empty()
                            ? src.per_cell_strength[local]
                            : src.strength;
        }

        // Optional normalization: scale so \int q dV = 1 within each shell
        for (int s = 0; s < static_cast<int>(msh_inp.shells.size()); ++s) {
            const std::string& src_id = msh_inp.shells[s].source;
            if (src_id == "none" || src_id.empty()) continue;
            auto it = src_idx.find(src_id);
            if (it == src_idx.end() || !msh_inp.distributed_sources[it->second].normalize) continue;
            double vol_integral = 0.0;
            for (int i = 0; i < msh_I; ++i)
                if (msh_shell_idx[i] == s) vol_integral += msh_q_dist(i) * msh_V(i);
            if (vol_integral > 0.0)
                for (int i = 0; i < msh_I; ++i)
                    if (msh_shell_idx[i] == s) msh_q_dist(i) /= vol_integral;
        }
    }

    // =======================================================================
    //  TASK 2 — ANGULAR QUADRATURE GRID
    // =======================================================================
    void compute_angular_grid()
    {
        msh_N_angle = msh_inp.quadrature.order;
        int N = msh_N_angle;

        // Copy direction cosines and weights (already sorted ascending)
        msh_mu.SetSize(N);
        msh_w_angle.SetSize(N);
        for (int m = 0; m < N; ++m) {
            msh_mu(m)      = msh_inp.quadrature.mu[m];
            msh_w_angle(m) = msh_inp.quadrature.w[m];
        }

        // Angular cell-edge cosines (Lec 25, Eq. 2):
        //   \mu_{-1/2} = -1,  \mu_{m+1/2} = \mu_{m-1/2} + w_m
        msh_mu_half.SetSize(N + 1);
        msh_mu_half(0) = -1.0;
        for (int m = 0; m < N; ++m)
            msh_mu_half(m + 1) = msh_mu_half(m) + msh_w_angle(m);

        // Alpha-coefficients (Lec 25, Eqs. 8-9):
        //   \alpha_{-1/2} = 0,  \alpha_{m+1/2} = \alpha_{m-1/2} − 2 \mu_m w_m
        msh_alpha.SetSize(N + 1);
        msh_alpha(0) = 0.0;
        for (int m = 0; m < N; ++m)
            msh_alpha(m + 1) = msh_alpha(m) - 2.0 * msh_mu(m) * msh_w_angle(m);

        // Beta weights for weighted-diamond in angle (Lec 25, Eq. 12):
        //   \beta_m = (\mu_m − \mu_{m-1/2}) / (\mu_{m+1/2} − \mu_{m-1/2})
        msh_beta.SetSize(N);
        for (int m = 0; m < N; ++m) {
            double dmu  = msh_mu_half(m + 1) - msh_mu_half(m);
            msh_beta(m) = (msh_mu(m) - msh_mu_half(m)) / dmu;
        }
    }

    // =======================================================================
    //  TOTAL SOURCE ASSEMBLY
    // =======================================================================
    // Returns Q(m, i): angular source for direction m, cell i (Lec 25, RHS Eq. 1)
    mfem::DenseMatrix compute_total_source() const
    {
        mfem::DenseMatrix Q(msh_N_angle, msh_I);
        Q = 0.0;
        for (int m = 0; m < msh_N_angle; ++m) {
            for (int i = 0; i < msh_I; ++i) {
                double q = 0.5 * msh_q_dist(i);                     // isotropic distributed
                q += 0.5 * msh_sigma_s0(i) * msh_phi_l(0, i);      // scatter k=0
                for (int k = 1; k <= msh_K_max; ++k)                // scatter k=1..K
                    q += (2.0*k + 1) / 2.0 * msh_sigma_sk(k-1, i)
                         * msh_phi_l(k, i) * Pn(k, msh_mu(m));
                Q(m, i) = q;
            }
        }
        return Q;
    }

    // =======================================================================
    //  BOUNDARY CONDITION HELPERS
    // =======================================================================

    // Compute msh_bc_scale so that the incoming half-range current = 1 when
    // bc.normalize = true.  Called once after compute_angular_grid().
    //
    //   J_in = sum_{mu_m < 0} (-mu_m) * w_m * psi_m
    //
    // For isotropic BC:  J_in = psi_in * sum_{mu_m < 0} (-mu_m) * w_m
    // For per_direction: J_in = sum_{m < N/2}  (-mu_m) * w_m * psi_m
    void compute_bc_scale()
    {
        if (!msh_inp.boundary.normalize) { msh_bc_scale = 1.0; return; }

        const auto& bc = msh_inp.boundary;
        double J_in = 0.0;
        int N_in = msh_N_angle / 2;   // negative-mu directions are m = 0..N/2-1
        for (int m = 0; m < N_in; ++m) {
            double psi_raw = bc.per_direction.empty() ? bc.isotropic_flux
                                                      : (m < static_cast<int>(bc.per_direction.size())
                                                         ? bc.per_direction[m] : 0.0);
            J_in += (-msh_mu(m)) * msh_w_angle(m) * std::max(psi_raw, 0.0);
        }
        if (J_in <= 0.0)
            throw std::runtime_error("bc.normalize=true but incoming partial current is zero");
        msh_bc_scale = 1.0 / J_in;
        std::cout << "  BC normalize: raw J_in = " << J_in
                  << "  -> psi scale = " << msh_bc_scale << "\n";
    }

    // Angular flux at the outer boundary for inward direction m (mu_m < 0).
    double bc_value(int m) const
    {
        const auto& bc = msh_inp.boundary;
        double val = 0.0;
        if (!bc.per_direction.empty()) {
            if (m < static_cast<int>(bc.per_direction.size()))
                val = bc.per_direction[m];
        } else {
            val = bc.isotropic_flux;
        }
        return std::max(val, 0.0) * msh_bc_scale;
    }

    // BC at the outer boundary for the starting direction (mu = -1).
    // Uses linear interpolation from the two most-negative quadrature directions
    // per project item 7d (reduces to psi_in for isotropic BC).
    double bc_value_start() const
    {
        const auto& bc = msh_inp.boundary;
        double psi1 = bc.per_direction.empty() ? bc.isotropic_flux
                      : (static_cast<int>(bc.per_direction.size()) > 0 ? bc.per_direction[0] : 0.0);
        double psi2 = bc.per_direction.empty() ? bc.isotropic_flux
                      : (static_cast<int>(bc.per_direction.size()) > 1 ? bc.per_direction[1] : psi1);
        psi1 = std::max(psi1, 0.0) * msh_bc_scale;
        psi2 = std::max(psi2, 0.0) * msh_bc_scale;

        if (msh_N_angle == 2) return psi1;   // S2: constant extrapolation

        // Linear interpolation to mu = -1 (= mu_{1/2})
        double mu1 = msh_mu(0), mu2 = msh_mu(1);
        double psi_start = psi1 * ((-1.0) - mu2) / (mu1 - mu2)
                         + psi2 * ((-1.0) - mu1) / (mu2 - mu1);
        return std::max(psi_start, 0.0);
    }

    // =======================================================================
    //  STARTING-DIRECTION SWEEP  (slab geometry at \mu = −1)
    // =======================================================================
    void starting_direction_sweep(const mfem::DenseMatrix& Q)
    {
        msh_psi_start_cell.SetSize(msh_I);
        msh_psi_start_cell = 0.0;

        // Source interpolated to \mu=-1 from m=0 and m=1 (Lec. 25, Eq. 15)
        mfem::Vector Qs(msh_I);
        double mu0 = msh_mu(0), mu1 = msh_mu(1);
        double dmu = mu1 - mu0;
        for (int i = 0; i < msh_I; ++i)
            Qs(i) = Q(0, i) * (mu1 + 1.0) / dmu - Q(1, i) * (mu0 + 1.0) / dmu;

        // Diamond-difference slab sweep inward at \mu=-1 (Lec. 25, Eq. 14):
        //   (\psi_out−\psi_in)/h + \Sigma_t*(\psi_in+\psi_out)/2 = Q_s
        //    \psi_out = (Q_s*h + (1−coeff)*\psi_in) / (1+coeff),  coeff = \Sigma_t*h/2
        double psi_in = bc_value_start();
        for (int i = msh_I - 1; i >= 0; --i) {
            double h       = msh_r_edge(i + 1) - msh_r_edge(i);
            double coeff   = 0.5 * msh_sigma_t(i) * h;
            double psi_out = (Qs(i) * h + (1.0 - coeff) * psi_in) / (1.0 + coeff);
            psi_out = std::max(psi_out, 0.0);
            msh_psi_start_cell(i) = 0.5 * (psi_in + psi_out);
            psi_in = psi_out;
        }
        msh_psi_start_origin = psi_in;
    }

    // =======================================================================
    //  TASK 3 — TRANSPORT SWEEP  (all angular and spatial cells)
    //
    //  Reference: Lec 25, Sec 3 "Spatial Discretization" + Sec 4 "Source Iteration"
    //
    //  Angular redistribution coefficients (common to both sweep directions):
    //    alpha_right = alpha_{m+1/2} = msh_alpha(m+1)
    //    alpha_left  = alpha_{m-1/2} = msh_alpha(m)
    //    ang_redist_diag   = dA[i] * alpha_right / (2 * quad_weight * wd_beta)
    //    ang_redist_source = (dA[i] / (2*quad_weight)) * (alpha_right*(1-wd_beta)/wd_beta + alpha_left)
    //
    //  NEGATIVE DIRECTIONS  m = 0..N/2-1  (dir_cosine < 0, inward sweep i = I-1..0):
    //    lhs_coeff = -2*dir_cosine*face_area_inner + ang_redist_diag + Sigma_t*cell_volume
    //    rhs       = Q*cell_volume - dir_cosine*(A[i+1]+A[i])*psi_sp_inflow
    //                + ang_redist_source*psi_ang_inflow
    //    DD:  psi_sp[i]   = 2*psi_cell_avg - psi_sp_inflow
    //    WD:  psi_ang_out = (psi_cell_avg - (1-wd_beta)*psi_ang_inflow) / wd_beta
    //
    //  POSITIVE DIRECTIONS  m = N/2..N-1  (dir_cosine > 0, outward sweep i = 0..I-1):
    //    psi_sp[0] = msh_psi_start_origin  (origin BC for all outward dirs)
    //    lhs_coeff = 2*dir_cosine*face_area_outer + ang_redist_diag + Sigma_t*cell_volume
    //    rhs       = Q*cell_volume + dir_cosine*(A[i+1]+A[i])*psi_sp_inflow
    //                + ang_redist_source*psi_ang_inflow
    //    DD:  psi_sp[i+1] = 2*psi_cell_avg - psi_sp_inflow
    //    WD:  psi_ang_out = (psi_cell_avg - (1-wd_beta)*psi_ang_inflow) / wd_beta
    // =======================================================================
    void transport_sweep(const mfem::DenseMatrix& Q)
    {
        starting_direction_sweep(Q);

        msh_psi_cell.SetSize(msh_I, msh_N_angle);
        msh_psi_cell = 0.0;
        msh_psi_bnd.SetSize(msh_N_angle);
        msh_psi_bnd = 0.0;

        // psi_ang_in[i]: angular inflow for current direction m
        // Initialized to starting-direction cell averages (\psî_{i,-1/2})
        mfem::Vector psi_ang_in(msh_I);
        for (int i = 0; i < msh_I; ++i)
            psi_ang_in(i) = msh_psi_start_cell(i);

        mfem::Vector psi_sp(msh_I + 1);   // spatial edge fluxes
        mfem::Vector psi_ang_out(msh_I);   // WD outflow  becomes inflow for m+1

        for (int m = 0; m < msh_N_angle; ++m) {
            double dir_cosine  = msh_mu(m);
            double quad_weight = msh_w_angle(m);
            double wd_beta     = msh_beta(m);         // WD angular interpolation weight
            double alpha_right = msh_alpha(m + 1);    // alpha_{m+1/2}
            double alpha_left  = msh_alpha(m);        // alpha_{m-1/2}

            psi_sp = 0.0;

            if (dir_cosine < 0.0) {
                //  Inward sweep: i = I-1 downto 0
                psi_sp(msh_I) = bc_value(m);

                for (int i = msh_I - 1; i >= 0; --i) {
                    double psi_ang_inflow  = psi_ang_in(i);     // angular inflow (from m-1)
                    double psi_sp_inflow   = psi_sp(i + 1);     // spatial inflow (outer face)
                    double cell_volume     = msh_V(i);
                    double face_area_diff  = msh_dA(i);          // A_{i+1/2} - A_{i-1/2}
                    double face_area_inner = msh_A(i);           // outflow face (inner) for mu<0

                    // Angular redistribution coefficients (Lec 25, Eq. before 20)
                    double ang_redist_diag   = face_area_diff * alpha_right / (2.0 * quad_weight * wd_beta);
                    double ang_redist_source = (face_area_diff / (2.0 * quad_weight))
                                              * (alpha_right * (1.0 - wd_beta) / wd_beta + alpha_left);

                    double lhs_coeff = -2.0 * dir_cosine * face_area_inner
                                       + ang_redist_diag + msh_sigma_t(i) * cell_volume;
                    double rhs = Q(m, i) * cell_volume
                                 - dir_cosine * (msh_A(i+1) + face_area_inner) * psi_sp_inflow
                                 + ang_redist_source * psi_ang_inflow;
                    double psi_cell_avg = std::max(rhs / lhs_coeff, 0.0);

                    msh_psi_cell(i, m) = psi_cell_avg;
                    psi_sp(i) = std::max(2.0 * psi_cell_avg - psi_sp_inflow, 0.0);  // DD closure

                    if (m < msh_N_angle - 1)
                        psi_ang_out(i) = std::max(
                            (psi_cell_avg - (1.0 - wd_beta) * psi_ang_inflow) / wd_beta, 0.0);
                }
                msh_psi_bnd(m) = bc_value(m);

            } else {
                //  Outward sweep: i = 0 to I-1
                psi_sp(0) = msh_psi_start_origin;  // origin BC for all positive dirs

                for (int i = 0; i < msh_I; ++i) {
                    double psi_ang_inflow  = psi_ang_in(i);     // angular inflow (from m-1)
                    double psi_sp_inflow   = psi_sp(i);         // spatial inflow (inner face)
                    double cell_volume     = msh_V(i);
                    double face_area_diff  = msh_dA(i);          // A_{i+1/2} - A_{i-1/2}
                    double face_area_outer = msh_A(i + 1);       // outflow face (outer) for mu>0

                    // Angular redistribution coefficients (Lec 25, Eq. before 20)
                    double ang_redist_diag   = face_area_diff * alpha_right / (2.0 * quad_weight * wd_beta);
                    double ang_redist_source = (face_area_diff / (2.0 * quad_weight))
                                              * (alpha_right * (1.0 - wd_beta) / wd_beta + alpha_left);

                    double lhs_coeff = 2.0 * dir_cosine * face_area_outer
                                       + ang_redist_diag + msh_sigma_t(i) * cell_volume;
                    double rhs = Q(m, i) * cell_volume
                                 + dir_cosine * (face_area_outer + msh_A(i)) * psi_sp_inflow
                                 + ang_redist_source * psi_ang_inflow;
                    double psi_cell_avg = std::max(rhs / lhs_coeff, 0.0);

                    msh_psi_cell(i, m) = psi_cell_avg;
                    psi_sp(i + 1) = std::max(2.0 * psi_cell_avg - psi_sp_inflow, 0.0);  // DD closure

                    if (m < msh_N_angle - 1)
                        psi_ang_out(i) = std::max(
                            (psi_cell_avg - (1.0 - wd_beta) * psi_ang_inflow) / wd_beta, 0.0);
                }
                msh_psi_bnd(m) = psi_sp(msh_I);
            }

            // Pass WD outflow as angular inflow for next direction
            if (m < msh_N_angle - 1)
                psi_ang_in = psi_ang_out;
        }
    }

    // =======================================================================
    //  SCALAR FLUX LEGENDRE MOMENTS
    //    \phi_k(i) = \Sigma_m  w_m * P_k(\mu_m) * \psi_{i,m}
    //  Weights sum to 2, so \phi_0(i) = \int_{-1}^{1} \psi_i d\mu.
    // =======================================================================
    void update_phi_moments()
    {
        msh_phi_l.SetSize(msh_K_max + 1, msh_I);
        msh_phi_l = 0.0;
        for (int k = 0; k <= msh_K_max; ++k) {
            for (int i = 0; i < msh_I; ++i) {
                double sum = 0.0;
                for (int m = 0; m < msh_N_angle; ++m)
                    sum += msh_w_angle(m) * Pn(k, msh_mu(m)) * msh_psi_cell(i, m);
                msh_phi_l(k, i) = sum;
            }
        }
    }

    // =======================================================================
    //  TASK 5 — DSA CORRECTION
    //
    //  Lec 19, Sec 2, Eq. 23 (isotropic) / Eq. 38 (anisotropic, phi-only).
    //
    //  After the transport sweep produces phi^{l+1/2}, solve:
    //    -nabla.(D nabla Delta_phi) + Sigma_a Delta_phi
    //        = Sigma_{s,0} (phi^{l+1/2} - phi^l)             [Eq. 23c / 38c]
    //  where the transport-corrected diffusion coefficient is:
    //    D = 1 / (3 (Sigma_t - Sigma_{s,1}))                 [Eq. 38c]
    //  (reduces to 1/(3 Sigma_t) for isotropic scattering where Sigma_{s,1}=0)
    //  Then update (Eq. 38d): phi^{l+1} = phi^{l+1/2} + Delta_phi
    //
    //  BCs (Lec 19, Sec 2.1): vacuum transport BC -> vacuum diffusion BC at r=R;
    //  reflecting (spherical symmetry) at r=0.
    // =======================================================================
    void apply_dsa(const mfem::Vector& phi_old)
    {
        mfem::Vector eps(msh_I), D(msh_I);
        for (int i = 0; i < msh_I; ++i) {
            eps(i) = msh_phi_l(0, i) - phi_old(i);
            // Transport-corrected D (Lec 19 Eq. 38c): D = 1/(3(Sigma_t - Sigma_{s,1}))
            // msh_sigma_sk row 0 stores Sigma_{s,1}; zero when purely isotropic.
            double sigma_s1 = (msh_K_max >= 1) ? msh_sigma_sk(0, i) : 0.0;
            D(i) = 1.0 / (3.0 * (msh_sigma_t(i) - sigma_s1));
        }

        mfem::Vector a(msh_I), b(msh_I), c(msh_I), d(msh_I);
        a = 0.0;  b = 0.0;  c = 0.0;  d = 0.0;

        for (int i = 0; i < msh_I; ++i) {
            d(i) = msh_sigma_s0(i) * eps(i) * msh_V(i);

            if (i == 0) {
                // Reflecting BC at origin: no left face
                double h_c_r    = msh_r_ctr(1) - msh_r_ctr(0);
                double D_face_r = 2.0 * D(0) * D(1) / (D(0) + D(1));
                double fac_r    = msh_A(1) * D_face_r / h_c_r;
                b(0) = msh_sigma_a(0) * msh_V(0) + fac_r;
                a(0) = 0.0;
                c(0) = -fac_r;

            } else if (i == msh_I - 1) {
                // Vacuum BC at r=R: half-cell Marshak approximation
                double h_c_l    = msh_r_ctr(i) - msh_r_ctr(i - 1);
                double D_face_l = 2.0 * D(i) * D(i-1) / (D(i) + D(i-1));
                double fac_l    = msh_A(i) * D_face_l / h_c_l;
                double h_last   = msh_r_edge(msh_I) - msh_r_edge(msh_I - 1);
                double fac_r    = msh_A(msh_I) * 2.0 * D(i) / h_last;
                b(i) = msh_sigma_a(i) * msh_V(i) + fac_l + fac_r;
                a(i) = -fac_l;
                c(i) = 0.0;

            } else {
                // Interior cell: harmonic-mean face diffusion coefficients
                double h_c_r    = msh_r_ctr(i + 1) - msh_r_ctr(i);
                double h_c_l    = msh_r_ctr(i) - msh_r_ctr(i - 1);
                double D_face_r = 2.0 * D(i) * D(i+1) / (D(i) + D(i+1));
                double D_face_l = 2.0 * D(i) * D(i-1) / (D(i) + D(i-1));
                double fac_r    = msh_A(i + 1) * D_face_r / h_c_r;
                double fac_l    = msh_A(i)     * D_face_l / h_c_l;
                b(i) = msh_sigma_a(i) * msh_V(i) + fac_l + fac_r;
                a(i) = -fac_l;
                c(i) = -fac_r;
            }
        }

        mfem::Vector dphi = thomas(a, b, c, d);
        for (int i = 0; i < msh_I; ++i)
            msh_phi_l(0, i) += dphi(i);
    }

    // =======================================================================
    //  TASK 4 — SOURCE ITERATION LOOP
    // =======================================================================
    void source_iteration()
    {
        // Initialize \phi_l to zero
        msh_phi_l.SetSize(msh_K_max + 1, msh_I);
        msh_phi_l = 0.0;

        bool   use_dsa = (msh_inp.solver.acceleration == "dsa");
        double tol     = msh_inp.solver.convergence_tolerance;
        int    maxIt   = msh_inp.solver.max_iterations;

        using Clock = std::chrono::steady_clock;
        using Ms    = std::chrono::duration<double, std::milli>;

        double rel_err = 1.0 + tol;
        int    iter    = 0;
        double t_source = 0, t_sweep = 0, t_phi = 0, t_dsa = 0;

        while (iter < maxIt && rel_err > tol) {
            // Save current scalar flux for convergence check
            mfem::Vector phi_old(msh_I);
            for (int i = 0; i < msh_I; ++i)
                phi_old(i) = msh_phi_l(0, i);

            auto ta = Clock::now();
            mfem::DenseMatrix Q = compute_total_source();
            auto tb = Clock::now();
            transport_sweep(Q);
            auto tc = Clock::now();
            update_phi_moments();
            auto td = Clock::now();

            if (use_dsa)
                apply_dsa(phi_old);
            auto te = Clock::now();

            t_source += Ms(tb - ta).count();
            t_sweep  += Ms(tc - tb).count();
            t_phi    += Ms(td - tc).count();
            t_dsa    += Ms(te - td).count();

            // L_\inf relative convergence check
            rel_err = 0.0;
            for (int i = 0; i < msh_I; ++i) {
                double num = std::abs(msh_phi_l(0, i) - phi_old(i));
                double den = std::max(msh_phi_l(0, i), 1e-14);
                rel_err = std::max(rel_err, num / den);
            }

            ++iter;
            if (iter % 10 == 0)
                std::cout << "  iter " << std::setw(4) << iter
                          << "  rel_err = " << std::scientific
                          << std::setprecision(4) << rel_err << "\n";
        }

        if (rel_err <= tol)
            std::cout << "Converged in " << iter
                      << " iterations (rel_err=" << rel_err << ")\n";
        else
            std::cout << "NOT converged after " << iter
                      << " iterations (rel_err=" << rel_err << ")\n";

        std::cout << "  -- iteration timing (total across " << iter << " iters) --\n"
                  << std::fixed << std::setprecision(3)
                  << "    compute_total_source: " << std::setw(10) << t_source << " ms\n"
                  << "    transport_sweep     : " << std::setw(10) << t_sweep  << " ms\n"
                  << "    update_phi_moments  : " << std::setw(10) << t_phi    << " ms\n";
        if (use_dsa)
            std::cout << "    apply_dsa           : " << std::setw(10) << t_dsa << " ms\n";
    }

    // =======================================================================
    //  TASK 6 — OUTPUT
    // =======================================================================
    void write_output()
    {
        fs::path dir = msh_inp.output.directory;
        fs::create_directories(dir);
        std::string d = dir.string();

        if (msh_inp.output.scalar_flux_csv)               write_scalar_flux_csv(d);
        if (msh_inp.output.scalar_flux_pdv)               write_scalar_flux_pdv(d);
        if (msh_inp.output.balance_table)             write_balance_table(d);
        if (msh_inp.output.angular_flux_boundary)     write_angular_flux_boundary(d);
        if (msh_inp.output.starting_direction_origin) write_starting_direction(d);
    }

    //  scalar_flux.csv 
    // Columns: cell, r_center_cm, phi_0 [, phi_1, ..., phi_K]
    void write_scalar_flux_csv(const std::string& dir)
    {
        std::ofstream f(dir + "/" + msh_inp.name + "scalar_flux.csv");
        f << std::scientific << std::setprecision(8);
        f << "cell,r_center_cm,phi_0";
        for (int k = 1; k <= msh_K_max; ++k) f << ",phi_" << k;
        f << "\n";
        for (int i = 0; i < msh_I; ++i) {
            f << i << "," << msh_r_ctr(i);
            for (int k = 0; k <= msh_K_max; ++k) f << "," << msh_phi_l(k, i);
            f << "\n";
        }
    }

    void write_scalar_flux_pdv(const std::string& dir)
    {
        std::ofstream f(dir + "/" + msh_inp.name + "scalar_flux.ult");
        f << std::scientific << std::setprecision(8);
        for (int k = 0; k <= msh_K_max; ++k)
        {
            f << "# phi-" << k << " vs. cell center" << "\n";
            for (int i = 0; i < msh_I; ++i)
            {
                f << msh_r_ctr(i) << " "  << msh_phi_l(k, i);
                f << "\n";
            }
            f << "\n\n";
        }
    }

    //  balance_table.txt 
    // source_rate = \Sigma_i q[i]*V[i]
    // abs_rate    = \Sigma_i \Sigma_a[i]*\phi_0[i]*V[i]
    // incident    = A[I] * \Sigma_{m<N/2} (−\mu_m)*w_m*bc_value(m)
    // leakage     = A[I] * \Sigma_m \mu_m*w_m*\psi_bnd[m]
    void write_balance_table(const std::string& dir)
    {
        double source_rate = 0.0, abs_rate = 0.0;
        for (int i = 0; i < msh_I; ++i) {
            source_rate += msh_q_dist(i) * msh_V(i);
            abs_rate    += msh_sigma_a(i) * msh_phi_l(0, i) * msh_V(i);
        }

        double incident = 0.0;
        for (int m = 0; m < msh_N_angle / 2; ++m)
            incident += (-msh_mu(m)) * msh_w_angle(m) * bc_value(m);
        incident *= msh_A(msh_I);

        // Outgoing (leakage): only positive-mu directions at outer boundary.
        // Correct balance: source_rate + incident = abs_rate + leakage
        double leakage = 0.0;
        for (int m = msh_N_angle / 2; m < msh_N_angle; ++m)
            leakage += msh_mu(m) * msh_w_angle(m) * msh_psi_bnd(m);
        leakage *= msh_A(msh_I);

        double balance_error = (source_rate + incident) - (abs_rate + leakage);
        double rel_err       = balance_error / std::max(source_rate + incident, 1e-30);

        std::ofstream f(dir + "/" + msh_inp.name + "balance_table.txt");
        f << std::scientific << std::setprecision(6);
        f << "source_rate    = " << source_rate  << "\n"
          << "abs_rate       = " << abs_rate      << "\n"
          << "incident       = " << incident      << "\n"
          << "leakage        = " << leakage       << "\n"
          << "balance_error  = " << balance_error << "\n"
          << "relative_error = " << rel_err       << "\n";
        std::cout << "[balance] relative_error = " << rel_err << "\n";
    }

    //  angular_flux_boundary.csv 
    // Columns: direction, mu_m, w_m, psi_boundary
    void write_angular_flux_boundary(const std::string& dir)
    {
        std::ofstream f(dir + "/" + msh_inp.name + "angular_flux_boundary.csv");
        f << std::scientific << std::setprecision(8);
        f << "direction,mu_m,w_m,psi_boundary\n";
        for (int m = 0; m < msh_N_angle; ++m)
            f << m << "," << msh_mu(m) << "," << msh_w_angle(m)
              << "," << msh_psi_bnd(m) << "\n";
    }

    //  starting_direction.csv 
    // Row 0: r=0, psi_start_origin;  Rows 1..I: r_center[i], psi_start_cell[i]
    void write_starting_direction(const std::string& dir)
    {
        std::ofstream f(dir + "/" + msh_inp.name + "starting_direction.csv");
        f << std::scientific << std::setprecision(8);
        f << "cell,r_center_cm,psi_start\n";
        f << 0 << "," << 0.0 << "," << msh_psi_start_origin << "\n";
        for (int i = 0; i < msh_I; ++i)
            f << (i + 1) << "," << msh_r_ctr(i) << "," << msh_psi_start_cell(i) << "\n";
    }
};
