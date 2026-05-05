"""
Compute dR/dsigma_a via first-order perturbation theory.

Formula (Lec 24):
    dR/dsigma_a = -<psi, psi†>
                = -sum_i  V_i  sum_m  w_m  psi_{i,m}  psi†_{i,m}

where
    psi_{i,m}   = cell-averaged forward angular flux   (from *angular_flux_cell.csv)
    psi†_{i,m}  = cell-averaged adjoint angular flux
                = hat_psi_{i, m_bar}  (adjoint output at reversed direction m_bar = N-1-m)
    V_i         = cell volume = (4pi/3)(r_{i+1/2}^3 - r_{i-1/2}^3)
    w_m         = Gauss-Legendre quadrature weight

Note: the inner product uses dmu (not dOmega = 2pi dmu) because the response
R = iint_{mu>0} psi mu dmu dA is defined with the cosine measure, matching
the solver's leakage convention (leakage = A sum_{m>0} mu_m w_m psi_bnd).

Usage:
    python compute_dR_dsigma_a.py \
        --fwd  2a/2a_angular_flux_cell.csv \
        --adj  2b/2b_angular_flux_cell.csv \
        --quad 2a/2a_angular_flux_boundary.csv \
        --ncells 50 --R 1.0
"""

import argparse
import csv
import math


def load_angular_flux_cell(path):
    """Return dict[(cell, direction)] = psi_cell."""
    data = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            key = (int(row["cell"]), int(row["direction"]))
            data[key] = float(row["psi_cell"])
    return data


def load_quadrature(path):
    """Return (mu, w) dicts keyed by direction index."""
    mu, w = {}, {}
    with open(path) as f:
        for row in csv.DictReader(f):
            m = int(row["direction"])
            mu[m] = float(row["mu_m"])
            w[m]  = float(row["w_m"])
    return mu, w


def cell_volumes(n_cells, R):
    """Spherical shell volumes V_i = (4pi/3)(r_{i+1/2}^3 - r_{i-1/2}^3)."""
    h = R / n_cells
    V = {}
    for i in range(n_cells):
        r_in  = i * h
        r_out = (i + 1) * h
        V[i] = (4 * math.pi / 3) * (r_out**3 - r_in**3)
    return V


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fwd",    required=True, help="Forward angular_flux_cell.csv")
    parser.add_argument("--adj",    required=True, help="Adjoint angular_flux_cell.csv")
    parser.add_argument("--quad",   required=True, help="angular_flux_boundary.csv (for w_m)")
    parser.add_argument("--ncells", type=int,   required=True, help="Number of spatial cells")
    parser.add_argument("--R",      type=float, required=True, help="Outer radius [cm]")
    args = parser.parse_args()

    psi_fwd = load_angular_flux_cell(args.fwd)
    psi_adj = load_angular_flux_cell(args.adj)
    mu, w   = load_quadrature(args.quad)
    V       = cell_volumes(args.ncells, args.R)

    N = len(mu)  # total number of directions (e.g. 16 for S16)

    # dR/dsigma_a = -sum_i V_i sum_m w_m psi_fwd[i,m] * psi_adj[i, N-1-m]
    #
    # psi†[i,m] = hat_psi[i, m_bar] where m_bar = N-1-m reverses the GL index,
    # mapping mu_m -> -mu_m (adjoint direction).
    total = 0.0
    for i in range(args.ncells):
        for m in range(N):
            m_bar  = N - 1 - m
            total += V[i] * w[m] * psi_fwd[(i, m)] * psi_adj[(i, m_bar)]

    dR_dsigma_a = -total

    print(f"N (SN order)     = {N}")
    print(f"Cells            = {args.ncells}")
    print(f"R                = {args.R} cm")
    print(f"dR/dsigma_a      = {dR_dsigma_a:.10f} cm")
    print()
    print("To predict delta_R for a perturbation delta_sigma_a:")
    print("    delta_R_pred = dR/dsigma_a * delta_sigma_a")


if __name__ == "__main__":
    main()
