#!/usr/bin/env python3
"""Generate plots for Problems 17a-d."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# ── 17a: spatial convergence ──────────────────────────────────────────────────
base = "17a"
fig, ax = plt.subplots(figsize=(6, 4))
for N, ls in [(50, '--'), (100, '-'), (200, ':')]:
    d = np.loadtxt(f"{base}/{N}_cells_scalar_flux.csv", delimiter=',', skiprows=1)
    ax.plot(d[:, 1], d[:, 2], ls=ls, label=f"N={N}")
ax.set_xlabel("r [cm]")
ax.set_ylabel(r"$\phi_0$ [cm$^{-2}$s$^{-1}$]")
ax.set_title("Problem 17a — Spatial convergence\n(S4, pure absorber, incident flux)")
ax.legend()
ax.grid(True, linestyle='--', color='gray', alpha=0.5)
fig.tight_layout()
fig.savefig(f"../docs/figs/plot_17a.png", dpi=150)
plt.close(fig)
print("17a done")

# ── 17b: compare with analytic φ = q/σ_a = 1 ─────────────────────────────────
base = "17b"
d = np.loadtxt(f"{base}/17b_scalar_flux.csv", delimiter=',', skiprows=1)
r = d[:, 1]
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].plot(r, d[:, 2], '-', label=r"S$_8$")
axes[0].axhline(1.0, color='k', ls='--', label="Analytic φ=1")
axes[0].set_xlabel("r [cm]")
axes[0].set_ylabel(r"$\phi_0$")
axes[0].set_title("Problem 17b — Scalar flux")
axes[0].set_ylim(0.9, 1.1)
axes[0].legend()
axes[1].semilogy(r, np.abs(d[:, 2] - 1.0), '-')
axes[1].set_xlabel("r [cm]")
axes[1].set_ylabel(r"$|\phi - \phi_{\rm analytic}|$")
axes[1].set_title("Problem 17b — Error vs analytic")
axes[0].grid(True, linestyle='--', color='gray', alpha=0.5)
axes[1].grid(True, linestyle='--', color='gray', alpha=0.5)
fig.tight_layout()
fig.savefig(f"../docs/figs/plot_17b.png", dpi=150)
plt.close(fig)
print("17b done")

# ── 17c: compare with 17a (100 cells, S4, isotropic scatter) ─────────────────
base_a = "17a"
base_c = "17c"
d_a = np.loadtxt(f"{base_a}/100_cells_scalar_flux.csv", delimiter=',', skiprows=1)
d_c = np.loadtxt(f"{base_c}/17c_scalar_flux.csv", delimiter=',', skiprows=1)
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(d_a[:, 1], d_a[:, 2], '-',  label=r"17a — isotropic (P0)")
ax.plot(d_c[:, 1], d_c[:, 2], '--', label=r"17c — P3 anisotropic")
ax.set_xlabel("r [cm]")
ax.set_ylabel(r"$\phi_0$ [cm$^{-2}$s$^{-1}$]")
ax.set_title("Problem 17c — Anisotropic vs isotropic scattering (S4)")
ax.legend()
ax.grid(True, linestyle='--', color='gray', alpha=0.5)
fig.tight_layout()
fig.savefig(f"{base_c}/plot_17c.png", dpi=150)
plt.close(fig)
print("17c done")

# ── 17d: compare with analytic diffusion solution ─────────────────────────────
base = "17d"
d = np.loadtxt(f"{base}/17d_scalar_flux.csv", delimiter=',', skiprows=1)
r = d[:, 1]

# Analytic: D = 1/(3*30) = 1/90, phi = C - (q/(6D))*r^2 = C - 15*r^2
# BC at r=1: -D*dphi/dr = <mu>*phi(1) => 1/3 = <mu>*(C-15)
# S8 <mu> = sum_{mu_m>0} mu_m w_m (biradian weights sum to 2)
mu_pos = np.array([0.183434642496, 0.525532409916, 0.796666477414, 0.960289856498])
w_pos  = np.array([0.362683783378, 0.313706645878, 0.222381034453, 0.101228536290])
mu_avg = np.dot(mu_pos, w_pos)
D = 1.0 / 90.0
# -D * d/dr(C - 15r^2)|_{r=1} = <mu> * phi(1)
# -D * (-30) = <mu> * (C - 15)
# 30D = <mu> * (C - 15)
# C - 15 = 30D / <mu>
C = 15.0 + 30.0 * D / mu_avg
r_fine = np.linspace(0, 1, 500)
phi_analytic = C - 15.0 * r_fine**2

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].plot(r, d[:, 2], '-',  label=r"S$_8$ + DSA")
axes[0].plot(r_fine, phi_analytic, 'k--', label=f"Diffusion (C={C:.4f})")
axes[0].set_xlabel("r [cm]")
axes[0].set_ylabel(r"$\phi_0$")
axes[0].set_title("Problem 17d — Scalar flux vs diffusion")
axes[0].legend()
axes[0].grid(True, linestyle='--', color='gray', alpha=0.5)

phi_at_r = C - 15.0 * r**2
axes[1].semilogy(r, np.abs(d[:, 2] - phi_at_r), '-')
axes[1].set_xlabel("r [cm]")
axes[1].set_ylabel(r"$|\phi - \phi_{\rm diffusion}|$")
axes[1].set_title("Problem 17d — Error vs analytic diffusion")
axes[1].grid(True, linestyle='--', color='gray', alpha=0.5)
fig.tight_layout()
fig.savefig(f"../docs/figs/plot_17d.png", dpi=150)
plt.close(fig)
print(f"17d done  (C={C:.6f}, <mu>={mu_avg:.6f})")
