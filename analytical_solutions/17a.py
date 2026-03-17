import numpy as np
from scipy import integrate

r_b     = 1.0   # cm
sigma_t = 1.0   # cm^-1
N_list  = [50, 100, 200]

def scalar_flux(r):
    if r == 0.0:
        return 4.0 * np.exp(-sigma_t * r_b)
    result, _ = integrate.quad(
        lambda mu: 2.0 * np.exp(-sigma_t * (r * mu + np.sqrt(r_b**2 - r**2 * (1 - mu**2)))),
        -1.0, 1.0, limit=200, epsabs=1e-12, epsrel=1e-12
    )
    return result

for N in N_list:
    dr      = r_b / N
    r_cells = np.array([(i + 0.5) * dr for i in range(N)])
    phi     = np.array([scalar_flux(r) for r in r_cells])

    with open(f"../examples/17a/phi_analytical_N{N}.ult", "w") as f:
        f.write("# phi-0 vs. cell center\n")
        for r, p in zip(r_cells, phi):
            f.write(f"{r:.8e} {p:.8e}\n")

    print(f"Written phi_analytical_N{N}.ult")