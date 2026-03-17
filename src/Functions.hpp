#pragma once

#include <cmath>
#include <mfem.hpp>
#include "axom/inlet.hpp"


// ---------------------------------------------------------------------------
//  DECLARATION: Gauss-Legendre quadrature on [-1, 1] via MFEM
//  Returns N points and weights; weights sum to 2.
// ---------------------------------------------------------------------------
void GaussLegendreQuadrature(int N,
                              std::vector<double>& abscissae,
                              std::vector<double>& weights);

double Pn(unsigned int n, double x);