#include "Functions.hpp"
// ---------------------------------------------------------------------------
//  Gauss-Legendre quadrature on [-1, 1] via MFEM
//  Returns N points and weights; weights sum to 2.
// ---------------------------------------------------------------------------
void GaussLegendreQuadrature(int N,
                              std::vector<double>& points,
                              std::vector<double>& weights)
{
    // MFEM's 1-D segment reference element spans [0,1].
    // An integration rule of order (2N-1) uses exactly N points and
    // integrates polynomials of degree up to 2N-1 exactly.
    const mfem::IntegrationRule& ir =
        mfem::IntRules.Get(mfem::Geometry::SEGMENT, 2 * N - 1);

    const int np = ir.GetNPoints();
    points.resize(np);
    weights.resize(np);

    for (int i = 0; i < np; ++i)
    {
        const mfem::IntegrationPoint& ip = ir.IntPoint(i);
        // Map from [0,1] -> [-1,1]:  x_new = 2*x - 1,  w_new = 2*w
        // MFEM weights on [0,1] sum to 1, so transformed weights sum to 2.
        points[i]  = 2.0 * ip.x - 1.0;
        weights[i] = 2.0 * ip.weight;
    }
}
double Pn(unsigned int n, double x)
{
    return std::legendre(n, x);
}