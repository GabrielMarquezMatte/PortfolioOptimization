#include "riskparity/objective.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace po::riskparity
{
namespace
{

// p = Sigma w, v = w'p, c_i = w_i p_i / v and r = c - b at one point.
struct Moments
{
    std::vector<double> product;
    double variance;
    std::vector<double> contributions;
    std::vector<double> gap;
};

double dot(std::span<const double> a, std::span<const double> b)
{
    double total = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        total += a[i] * b[i];
    }
    return total;
}

} // namespace

// Dot2 of Ogita, Rump and Oishi (2005): each product and each partial sum is split exactly into its rounded
// value and its error (TwoProduct by fma, TwoSum), and the errors are added back once, so the result is as
// accurate as if summed in twice the working precision and then rounded. The splits are exact only when every
// operation is rounded as written: no contraction into fma and no reassociation.
double compensated_dot(std::span<const double> a, std::span<const double> b)
{
    double sum = 0.0;
    double correction = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        double product = a[i] * b[i];
        double product_error = std::fma(a[i], b[i], -product);
        double total = sum + product;
        double carried = total - sum;
        double sum_error = (sum - (total - carried)) + (product - carried);
        sum = total;
        correction += sum_error + product_error;
    }
    return sum + correction;
}

// The entries of Sigma x cancel by up to cond(Sigma), and an error of eps cond(Sigma) in p = Sigma w passes
// into grad f at the same relative size, above the KKT bands: each entry is a compensated dot product.
std::vector<double> covariance_times(const Data& data, std::span<const double> x)
{
    std::size_t n = data.asset_count();
    std::vector<double> product(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        product[i] = compensated_dot(data.covariance().subspan(i * n, n), x);
    }
    return product;
}

namespace
{

Moments moments_of(const Data& data, std::span<const double> budget, std::span<const double> w)
{
    std::size_t n = data.asset_count();
    Moments moments{covariance_times(data, w), 0.0, std::vector<double>(n), std::vector<double>(n)};
    moments.variance = dot(w, moments.product);
    for (std::size_t i = 0; i < n; ++i)
    {
        moments.contributions[i] = w[i] * moments.product[i] / moments.variance;
        moments.gap[i] = moments.contributions[i] - budget[i];
    }
    return moments;
}

// J_ij = (delta_ij p_i + w_i Sigma_ij - 2 c_i p_j) / v.
std::vector<double> jacobian_of(const Data& data, std::span<const double> w, const Moments& moments)
{
    std::size_t n = data.asset_count();
    std::span<const double> covariance = data.covariance();
    std::vector<double> jacobian(n * n);
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            double cross = w[i] * covariance[i * n + j] - 2.0 * moments.contributions[i] * moments.product[j];
            jacobian[i * n + j] = cross / moments.variance;
        }
        jacobian[i * n + i] += moments.product[i] / moments.variance;
    }
    return jacobian;
}

std::vector<double> gradient_of(std::span<const double> jacobian, std::span<const double> gap)
{
    std::size_t n = gap.size();
    std::vector<double> gradient(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            gradient[j] += 2.0 * jacobian[i * n + j] * gap[i];
        }
    }
    return gradient;
}

// a = sum_i r_i grad(w_i p_i) = r o p + Sigma (w o r).
std::vector<double> numerator_gradient(const Data& data, std::span<const double> w, const Moments& moments)
{
    std::size_t n = w.size();
    std::vector<double> weighted_gap(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        weighted_gap[i] = w[i] * moments.gap[i];
    }
    std::vector<double> numerator = covariance_times(data, weighted_gap);
    for (std::size_t i = 0; i < n; ++i)
    {
        numerator[i] += moments.gap[i] * moments.product[i];
    }
    return numerator;
}

// sum_i r_i hess(c_i) of c_i = N_i / v, N_i = w_i p_i, with a and k = c'r:
//   (diag(r) Sigma + Sigma diag(r)) / v - 2 (a p' + p a') / v^2 - 2 k Sigma / v + 8 k p p' / v^2.
std::vector<double> second_order_term(const Data& data, std::span<const double> w, const Moments& moments)
{
    std::size_t n = w.size();
    std::span<const double> covariance = data.covariance();
    std::span<const double> p = moments.product;
    std::span<const double> r = moments.gap;
    std::vector<double> a = numerator_gradient(data, w, moments);
    double v = moments.variance;
    double k = dot(moments.contributions, r);
    std::vector<double> term(n * n);
    for (std::size_t j = 0; j < n; ++j)
    {
        for (std::size_t l = 0; l < n; ++l)
        {
            double sigma = covariance[j * n + l];
            double gap_curvature = (r[j] + r[l] - 2.0 * k) * sigma / v;
            double cross = -2.0 * (a[j] * p[l] + p[j] * a[l]) / (v * v);
            double outer = 8.0 * k * p[j] * p[l] / (v * v);
            term[j * n + l] = gap_curvature + cross + outer;
        }
    }
    return term;
}

std::vector<double> hessian_of(const Data& data, std::span<const double> w, const Moments& moments,
                               std::span<const double> jacobian)
{
    std::vector<double> hessian = gauss_newton(jacobian, w.size());
    std::vector<double> term = second_order_term(data, w, moments);
    for (std::size_t k = 0; k < hessian.size(); ++k)
    {
        hessian[k] += 2.0 * term[k];
    }
    return hessian;
}

} // namespace

Linearization linearize(const Data& data, std::span<const double> budget, std::span<const double> w,
                        bool with_hessian)
{
    Moments moments = moments_of(data, budget, w);
    Linearization linearization{0.0, {}, jacobian_of(data, w, moments), {}};
    for (double gap : moments.gap)
    {
        linearization.value += gap * gap;
    }
    linearization.gradient = gradient_of(linearization.jacobian, moments.gap);
    if (with_hessian)
    {
        linearization.hessian = hessian_of(data, w, moments, linearization.jacobian);
    }
    return linearization;
}

namespace
{

constexpr std::size_t kPanel = 4;
// Below this size, packing the panels costs more than the products that the symmetry saves
// (docs/performance/benchmarks.md).
constexpr std::size_t kSymmetricFrom = 16;

// Entry (j, k) of 2 J'J is the sum over i of (2 J_ij) J_ik, from +0 in increasing i.
std::vector<double> full_gauss_newton(std::span<const double> jacobian, std::size_t n)
{
    std::vector<double> curvature(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            double scaled = 2.0 * jacobian[i * n + j];
            for (std::size_t k = 0; k < n; ++k)
            {
                curvature[j * n + k] += scaled * jacobian[i * n + k];
            }
        }
    }
    return curvature;
}

// 2 J_ij is exact for every entry, so (2 J_ij) J_ik and (2 J_ik) J_ij are the same rounding of 2 J_ij J_ik
// and entry (k, j) of 2 J'J has the bits of entry (j, k). An entry above max / 2 (or one that is not finite)
// breaks it.
bool doubles_exactly(std::span<const double> jacobian)
{
    constexpr double kHalfLargest = std::numeric_limits<double>::max() / 2.0;
    return std::ranges::all_of(jacobian, [](double entry) { return std::fabs(entry) <= kHalfLargest; });
}

// J by panels of kPanel columns: panel p holds J_i,(p kPanel + b) at i kPanel + b, with zeros past column n.
std::vector<double> column_panels(std::span<const double> jacobian, std::size_t n)
{
    std::size_t count = (n + kPanel - 1) / kPanel;
    std::vector<double> panels(count * n * kPanel, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            panels[((j / kPanel) * n + i) * kPanel + j % kPanel] = jacobian[i * n + j];
        }
    }
    return panels;
}

// 2 J_ij of a panel.
void double_panel(std::span<const double> panel, std::span<double> doubled)
{
    for (std::size_t e = 0; e < panel.size(); ++e)
    {
        doubled[e] = 2.0 * panel[e];
    }
}

// The entries (j, k) of the rows of the doubled panel and the columns of panel, each summed from +0 in
// increasing i as full_gauss_newton sums it, stored at (j, k) and at (k, j).
void add_block(std::span<const double> doubled, std::span<const double> panel, std::size_t first_row,
               std::size_t first_column, std::span<double> curvature)
{
    std::size_t n = panel.size() / kPanel;
    std::array<double, kPanel * kPanel> sums{};
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t a = 0; a < kPanel; ++a)
        {
            double scale = doubled[i * kPanel + a];
            for (std::size_t b = 0; b < kPanel; ++b)
            {
                sums[a * kPanel + b] += scale * panel[i * kPanel + b];
            }
        }
    }
    std::size_t rows = std::min(kPanel, n - first_row);
    std::size_t columns = std::min(kPanel, n - first_column);
    for (std::size_t a = 0; a < rows; ++a)
    {
        for (std::size_t b = 0; b < columns; ++b)
        {
            curvature[(first_row + a) * n + first_column + b] = sums[a * kPanel + b];
            curvature[(first_column + b) * n + first_row + a] = sums[a * kPanel + b];
        }
    }
}

// The blocks of kPanel x kPanel entries on and above the diagonal, mirrored below it.
std::vector<double> symmetric_gauss_newton(std::span<const double> jacobian, std::size_t n)
{
    std::vector<double> curvature(n * n, 0.0);
    std::vector<double> panels = column_panels(jacobian, n);
    std::vector<double> doubled(kPanel * n);
    std::size_t panel_size = kPanel * n;
    std::size_t count = (n + kPanel - 1) / kPanel;
    std::span<const double> all(panels);
    for (std::size_t p = 0; p < count; ++p)
    {
        double_panel(all.subspan(p * panel_size, panel_size), doubled);
        for (std::size_t q = p; q < count; ++q)
        {
            add_block(doubled, all.subspan(q * panel_size, panel_size), p * kPanel, q * kPanel, curvature);
        }
    }
    return curvature;
}

} // namespace

std::vector<double> gauss_newton(std::span<const double> jacobian, std::size_t n)
{
    if (n < kSymmetricFrom || !doubles_exactly(jacobian))
    {
        return full_gauss_newton(jacobian, n);
    }
    return symmetric_gauss_newton(jacobian, n);
}

// With s = step, q = Sigma s and N_i = w_i p_i: dN_i = s_i p_i + w_i q_i + s_i q_i, dv = 2 p's + s'q,
// dc_i = (dN_i v - N_i dv) / (v (v + dv)) and df = sum_i dc_i (2 r_i + dc_i).
Change value_change(const Data& data, std::span<const double> budget, std::span<const double> w,
                    std::span<const double> step)
{
    Moments moments = moments_of(data, budget, w);
    std::vector<double> q = covariance_times(data, step);
    Change change{0.0, 2.0 * dot(moments.product, step) + dot(step, q)};
    double v = moments.variance;
    double updated = v + change.variance;
    for (std::size_t i = 0; i < w.size(); ++i)
    {
        double numerator = w[i] * moments.product[i];
        double numerator_change = step[i] * moments.product[i] + w[i] * q[i] + step[i] * q[i];
        double contribution_change = (numerator_change * v - numerator * change.variance) / (v * updated);
        change.value += contribution_change * (2.0 * moments.gap[i] + contribution_change);
    }
    return change;
}

} // namespace po::riskparity
