#include "linalg/qr.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace po::linalg
{

namespace
{

double column_norm(std::span<const double> a, std::size_t rows, std::size_t cols, std::size_t j)
{
    double sum = 0.0;
    for (std::size_t i = 0; i < rows; ++i)
    {
        sum += a[i * cols + j] * a[i * cols + j];
    }
    return std::sqrt(sum);
}

// The reflector that maps column k of work, from row k down, onto -sign(x_k) |x| e_k; R_kk is that value.
double make_reflector(const std::vector<double>& work, Qr& qr, std::size_t k)
{
    std::size_t cols = qr.cols;
    std::span<double> v = std::span<double>(qr.reflectors).subspan(k * qr.rows, qr.rows);
    double norm_squared = 0.0;
    for (std::size_t i = k; i < qr.rows; ++i)
    {
        v[i] = work[i * cols + k];
        norm_squared += v[i] * v[i];
    }
    double norm = std::sqrt(norm_squared);
    if (norm == 0.0)
    {
        return 0.0;
    }
    double diagonal = -norm;
    if (v[k] < 0.0)
    {
        diagonal = norm;
    }
    v[k] -= diagonal;
    double length = std::sqrt(2.0 * norm * (norm + std::fabs(work[k * cols + k])));
    for (std::size_t i = k; i < qr.rows; ++i)
    {
        v[i] /= length;
    }
    return diagonal;
}

void reflect(std::span<const double> v, std::size_t k, std::span<double> x, std::size_t stride)
{
    double projection = 0.0;
    for (std::size_t i = k; i < v.size(); ++i)
    {
        projection += v[i] * x[i * stride];
    }
    for (std::size_t i = k; i < v.size(); ++i)
    {
        x[i * stride] -= 2.0 * projection * v[i];
    }
}

void solve_upper(const Qr& qr, std::span<double> b)
{
    std::size_t n = qr.cols;
    for (std::size_t i = n; i-- > 0;)
    {
        double sum = b[i];
        for (std::size_t k = i + 1; k < n; ++k)
        {
            sum -= qr.upper[i * n + k] * b[k];
        }
        b[i] = sum / qr.upper[i * n + i];
    }
}

} // namespace

Qr householder_qr(std::span<const double> a, std::size_t rows, std::size_t cols)
{
    Qr qr{rows, cols, std::vector<double>(cols * rows, 0.0), std::vector<double>(cols * cols, 0.0), 0.0};
    std::vector<double> work(a.begin(), a.end());
    for (std::size_t j = 0; j < cols; ++j)
    {
        qr.largest_column_norm = std::max(qr.largest_column_norm, column_norm(a, rows, cols, j));
    }
    for (std::size_t k = 0; k < cols; ++k)
    {
        double diagonal = make_reflector(work, qr, k);
        std::span<const double> v = std::span<const double>(qr.reflectors).subspan(k * rows, rows);
        for (std::size_t j = k + 1; j < cols; ++j)
        {
            reflect(v, k, std::span<double>(work).subspan(j), cols);
        }
        qr.upper[k * cols + k] = diagonal;
        for (std::size_t j = k + 1; j < cols; ++j)
        {
            qr.upper[k * cols + j] = work[k * cols + j];
        }
    }
    return qr;
}

bool has_full_column_rank(const Qr& qr)
{
    double size = static_cast<double>(std::max(qr.rows, qr.cols));
    double threshold = kRankFactor * size * std::numeric_limits<double>::epsilon() * qr.largest_column_norm;
    for (std::size_t j = 0; j < qr.cols; ++j)
    {
        if (!(std::fabs(qr.upper[j * qr.cols + j]) > threshold))
        {
            return false;
        }
    }
    return true;
}

void apply_transposed(const Qr& qr, std::span<double> b)
{
    for (std::size_t k = 0; k < qr.cols; ++k)
    {
        reflect(std::span<const double>(qr.reflectors).subspan(k * qr.rows, qr.rows), k, b, 1);
    }
}

std::vector<double> least_squares(const Qr& qr, std::vector<double> b)
{
    apply_transposed(qr, b);
    b.resize(qr.cols);
    solve_upper(qr, b);
    return b;
}

std::vector<double> lower_factor(const Qr& qr)
{
    std::size_t n = qr.cols;
    std::vector<double> lower(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        double sign = 1.0;
        if (qr.upper[i * n + i] < 0.0)
        {
            sign = -1.0;
        }
        for (std::size_t j = i; j < n; ++j)
        {
            lower[j * n + i] = sign * qr.upper[i * n + j];
        }
    }
    return lower;
}

} // namespace po::linalg
