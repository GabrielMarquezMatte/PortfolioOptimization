#pragma once

#include <cstddef>
#include <span>
#include <vector>

namespace po::linalg
{

// Householder QR of a rows x cols matrix a (row-major, rows >= cols): a = Q R, with Q orthogonal, kept as its
// cols reflectors, and R upper triangular. Q is never formed.
struct Qr
{
    std::size_t rows = 0;
    std::size_t cols = 0;
    // Reflector k is the unit vector v of H_k = I - 2 v v', zero above entry k (cols x rows, row-major); a
    // zero vector is the identity.
    std::vector<double> reflectors;
    // R, cols x cols, row-major.
    std::vector<double> upper;
    double largest_column_norm = 0.0;
};

// The rank decision: a has full column rank when every |R_jj| exceeds
// kRankFactor * max(rows, cols) * DBL_EPSILON * the largest column norm of a. The rounding of Householder QR
// perturbs R by a small multiple of max(rows, cols) eps |a|, so a diagonal entry below the threshold is that
// rounding and the column depends on the earlier ones to working precision.
inline constexpr double kRankFactor = 16.0;

Qr householder_qr(std::span<const double> a, std::size_t rows, std::size_t cols);

bool has_full_column_rank(const Qr& qr);

// b (rows entries) becomes Q'b.
void apply_transposed(const Qr& qr, std::span<double> b);

// x (rows entries) becomes Q x: the reflectors applied last to first.
void apply(const Qr& qr, std::span<double> x);

// b (cols entries) becomes the x of R x = b; requires full column rank.
void solve_upper(const Qr& qr, std::span<double> b);

// b (cols entries) becomes the x of R'x = b; requires full column rank.
void solve_upper_transposed(const Qr& qr, std::span<double> b);

// The x (cols entries) that minimizes |a x - b|; requires full column rank.
std::vector<double> least_squares(const Qr& qr, std::vector<double> b);

// R' with each column signed so that the diagonal is nonnegative: the lower factor L of a'a = L L', computed
// without forming a'a.
std::vector<double> lower_factor(const Qr& qr);

} // namespace po::linalg
