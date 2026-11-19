#pragma once

namespace po::portable
{

// The same bits on every platform: these use only +, -, *, / and sqrt, which IEEE 754 requires correctly
// rounded, so the result does not depend on the platform's libm.

// exp(x) for finite x, with an error below 1 ulp (fdlibm 5.3 e_exp.c). Overflows to +inf above
// 7.09782712893383973096e+02 and underflows to +0 below -7.45133219101941108420e+02.
double exp(double x);

// log(x) for finite x > 0, subnormals included, with an error below 1 ulp (fdlibm 5.3 e_log.c).
double log(double x);

// sqrt(a^2 + b^2) for finite a and b, with no intermediate overflow or underflow: both are scaled by the
// power of two that brings the larger to [1, 2), and the result is scaled back.
double hypot(double a, double b);

} // namespace po::portable
