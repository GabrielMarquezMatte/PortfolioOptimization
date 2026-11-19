#pragma once

namespace po::portable
{

// The same bits on every platform: these use only +, -, *, / and sqrt, which IEEE 754 requires correctly
// rounded, so the result does not depend on the platform's libm. Outside the domain each returns the special
// value of C99 Annex F, as libm does.

// exp(x) with an error below 1 ulp for finite x (fdlibm 5.3 e_exp.c). Overflows to +inf above
// 7.09782712893383973096e+02 and underflows to +0 below -7.45133219101941108420e+02; exp(+inf) = +inf,
// exp(-inf) = +0 and exp(NaN) is NaN.
double exp(double x);

// log(x) with an error below 1 ulp for finite x > 0, subnormals included (fdlibm 5.3 e_log.c);
// log(+-0) = -inf, log(+inf) = +inf, and NaN for x < 0 and for NaN.
double log(double x);

// sqrt(a^2 + b^2) for finite a and b, with no intermediate overflow or underflow: both are scaled by the
// power of two that brings the larger to [1, 2), and the result is scaled back. +inf when either is infinite,
// even with the other NaN; NaN when either is NaN and neither is infinite.
double hypot(double a, double b);

} // namespace po::portable
