// exp and log are ported from fdlibm 5.3 (https://www.netlib.org/fdlibm/), whose notices follow.
//
// e_exp.c:
// ====================================================
// Copyright (C) 2004 by Sun Microsystems, Inc. All rights reserved.
//
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
//
// e_log.c:
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunSoft, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================

#include "core/portable_math.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

// fdlibm's error bounds hold only with every operation rounded as written. GCC and Clang get
// -ffp-contract=off from the build.
#ifdef _MSC_VER
#pragma fp_contract(off)
#endif

namespace po::portable
{
namespace
{

constexpr double from_bits(std::uint64_t bits)
{
    return std::bit_cast<double>(bits);
}

constexpr double kInfinity = std::numeric_limits<double>::infinity();

// The constants are the hexadecimal words of fdlibm, which it states are the intended values.
constexpr double kLn2Hi = from_bits(0x3FE62E42'FEE00000);
constexpr double kLn2Lo = from_bits(0x3DEA39EF'35793C76);

constexpr double kExpOverflowThreshold = from_bits(0x40862E42'FEFA39EF);
constexpr double kExpUnderflowThreshold = from_bits(0xC0874910'D52D3051);
constexpr double kInvLn2 = from_bits(0x3FF71547'652B82FE);
constexpr double kTwoToMinus1000 = from_bits(0x01700000'00000000);
constexpr std::array<double, 2> kHalf = {0.5, -0.5};
constexpr std::array<double, 2> kSignedLn2Hi = {kLn2Hi, -kLn2Hi};
constexpr std::array<double, 2> kSignedLn2Lo = {kLn2Lo, -kLn2Lo};
constexpr double kP1 = from_bits(0x3FC55555'5555553E);
constexpr double kP2 = from_bits(0xBF66C16C'16BEBD93);
constexpr double kP3 = from_bits(0x3F11566A'AF25DE2C);
constexpr double kP4 = from_bits(0xBEBBBD41'C5D26BF1);
constexpr double kP5 = from_bits(0x3E663769'72BEA4D0);

constexpr double kTwo54 = from_bits(0x43500000'00000000);
constexpr double kOneThird = from_bits(0x3FD55555'55555555);
constexpr double kLg1 = from_bits(0x3FE55555'55555593);
constexpr double kLg2 = from_bits(0x3FD99999'9997FA04);
constexpr double kLg3 = from_bits(0x3FD24924'94229359);
constexpr double kLg4 = from_bits(0x3FCC71C5'1D8E78AF);
constexpr double kLg5 = from_bits(0x3FC74664'96CB03DE);
constexpr double kLg6 = from_bits(0x3FC39A09'D078C69F);
constexpr double kLg7 = from_bits(0x3FC2F112'DF3E5244);

std::uint32_t high_word(double x)
{
    return static_cast<std::uint32_t>(std::bit_cast<std::uint64_t>(x) >> 32);
}

double with_high_word(double x, std::uint32_t high)
{
    std::uint64_t low = std::bit_cast<std::uint64_t>(x) & 0xFFFFFFFF;
    return std::bit_cast<double>((static_cast<std::uint64_t>(high) << 32) | low);
}

// Adds k to the exponent of y; exact while the result stays normal.
double add_to_exponent(double y, int k)
{
    std::uint64_t shift = static_cast<std::uint64_t>(static_cast<std::int64_t>(k)) << 52;
    return std::bit_cast<double>(std::bit_cast<std::uint64_t>(y) + shift);
}

struct ExpReduction
{
    double hi;
    double lo;
    int k;
};

// x = k ln2 + (hi - lo) for 0.5 ln2 < |x| <= 709.78.
ExpReduction reduce_exp_argument(double x, std::uint32_t abs_high, std::size_t sign)
{
    if (abs_high < 0x3FF0A2B2)
    {
        return {x - kSignedLn2Hi[sign], kSignedLn2Lo[sign], 1 - static_cast<int>(sign + sign)};
    }
    int k = static_cast<int>(kInvLn2 * x + kHalf[sign]);
    double t = k;
    // t * kLn2Hi is exact: kLn2Hi has 32 significant bits and |k| < 2^11.
    return {x - t * kLn2Hi, t * kLn2Lo, k};
}

// r - r^2 (P1 + P2 r^2 + ... + P5 r^8), the R1(r) of fdlibm.
double exp_correction(double r)
{
    double t = r * r;
    return r - t * (kP1 + t * (kP2 + t * (kP3 + t * (kP4 + t * kP5))));
}

double scale_by_power_of_two(double y, int k)
{
    if (k >= -1021)
    {
        return add_to_exponent(y, k);
    }
    return add_to_exponent(y, k + 1000) * kTwoToMinus1000;
}

double log_near_one(double f, int k)
{
    double dk = k;
    if (f == 0.0 && k == 0)
    {
        return 0.0;
    }
    if (f == 0.0)
    {
        return dk * kLn2Hi + dk * kLn2Lo;
    }
    double r = f * f * (0.5 - kOneThird * f);
    if (k == 0)
    {
        return f - r;
    }
    return dk * kLn2Hi - ((r - dk * kLn2Lo) - f);
}

// log(1 + f) = f - (hfsq - s (hfsq + R)), the form fdlibm uses for its largest |f|.
double log_with_half_square(double f, int k, double s, double r)
{
    double dk = k;
    double hfsq = 0.5 * f * f;
    if (k == 0)
    {
        return f - (hfsq - s * (hfsq + r));
    }
    return dk * kLn2Hi - ((hfsq - (s * (hfsq + r) + dk * kLn2Lo)) - f);
}

// log(2^k (1 + f)) for sqrt(2)/2 <= 1 + f < sqrt(2) and |f| >= 2^-20; fraction_high is the high fraction
// word of the argument before normalization.
double log_of_reduced(double f, int k, std::int32_t fraction_high)
{
    double s = f / (2.0 + f);
    double dk = k;
    double z = s * s;
    double w = z * z;
    double t1 = w * (kLg2 + w * (kLg4 + w * kLg6));
    double t2 = z * (kLg1 + w * (kLg3 + w * (kLg5 + w * kLg7)));
    double r = t2 + t1;
    if (fraction_high >= 0x6147A && fraction_high <= 0x6B851)
    {
        return log_with_half_square(f, k, s, r);
    }
    if (k == 0)
    {
        return f - s * (f - r);
    }
    return dk * kLn2Hi - ((s * (f - r) - dk * kLn2Lo) - f);
}

// log(2^k x) for a normal x > 0.
double log_of_normal(double x, int k)
{
    auto high = static_cast<std::int32_t>(high_word(x));
    k += (high >> 20) - 1023;
    std::int32_t fraction_high = high & 0x000FFFFF;
    // Picks x or x/2, whichever puts the mantissa 1 + f in [sqrt(2)/2, sqrt(2)).
    std::int32_t halve = (fraction_high + 0x95F64) & 0x100000;
    double normalized = with_high_word(x, static_cast<std::uint32_t>(fraction_high | (halve ^ 0x3FF00000)));
    k += halve >> 20;
    double f = normalized - 1.0;
    if ((0x000FFFFF & (2 + fraction_high)) < 3)
    {
        return log_near_one(f, k);
    }
    return log_of_reduced(f, k, fraction_high);
}

// 2^exponent for -1074 <= exponent <= 1023.
double power_of_two(int exponent)
{
    if (exponent < -1022)
    {
        return from_bits(std::uint64_t{1} << (exponent + 1074));
    }
    return from_bits(static_cast<std::uint64_t>(exponent + 1023) << 52);
}

// 2^extra_exponent sqrt(larger^2 + smaller^2) for a normal larger >= smaller >= 0.
double hypot_of_normal(double larger, double smaller, int extra_exponent)
{
    int exponent = static_cast<int>(std::bit_cast<std::uint64_t>(larger) >> 52) - 1023;
    double scale = power_of_two(-exponent);
    double x = larger * scale;
    double y = smaller * scale;
    return std::sqrt(x * x + y * y) * power_of_two(exponent + extra_exponent);
}

// C99 Annex F outside (0, +inf): -inf at +-0, NaN below 0 and at NaN, +inf at +inf.
double log_outside_domain(double x)
{
    if (x == 0.0)
    {
        return -kInfinity;
    }
    if (x < 0.0)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return x + x;
}

// C99 Annex F for |a| or |b| infinite or NaN: +inf when either is infinite, even with the other NaN; NaN
// otherwise.
double hypot_outside_domain(double abs_a, double abs_b)
{
    if (abs_a == kInfinity || abs_b == kInfinity)
    {
        return kInfinity;
    }
    return abs_a + abs_b;
}

} // namespace

double exp(double x)
{
    if (std::isnan(x))
    {
        return x + x;
    }
    std::uint32_t high = high_word(x);
    std::size_t sign = high >> 31;
    std::uint32_t abs_high = high & 0x7FFFFFFF;
    if (x > kExpOverflowThreshold)
    {
        return std::numeric_limits<double>::infinity();
    }
    if (x < kExpUnderflowThreshold)
    {
        return 0.0;
    }
    if (abs_high < 0x3E300000)
    {
        return 1.0 + x;
    }
    if (abs_high <= 0x3FD62E42)
    {
        double c = exp_correction(x);
        return 1.0 - ((x * c) / (c - 2.0) - x);
    }
    ExpReduction reduced = reduce_exp_argument(x, abs_high, sign);
    double r = reduced.hi - reduced.lo;
    double c = exp_correction(r);
    double y = 1.0 - ((reduced.lo - (r * c) / (2.0 - c)) - reduced.hi);
    return scale_by_power_of_two(y, reduced.k);
}

double log(double x)
{
    if (!(x > 0.0 && x < kInfinity))
    {
        return log_outside_domain(x);
    }
    if (x < 0x1p-1022)
    {
        return log_of_normal(x * kTwo54, -54);
    }
    return log_of_normal(x, 0);
}

double hypot(double a, double b)
{
    double abs_a = std::abs(a);
    double abs_b = std::abs(b);
    if (!(abs_a < kInfinity && abs_b < kInfinity))
    {
        return hypot_outside_domain(abs_a, abs_b);
    }
    double larger = std::max(abs_a, abs_b);
    double smaller = std::min(abs_a, abs_b);
    if (larger == 0.0)
    {
        return 0.0;
    }
    if (larger < 0x1p-1022)
    {
        return hypot_of_normal(larger * kTwo54, smaller * kTwo54, -54);
    }
    return hypot_of_normal(larger, smaller, 0);
}

} // namespace po::portable
