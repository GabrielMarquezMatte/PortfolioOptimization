#include "core/portable_math.hpp"

#include <gtest/gtest.h>

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>

namespace
{

struct Unary
{
    double x;
    std::uint64_t bits;
};

struct Binary
{
    double a;
    double b;
    std::uint64_t bits;
};

std::uint64_t bits_of(double x)
{
    return std::bit_cast<std::uint64_t>(x);
}

// Keeps the compiler from evaluating a call with constant arguments at compile time, where GCC rounds
// std::fma with MPFR instead of calling the platform's fma.
double opaque(double x)
{
    volatile double stored = x;
    return stored;
}

// The committed bits are the output of the MinGW build. Each is checked against a 60-digit reference,
//     python -c "from mpmath import mp, mpf, exp, log, sqrt; mp.dps = 60; print(exp(mpf(X)))"
// (log(mpf(X)), sqrt(mpf(A)**2 + mpf(B)**2) for the others), at a distance of |result - exact| in units of
// 2^(max(floor(log2 |exact|), -1022) - 52). The largest distance in these tables is 0.6745 ulp, at exp(1);
// 20000 random arguments of each kind reach 0.8052 ulp for exp, 0.8102 for log and 1.1062 for hypot.
constexpr std::array<Unary, 11> kExp = {{
    {0.3, 0x3FF59905'8C8C1A96},
    {-0.3, 0x3FE7B4C8'69C37C05},
    {0.5, 0x3FFA6129'8E1E069C},
    {1.0, 0x4005BF0A'8B14576A},
    {-0.9, 0x3FDA053C'C0086E1E},
    {2.5, 0x40285D6F'D931E0BB},
    {-20.5, 0x3E157A3A'FEED00AB},
    {100.0, 0x48F3494A'9B171BF5},
    {700.0, 0x7F0D945D'F4F8EC8E},
    {-740.0, 0x00000000'00000055},
    {1e-10, 0x3FF00000'0006DF38},
}};

constexpr std::array<Unary, 11> kLog = {{
    {0.1, 0xC0026BB1'BBB55515},
    {0.5, 0xBFE62E42'FEFA39EF},
    {0.75, 0xBFD26962'1134DB92},
    {1.4, 0x3FD588C2'D913348E},
    {2.0, 0x3FE62E42'FEFA39EF},
    {3.0, 0x3FF193EA'7AAD030A},
    {10.0, 0x40026BB1'BBB55516},
    {1 + 0x1p-30, 0x3E0FFFFF'FFC00000},
    {1e-300, 0xC0859634'47F87FB5},
    {1e300, 0x40859634'47F87FB5},
    {5e-320, 0xC086F9BE'0F7CEE66},
}};

constexpr std::array<Binary, 10> kHypot = {{
    {1.0, 1.0, 0x3FF6A09E'667F3BCD},
    {1.0, 2.0, 0x4001E377'9B97F4A8},
    {0.1, 0.2, 0x3FCC9F25'C5BFEDDA},
    {-2.5, 7.1, 0x401E1BF0'428EA595},
    {1e10, 1.0, 0x4202A05F'20000000},
    {1e-200, 3e-200, 0x16835D52'44B69495},
    {1e308, 1e308, 0x7FE92C80'954C51F5},
    {0x1p-1060, 0x3p-1060, 0x00000000'0000CA63},
    {0.3, 0.4, 0x3FE00000'00000000},
    {5.0, 12.0, 0x402A0000'00000000},
}};

TEST(PortableMath, ExpHasTheCommittedBits)
{
    for (const Unary& entry : kExp)
    {
        EXPECT_EQ(bits_of(po::portable::exp(entry.x)), entry.bits) << std::hexfloat << entry.x;
    }
}

TEST(PortableMath, LogHasTheCommittedBits)
{
    for (const Unary& entry : kLog)
    {
        EXPECT_EQ(bits_of(po::portable::log(entry.x)), entry.bits) << std::hexfloat << entry.x;
    }
}

TEST(PortableMath, HypotHasTheCommittedBits)
{
    for (const Binary& entry : kHypot)
    {
        EXPECT_EQ(bits_of(po::portable::hypot(entry.a, entry.b)), entry.bits)
            << std::hexfloat << entry.a << ", " << entry.b;
    }
}

TEST(PortableMath, HypotScalesExtremeArguments)
{
    EXPECT_EQ(bits_of(po::portable::hypot(0x3p-1074, 0x4p-1074)), bits_of(0x5p-1074));
    EXPECT_EQ(bits_of(po::portable::hypot(0x1p-1074, 0x1p-1074)), bits_of(0x1p-1074));
    EXPECT_EQ(bits_of(po::portable::hypot(1e-300, 1e300)), bits_of(1e300));
    EXPECT_EQ(bits_of(po::portable::hypot(0.0, 0.0)), bits_of(0.0));
    EXPECT_EQ(bits_of(po::portable::hypot(0.0, -3.0)), bits_of(3.0));
    EXPECT_EQ(bits_of(po::portable::hypot(3.0, 4.0)), bits_of(5.0));
    EXPECT_EQ(bits_of(po::portable::hypot(0x1p1000, 0x1p1000)), bits_of(0x1.6a09e667f3bcdp1000));
}

// The 60-digit references are
//     log(2^-1074) = -744.440071921381262314107298446081634113087144302914142926,
//     log(2^-1022) = -708.396418532264106224411228130256452573161137316180869711;
// the committed bits are the nearest double to each, 0.3890 and 0.2417 ulp away.
TEST(PortableMath, LogHandlesSubnormalsAndOne)
{
    EXPECT_EQ(bits_of(po::portable::log(1.0)), bits_of(0.0));
    EXPECT_EQ(bits_of(po::portable::log(std::numeric_limits<double>::denorm_min())),
              std::uint64_t{0xC0874385'446D71C3});
    EXPECT_EQ(bits_of(po::portable::log(0x1p-1022)), std::uint64_t{0xC086232B'DD7ABCD2});
}

// Each exact result needs the product unrounded: a * b + c rounded twice gives another value in all four.
TEST(PortableMath, FmaIsCorrectlyRounded)
{
    EXPECT_EQ(bits_of(std::fma(opaque(1 + 0x1p-52), opaque(1 - 0x1p-52), opaque(-1.0))), bits_of(-0x1p-104));
    EXPECT_EQ(bits_of(std::fma(opaque(0x1.0000001p0), opaque(0x1.0000001p0), opaque(-0x1.0000002p0))),
              bits_of(0x1p-56));
    // The product is a tie that the tiny addend breaks upward.
    EXPECT_EQ(bits_of(std::fma(opaque(1 + 0x1p-30), opaque(1 + 0x1p-23), opaque(0x1p-100))),
              bits_of(0x1.0000020400001p0));
    // The product is a tie at an odd last bit, so rounding it alone goes up; the tiny negative addend sends
    // it down.
    EXPECT_EQ(bits_of(std::fma(opaque(1 + 0x1p-30), opaque(1 + 0x1p-22 + 0x1p-23), opaque(-0x1p-100))),
              bits_of(0x1.0000060400001p0));
}

} // namespace
