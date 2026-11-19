#include <gtest/gtest.h>

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "random_covariance.hpp"

namespace
{

std::uint64_t bits_of(double x)
{
    return std::bit_cast<std::uint64_t>(x);
}

// The committed values are the output of the MinGW build; every platform must reproduce them, since the test
// instances are drawn from these generators.
constexpr std::array<std::uint64_t, 8> kNormals = {{
    0x3FF4B37D'0B4DAA8F,
    0x3FD97876'2F4A05DF,
    0x3FF1E599'FB914976,
    0xBFF7E03E'64E8AB39,
    0xBFD3C385'727C930F,
    0x3FF69D26'7030B33E,
    0x3FF97DCE'06BD632B,
    0x3FDCD46B'9F4117B9,
}};

constexpr std::uint64_t kExp10 = 0x3E60FA33'89D6EB32;

constexpr std::array<std::uint64_t, 5> kIndices = {{7, 14, 0, 3, 3}};

constexpr std::array<std::int64_t, 5> kIntegers = {{7, 0, 0, 1, 0}};

constexpr std::array<std::uint64_t, 16> kCovariance = {{
    0x40C0AA4E'6A8E260D, 0x40FD1849'2E26358D, 0xC0D8D191'CC4BB00C, 0x40FECAF2'6649C1A8,
    0x40FD1849'2E26358D, 0x413BB485'C06A5B68, 0xC119D124'68C2C907, 0x413D9F57'FEBF46F8,
    0xC0D8D191'CC4BB00C, 0xC119D124'68C2C907, 0x40F9F91B'D0740D83, 0xC11BDD3E'A3C37A78,
    0x40FECAF2'6649C1A8, 0x413D9F57'FEBF46F8, 0xC11BDD3E'A3C37A78, 0x413FB550'3E12C436,
}};

// The first 8 normals of seed 42, exp10(-7.5), then 5 indices below 17 and 5 integers in [-8, 8] from one engine of
// seed 42, and random_covariance(4, 7, 1e8), all by bits.
TEST(Generators, HaveTheCommittedBits)
{
    std::mt19937_64 normals(42);
    for (std::size_t k = 0; k < kNormals.size(); ++k)
    {
        EXPECT_EQ(bits_of(po::testing::standard_normal(normals)), kNormals[k]) << "normal " << k;
    }
    EXPECT_EQ(bits_of(po::testing::exp10(-7.5)), kExp10);
    std::mt19937_64 engine(42);
    for (std::size_t k = 0; k < kIndices.size(); ++k)
    {
        EXPECT_EQ(po::testing::uniform_index(engine, 17), kIndices[k]) << "index " << k;
    }
    for (std::size_t k = 0; k < kIntegers.size(); ++k)
    {
        EXPECT_EQ(po::testing::uniform_integer(engine, -8, 8), kIntegers[k]) << "integer " << k;
    }
    std::vector<double> covariance = po::testing::random_covariance(4, 7, 1e8);
    ASSERT_EQ(covariance.size(), kCovariance.size());
    for (std::size_t k = 0; k < kCovariance.size(); ++k)
    {
        EXPECT_EQ(bits_of(covariance[k]), kCovariance[k]) << "covariance " << k;
    }
}

} // namespace
