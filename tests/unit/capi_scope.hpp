#pragma once

#include <cstddef>
#include <cstdint>
#include <gtest/gtest.h>
#include <memory>
#include <portopt/portopt.h>
#include <span>
#include <vector>

namespace po::testing
{

struct SpecFree
{
    void operator()(po_spec* handle) const
    {
        po_spec_free(handle);
    }
};

struct DataFree
{
    void operator()(po_data* handle) const
    {
        po_data_free(handle);
    }
};

struct ResultFree
{
    void operator()(po_result* handle) const
    {
        po_result_free(handle);
    }
};

using SpecPtr = std::unique_ptr<po_spec, SpecFree>;
using DataPtr = std::unique_ptr<po_data, DataFree>;
using ResultPtr = std::unique_ptr<po_result, ResultFree>;

inline SpecPtr make_spec(int32_t objective = PO_OBJECTIVE_MIN_VARIANCE)
{
    po_spec* raw = nullptr;
    EXPECT_EQ(po_spec_create(objective, &raw), PO_OK);
    return SpecPtr(raw);
}

inline DataPtr make_data(std::size_t n, std::span<const double> mu, std::span<const double> covariance)
{
    po_data* raw = nullptr;
    const double* mu_pointer = mu.empty() ? nullptr : mu.data();
    const double* covariance_pointer = covariance.empty() ? nullptr : covariance.data();
    EXPECT_EQ(po_data_create(n, mu_pointer, covariance_pointer, &raw), PO_OK);
    return DataPtr(raw);
}

inline ResultPtr solve_ok(const po_spec* spec, const po_data* data)
{
    po_result* raw = nullptr;
    EXPECT_EQ(po_solve(spec, data, &raw), PO_OK);
    return ResultPtr(raw);
}

inline std::vector<double> weights_of(const po_result* result)
{
    std::size_t count = 0;
    EXPECT_EQ(po_result_weight_count(result, &count), PO_OK);
    std::vector<double> weights(count);
    EXPECT_EQ(po_result_weights(result, weights.data(), count), PO_OK);
    return weights;
}

inline std::int32_t status_of(const po_result* result)
{
    std::int32_t status = -1;
    EXPECT_EQ(po_result_status(result, &status), PO_OK);
    return status;
}

} // namespace po::testing
