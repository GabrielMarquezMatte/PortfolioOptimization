#include <algorithm>
#include <cstddef>
#include <limits>
#include <span>
#include <string>
#include <vector>

#include "capi/error.hpp"
#include "capi/handles.hpp"
#include "estimation/bayes_stein.hpp"
#include "estimation/black_litterman.hpp"
#include "estimation/moments.hpp"
#include "estimation/posterior.hpp"

using po::capi::from_failure;
using po::capi::guard;
using po::capi::null_argument;
using po::capi::null_handle;
using po::capi::set_error;

namespace
{

po_status check_d(const po_data* d)
{
    if (d == nullptr)
    {
        return null_handle("d");
    }
    return PO_OK;
}

po_status check_d_and_v(const po_data* d, const po_views* v)
{
    po_status status = check_d(d);
    if (status != PO_OK)
    {
        return status;
    }
    if (v == nullptr)
    {
        return null_handle("v");
    }
    return PO_OK;
}

// The span of count values of a caller array; the pointer was checked, and count is a size the caller owes.
std::span<const double> array_of(const double* values, std::size_t count)
{
    return std::span<const double>(values, count);
}

// An output vector of the core, copied to the caller only after the core succeeded.
void copy_out(const std::vector<double>& values, double* out)
{
    std::ranges::copy(values, out);
}

} // namespace

extern "C" po_status po_sample_moments(size_t t, size_t n, const double* returns, double* out_mu, double* out_cov)
{
    return guard(
        [&]
        {
            if (returns == nullptr)
            {
                return null_argument("returns");
            }
            if (out_mu == nullptr)
            {
                return null_argument("out_mu");
            }
            if (out_cov == nullptr)
            {
                return null_argument("out_cov");
            }
            if (n != 0 && t > std::numeric_limits<std::size_t>::max() / n)
            {
                return set_error(PO_INVALID_ARGUMENT, "returns size is too large");
            }
            auto moments = po::estimation::sample_moments(array_of(returns, t * n), t, n);
            if (!moments.has_value())
            {
                return from_failure(moments.error());
            }
            copy_out(moments->mu, out_mu);
            copy_out(moments->covariance, out_cov);
            return PO_OK;
        });
}

extern "C" po_status po_bayes_stein(const po_data* d, size_t t, double* out_mu, double* out_grand_mean,
                                    double* out_weight)
{
    return guard(
        [&]
        {
            po_status status = check_d(d);
            if (status != PO_OK)
            {
                return status;
            }
            if (out_mu == nullptr)
            {
                return null_argument("out_mu");
            }
            if (out_grand_mean == nullptr)
            {
                return null_argument("out_grand_mean");
            }
            if (out_weight == nullptr)
            {
                return null_argument("out_weight");
            }
            auto shrunk = po::estimation::bayes_stein(d->data, t);
            if (!shrunk.has_value())
            {
                return from_failure(shrunk.error());
            }
            copy_out(shrunk->mu, out_mu);
            *out_grand_mean = shrunk->grand_mean;
            *out_weight = shrunk->weight;
            return PO_OK;
        });
}

extern "C" po_status po_bl_equilibrium(const po_data* d, const double* w_mkt, size_t n, double delta, double* out_pi)
{
    return guard(
        [&]
        {
            po_status status = check_d(d);
            if (status != PO_OK)
            {
                return status;
            }
            if (w_mkt == nullptr)
            {
                return null_argument("w_mkt");
            }
            if (out_pi == nullptr)
            {
                return null_argument("out_pi");
            }
            if (n != d->data.asset_count())
            {
                return set_error(PO_INVALID_ARGUMENT, "n must be " + std::to_string(d->data.asset_count()) +
                                                          ", got " + std::to_string(n));
            }
            auto pi = po::estimation::equilibrium(d->data, array_of(w_mkt, n), delta);
            if (!pi.has_value())
            {
                return from_failure(pi.error());
            }
            copy_out(*pi, out_pi);
            return PO_OK;
        });
}

extern "C" po_status po_bl_omega_he_litterman(const po_data* d, const po_views* v, double tau, double* out_omega)
{
    return guard(
        [&]
        {
            po_status status = check_d_and_v(d, v);
            if (status != PO_OK)
            {
                return status;
            }
            if (out_omega == nullptr)
            {
                return null_argument("out_omega");
            }
            auto omega = po::estimation::omega_he_litterman(d->data, v->views, tau);
            if (!omega.has_value())
            {
                return from_failure(omega.error());
            }
            copy_out(*omega, out_omega);
            return PO_OK;
        });
}

extern "C" po_status po_bl_omega_idzorek(const po_data* d, const po_views* v, double tau, const double* confidence,
                                         double* out_omega)
{
    return guard(
        [&]
        {
            po_status status = check_d_and_v(d, v);
            if (status != PO_OK)
            {
                return status;
            }
            if (confidence == nullptr)
            {
                return null_argument("confidence");
            }
            if (out_omega == nullptr)
            {
                return null_argument("out_omega");
            }
            auto omega = po::estimation::omega_idzorek(d->data, v->views, tau,
                                                       array_of(confidence, v->views.view_count()));
            if (!omega.has_value())
            {
                return from_failure(omega.error());
            }
            copy_out(*omega, out_omega);
            return PO_OK;
        });
}

extern "C" po_status po_bl_posterior(const po_data* d, const po_views* v, double tau, double* out_mu, double* out_cov)
{
    return guard(
        [&]
        {
            po_status status = check_d_and_v(d, v);
            if (status != PO_OK)
            {
                return status;
            }
            if (out_mu == nullptr)
            {
                return null_argument("out_mu");
            }
            auto posterior = po::estimation::posterior(d->data, v->views, tau, out_cov != nullptr);
            if (!posterior.has_value())
            {
                return from_failure(posterior.error());
            }
            copy_out(posterior->mu, out_mu);
            if (out_cov != nullptr)
            {
                copy_out(posterior->covariance, out_cov);
            }
            return PO_OK;
        });
}
