#include <portopt/portopt.h>

int po_c_header_check(void)
{
    po_status (*set_target_volatility)(po_spec*, double) = po_spec_set_target_volatility;
    po_status (*set_risk_budget)(po_spec*, const double*, size_t) = po_spec_set_risk_budget;
    po_status (*risk_contributions)(const po_data*, const double*, size_t, double*) = po_risk_contributions;
    (void)set_target_volatility;
    (void)set_risk_budget;
    (void)risk_contributions;
    return po_abi_version();
}
