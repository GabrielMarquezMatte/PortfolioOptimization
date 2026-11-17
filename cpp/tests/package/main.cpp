// A downstream project: the maximum Sharpe ratio portfolio of three assets and cash, through the C++ wrapper.
#include <cmath>
#include <cstdio>
#include <portopt/portopt.hpp>
#include <vector>

int main()
{
    const std::vector<double> mu{0.08, 0.12, 0.10};
    const std::vector<double> cov{0.040, 0.006, 0.010, 0.006, 0.090, 0.012, 0.010, 0.012, 0.060};
    po::Expected<po::Data> data = po::Data::create(3, mu, cov);
    if (!data)
    {
        std::fprintf(stderr, "FAIL: %s\n", data.error().message.c_str());
        return 1;
    }
    if (auto cash = data->set_risk_free(0.02); !cash)
    {
        std::fprintf(stderr, "FAIL: %s\n", cash.error().message.c_str());
        return 1;
    }
    const po::Expected<po::Result> result =
        po::solve(po::Spec{.objective = po::Objective::max_sharpe, .bounds = po::Bounds{0.0, 1.0}}, *data);
    if (!result)
    {
        std::fprintf(stderr, "FAIL: %s\n", result.error().message.c_str());
        return 1;
    }
    double sum = 0.0;
    for (const double weight : result->weights)
    {
        sum += weight;
    }
    if (result->status != po::SolveStatus::optimal || std::abs(sum - 1.0) > 1e-12)
    {
        std::fprintf(stderr, "FAIL: status %d, weights sum to %.17g\n", static_cast<int>(result->status), sum);
        return 1;
    }
    std::printf("PortOpt %.*s, maximum Sharpe weights:", static_cast<int>(po::version().size()), po::version().data());
    for (const double weight : result->weights)
    {
        std::printf(" %.6f", weight);
    }
    std::printf("\n");
    return 0;
}
