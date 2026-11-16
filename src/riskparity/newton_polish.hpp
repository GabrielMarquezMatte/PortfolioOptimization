#pragma once

#include "riskparity/kkt.hpp"

namespace po::riskparity
{

struct PolishResult
{
    // The weights of the start when no step was taken, those it names at a bound held there exactly. Its
    // multipliers are always the settled ones: least squares over the free weights for the equalities, the
    // stationarity rows for the bounds.
    KktPoint point;
    // largest(residual_of(point)): the KKT residual in units of its noise.
    double residual;
    // Newton steps taken to reach point.
    int steps;
    // The residual of point is below that of the start as it was given.
    bool improved;
};

// Newton on the KKT system of the equalities and of the active bounds, with the exact Hessian of the
// Lagrangian, from a point near a KKT point. The active bounds are those the bound multipliers of start name:
// their weights are held at the bound and their multipliers are read from the stationarity. At every point
// the multipliers of the equalities are the least-squares ones over the free weights. A weight that a step
// takes out of the box is held at the bound it crossed, and a bound whose multiplier comes out with the wrong
// sign is released. That settling is done at the start too, whether or not a step follows. A step is
// accepted only when it decreases the KKT residual, feasibility included, measured in units of its noise; the
// first step that does not ends the polish, and so does a residual within its noise. A local method that
// converges to any stationary point: the start decides which.
PolishResult newton_polish(const KktProblem& problem, const KktPoint& start);

} // namespace po::riskparity
