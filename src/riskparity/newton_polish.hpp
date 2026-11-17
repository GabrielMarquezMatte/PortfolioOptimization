#pragma once

#include "riskparity/kkt.hpp"

namespace po::riskparity
{

// One Newton solve of the KKT system at the weights of a point, which stay where they are, for its
// multipliers: the active bounds are those its bound multipliers name, the equalities those whose gradients
// over the free weights are independent.
struct MultiplierCorrection
{
    // The point to report: the weights of the point and, when accepted, its multipliers plus the correction:
    // lambda + dlambda for the equalities solved for, and for a held weight
    // beta_j = (grad f - sum_k lambda_k grad h_k)_j + (H dw)_j at the new lambda, H the Hessian of the
    // Lagrangian and dw the step of the free weights. Otherwise the multipliers of the point.
    KktPoint reported_point;
    // |dlambda|, and |dlambda_2| of a second solve from the corrected multipliers; NaN where that solve was
    // not made.
    double change;
    double second_change;
    // Both systems have full rank to working precision, the second correction is at most half the first (the
    // solve is in its region of convergence), the corrected multipliers leave a stationarity within
    // kKktNoiseFactor of its noise (where the system is so ill-conditioned that the correction is rounding,
    // it is not), and every held beta keeps the sign of its bound.
    bool accepted;
};

MultiplierCorrection corrected_multipliers(const KktProblem& problem, const KktPoint& point);

struct PolishResult
{
    // The weights of the start when no step was taken, those it names at a bound held there exactly. Its
    // multipliers are always the settled ones: least squares over the free weights for the equalities, the
    // stationarity rows for the bounds. The KKT gate checks this point; its multipliers are not the ones to
    // report.
    KktPoint least_squares_point;
    // corrected_multipliers at least_squares_point: correction.reported_point is what to report.
    MultiplierCorrection correction;
    // largest(residual_of(least_squares_point)): the KKT residual in units of its noise.
    double residual;
    // Newton steps taken to reach least_squares_point.
    int steps;
    // The residual of least_squares_point is below that of the start as it was given.
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
// converges to any stationary point: the start decides which. At the point it ends at, the multipliers are
// corrected once more by corrected_multipliers, which neither moves the weights nor decides any step.
PolishResult newton_polish(const KktProblem& problem, const KktPoint& start);

} // namespace po::riskparity
