# Viva notes

## The ray-distance idea

A telescope records a direction. If its measured unit direction is u, every point on that line has the form R + λu, where R is the heliocentric telescope position and λ is an unknown distance. For a predicted Eros position p, the closest point on the line is R + u·(p-R)u. The leftover sideways vector is (I-uuᵀ)(p-R). Its norm is the miss distance.

Tiny 2D example: let u=(1,0), R=(0,0), and p=(3,4). The predicted point is 3 units along the ray and 4 units sideways. The perpendicular miss is 4, not the 5-unit distance from the telescope.

## Why the problem is convex

For a fixed measured direction, the projector P=I-uuᵀ is known. The path p(τ)=r₀+v₀τ+½aτ²+… is linear in its unknown coefficients. Stacking the projected residuals gives Ac-b; the objective is ||Ac-b||², a convex quadratic with Hessian 2AᵀA. Least squares returns a global minimizer without an initial guess.

## What changed and what stayed

This approximation omits light-travel time, the gravity ODE, and the bound-orbit constraint. The hard solution direction has 0.47 arcsec RMS and 2.04 arcsec maximum pointwise separation in the current fit, so the roughly 1 arcsec statement is an RMS-scale check rather than a per-observation maximum. It retains ERFA Earth position and the original observatory offset. The easy state is therefore a polynomial fit, while the hard engine propagates a physical two-body orbit and evaluates light-time.

## Measured price of realism

Computed by python -m easy.run_all from 343 observations; all values come from the current run. See results/summary.md and results/experiment.json for the full degree sweep. In the latest run, Easy B's own ray residual is about 9.98 arcsec/axis, while its state scores about 1,841 arcsec/axis when transplanted into the original gravity and light-time engine. The original hard fit is about 0.333 arcsec/axis. This is a large and revealing mismatch: a local degree-2 polynomial does not capture the physical trajectory well enough across this interval.

## Likely questions

- Why project with I-uuᵀ? It removes the component parallel to the measured ray, leaving only sideways error.
- Why is the solve global? It is least squares on a linear map, so the objective is a convex quadratic.
- Why no initial guess? lstsq solves the linear problem directly.
- Why does transplant error grow? The polynomial model does not satisfy the gravity ODE or light-time model.
- What does the Hessian show? Its eigenvalues are nonnegative up to floating-point roundoff; it certifies convexity, not model accuracy.
- What are the limitations? This model has no range observation, no physical dynamics, and no constraint against implausible states. Higher degree can fit the observed arc while giving unstable extrapolation or state estimates.

## Honest interpretation

The computed degree-2 transplant error is very large. That is not a tuning target; it is the cost observed for the requested simplified model and this fit window. Degree and weighting change the numerical approximation but do not restore the omitted physical assumptions.
