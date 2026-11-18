# Slide text

## Slide 1: The idea: a telescope measures a ray

- A telescope reports the direction to Eros.
- It does not report how far away Eros is.
- We compare the predicted position to the whole ray.

Figure: easy/figures/ray_geometry.png

## Slide 2: Measure the miss as a distance

- uᵢ = (cos δ cos α, cos δ sin α, sin δ)
- Pᵢ = I − uᵢuᵢᵀ removes the along-ray part.
- ||Pᵢ(pᵢ − Rᵢ)|| is the sideways distance in au.

Figure: easy/figures/ray_geometry.png

## Slide 3: Model Eros’s path as a polynomial

- A: degree 1, 6 unknowns (position and velocity)
- B: degree 2, 9 unknowns (headline model)
- C: degree 3, 12 unknowns
- Degrees 1–5 were also evaluated.

Figure: easy/figures/degree_sweep.png

## Slide 4: Easy problem formulation

- minimize Σᵢ ||Pᵢ(pᵢ − Rᵢ)||²
- pᵢ = r₀ + v₀τᵢ + ½aτᵢ² + ⅙jτᵢ³ + …
- Dropped: light-time, gravity ODE, bound-orbit constraint.
- Kept: pyerfa Earth position and telescope offset Rᵢ=r₁+r₂.

Figure: easy/figures/ray_geometry.png

## Slide 5: Why it is convex, with no initial guess

- The unknown coefficients enter linearly.
- The objective is a sum of squared linear residuals.
- A linear least-squares solve gives the global minimum.
- There is no starting guess to choose.

Figure: easy/figures/convexity_line.png

## Slide 6: Easy results: the full-engine check matters

- Easy B own RMSE: 9.984 arcsec/axis
- Easy B transplanted RMSE: 1841.0 arcsec/axis
- Hard fit RMSE: 0.333 arcsec/axis
- Polynomial accuracy does not guarantee a physical orbit.

Figure: easy/figures/degree_sweep.png

## Slide 7: Price of Realism II

- Easy B assumes a local polynomial; hard fit enforces gravity and light-time.
- State offset: |Δr₀| = 0.278 au, |Δv₀| = 0.0017 au/day
- Elements: Easy B a=1.694 au, e=0.299; hard a=1.458 au, e=0.223.
- Transplant RMSE: Easy B 1841.0 vs hard 0.333 arcsec/axis.

Figure: easy/figures/price_of_realism_2.png
