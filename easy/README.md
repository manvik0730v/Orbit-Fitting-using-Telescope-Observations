# Easy convex orbit fit

Run from the repository root with python -m easy.run_all. The implementation imports the project observation pipeline, UTC-to-TT convention, ERFA Earth position, station offset, and hard-engine scoring. The default fit window contains 343 observations from 2025-10-15 through 2025-12-15.

The solve minimizes squared perpendicular distances to the measured rays. It uses a time-scaled polynomial basis and numpy.linalg.lstsq; reported state vectors are converted back to au and au/day. Results, plots, and the merged deck are written under this directory.

The transplant metric sends each fitted (r0, v0) into the full original engine, including gravity integration and light-time iteration. It is expected to differ from the local polynomial residual because the polynomial does not enforce physical dynamics.

Dependencies are listed in requirements-easy.txt. Run checks with python -m unittest easy.tests.test_easy -v.
