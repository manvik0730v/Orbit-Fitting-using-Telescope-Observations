"""Linear least-squares polynomial ray fit."""
from dataclasses import dataclass
import time
import numpy as np
from .geometry import design_matrix, coefficients_to_state
@dataclass
class Fit:
    """Fitted polynomial coefficients and numerical diagnostics."""
    degree: int
    coefficients: np.ndarray
    residual: np.ndarray
    seconds: float
    rank: int
    singular_values: np.ndarray
    A: np.ndarray
def fit(tau: np.ndarray, u: np.ndarray, R: np.ndarray, degree: int, weights: np.ndarray|None=None) -> Fit:
    """Solve the unconstrained convex ray-distance least-squares problem."""
    A,b=design_matrix(tau,u,R,degree)
    if weights is not None:
        w=np.repeat(np.sqrt(weights),3); A=A*w[:,None]; b=b*w
    start=time.perf_counter(); c,_,rank,s=np.linalg.lstsq(A,b,rcond=None); elapsed=time.perf_counter()-start
    return Fit(degree,c,A@c-b,elapsed,rank,s,A)
def state(result: Fit):
    """Extract physical position and velocity at the reference epoch."""
    return coefficients_to_state(result.coefficients,result.degree)
