"""Ray projection geometry for measured telescope directions."""
import numpy as np
from math import factorial
def unit_directions(ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    """Convert equatorial RA/Dec to unit Cartesian sight lines."""
    return np.stack((np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)), axis=-1)
def projectors(u: np.ndarray) -> np.ndarray:
    """Project vectors onto the plane perpendicular to each ray."""
    return np.eye(3)[None, :, :] - u[:, :, None] * u[:, None, :]
def design_matrix(tau: np.ndarray, u: np.ndarray, R: np.ndarray, degree: int, scale: float=30.0):
    """Build scaled linear system for perpendicular position residuals."""
    P=projectors(u); x=tau/scale; n=len(tau); k=degree+1
    blocks=[]
    for q in range(k):
        fac=1.0 / float(factorial(q))
        blocks.append(P * (x**q*fac)[:,None,None])
    A=np.concatenate(blocks,axis=2).reshape(3*n,3*k)
    b=(P@R[:,:,None]).reshape(3*n)
    return A,b
def coefficients_to_state(c: np.ndarray, degree: int, scale: float=30.0):
    """Undo scaled-time coefficients; return physical r0 and v0."""
    z=c.reshape(degree+1,3)
    return z[0].copy(), z[1].copy()/scale
def evaluate(c: np.ndarray, tau: np.ndarray, degree: int, scale: float=30.0):
    """Evaluate polynomial path in physical units."""
    z=c.reshape(degree+1,3); x=tau/scale
    return sum(z[q]*x[:,None]**q/float(factorial(q)) for q in range(degree+1))
