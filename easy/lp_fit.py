"""Absolute-residual linear programming fit."""
import numpy as np
from scipy.optimize import linprog
from .geometry import design_matrix
def fit_l1(tau: np.ndarray,u: np.ndarray,R: np.ndarray,degree: int=2):
    """Minimize sum of absolute stacked perpendicular components as an LP."""
    A,b=design_matrix(tau,u,R,degree); m,n=A.shape
    c=np.r_[np.zeros(n),np.ones(m)]
    aub=np.block([[A,-np.eye(m)],[-A,-np.eye(m)]])
    bub=np.r_[b,-b]
    res=linprog(c,A_ub=aub,b_ub=bub,bounds=[(None,None)]*n+[(0,None)]*m,method='highs')
    if not res.success: raise RuntimeError(res.message)
    return res.x[:n],float(res.fun),res
