"""Numerical checks for ray geometry, convex fits, and the hard baseline."""
import unittest
import numpy as np
from scipy.linalg import solve
from _process_data import get_arrays
from d_relativePosition import R1,R2,R5
from f_method1 import miss,t0,gauss
from scipy.optimize import least_squares
from easy.geometry import unit_directions,projectors
from easy.convex_fit import fit
from easy.lp_fit import fit_l1
from easy.config import FIT_START,FIT_END
class EasyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=get_arrays(FIT_START,FIT_END); d=cls.d
        cls.tau=d['jd_utc']+69.184/86400-t0
        cls.u=unit_directions(d['ra'],d['dec']); cls.R=R1(d['jd_utc'])+R2(d['jd_utc'],d['d'],d['h'],d['L'])
    def test_geometry_and_station(self):
        self.assertTrue(np.allclose(np.linalg.norm(self.u,axis=1),1))
        self.assertTrue(np.allclose(np.linalg.norm(self.R,axis=1),1,atol=.02))
        for i in [0,100,342]: self.assertTrue(np.allclose(self.R[i],R1(self.d['jd_utc'][i])+R2(self.d['jd_utc'][i:i+1],self.d['d'][i:i+1],self.d['h'][i:i+1],self.d['L'][i:i+1])[0]))
        P=projectors(self.u)
        self.assertTrue(np.allclose(P@P,P,atol=1e-12))
        self.assertTrue(np.allclose((P@self.u[:,:,None])[:,:,0],0,atol=1e-12))
    def test_convex_solution_and_hessian(self):
        f=fit(self.tau,self.u,self.R,2); H=2*f.A.T@f.A
        self.assertGreaterEqual(np.linalg.eigvalsh(H).min(),-1e-9)
        b=np.linalg.lstsq(f.A, f.A@f.coefficients-f.residual,rcond=None)[0]
        closed=solve(f.A.T@f.A,f.A.T@(f.A@f.coefficients-f.residual),assume_a='pos')
        self.assertTrue(np.allclose(f.coefficients,closed,rtol=1e-6,atol=1e-7))
    def test_hard_direction(self):
        d=self.d; n=len(d['ra']); y0=min(gauss(0,n//2,n-1),key=lambda y:np.sum(miss(y)**2))
        y=least_squares(miss,y0,x_scale='jac').x
        v=R5(d['jd_utc'],y[:3],y[3:],t0,d['d'],d['h'],d['L'])
        pred=v/np.linalg.norm(v,axis=1)[:,None]
        sep=np.arccos(np.clip(np.sum(pred*self.u,axis=1),-1,1))*206264.806
        self.assertLess(float(np.sqrt(np.mean(sep**2))),1.0)
    def test_lp_l1_dominates_ls(self):
        f=fit(self.tau,self.u,self.R,2)
        coeff,obj,_=fit_l1(self.tau,self.u,self.R,2)
        ls_l1=float(np.abs(f.residual).sum())
        self.assertLessEqual(obj,ls_l1+1e-6)
if __name__=='__main__': unittest.main(verbosity=2)
