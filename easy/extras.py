"""Optional angular weighting, L1 fit, and numerical convexity check."""
import json
import numpy as np
from .convex_fit import fit,state
from .geometry import evaluate
from .metrics import own_rmse,transplant_rmse
from .lp_fit import fit_l1
def run_extras(out,d,tau,u,fits):
    """Compute requested optional variants and append them to experiment JSON."""
    base=fits[2]; p=evaluate(base.coefficients,tau,2); rho=np.maximum(np.linalg.norm(p-d['R'],axis=1),1e-10)
    weighted=[]
    for _ in range(2):
        f=fit(tau,u,d['R'],2,weights=1/rho); rho=np.maximum(np.linalg.norm(evaluate(f.coefficients,tau,2)-d['R'],axis=1),1e-10)
        r,v=state(f); own,_=own_rmse(d,tau,f.coefficients,2); trans,_=transplant_rmse(d,r,v,out['t0_jd_tt'])
        weighted.append({'own_rmse_arcsec_axis':own,'transplant_rmse_arcsec_axis':trans,'r0_au':r.tolist(),'v0_au_day':v.tolist()})
    c,obj,res=fit_l1(tau,u,d['R'],2); ls_l1=float(np.abs(base.residual).sum())
    rng=np.random.default_rng(433); A=base.A; b=A@base.coefficients-base.residual
    def fval(x): return float(np.linalg.norm(A@x-b)**2)
    violations=0
    for _ in range(1000):
        x=rng.normal(size=len(c)); y=rng.normal(size=len(c)); mid=(x+y)/2
        if fval(mid) > (fval(x)+fval(y))/2+1e-9*max(fval(x),fval(y),1): violations+=1
    out['weighted_two_pass']=weighted
    out['l1_degree2']={'objective':obj,'least_squares_l1_objective':ls_l1,'lp_leq_ls':obj<=ls_l1+1e-6}
    out['convexity_midpoint_tests']={'n':1000,'violations':violations}
    json.dump(out,open('easy/results/experiment.json','w'),indent=2)
    with open('easy/results/summary.md','a') as f:
        f.write(f"\n\nAngular weighting (two passes): own RMSE {weighted[-1]['own_rmse_arcsec_axis']:.6g}, transplant RMSE {weighted[-1]['transplant_rmse_arcsec_axis']:.6g} arcsec/axis.\n")
        f.write(f"Degree-2 LP L1 objective: {obj:.8g}; LS L1 objective: {ls_l1:.8g}; LP <= LS: {obj<=ls_l1+1e-6}.\n")
        f.write(f"Random midpoint tests: {violations} violations out of 1000.\n")
    return out
