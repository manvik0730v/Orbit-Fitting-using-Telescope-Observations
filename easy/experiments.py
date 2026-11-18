"""Run easy polynomial fits, hard baseline, and write reproducible reports."""
import json, os, time
import numpy as np
from scipy.optimize import least_squares
from _process_data import get_arrays
from d_relativePosition import R1,R2
from f_method1 import miss
from c_Position import simulate
from .config import FIT_START,FIT_END,t0,DEGREES
from .data import fit_arrays
from .geometry import unit_directions
from .convex_fit import fit,state
from .metrics import own_rmse,transplant_rmse,elements
def run():
    """Run degree sweep and hard comparison; save JSON and summary."""
    d=fit_arrays(); n=len(d['ra']); tau=d['jd_utc']+69.184/86400-t0
    u=unit_directions(d['ra'],d['dec']); d['R']=R1(d['jd_utc'])+R2(d['jd_utc'],d['d'],d['h'],d['L'])
    rows=[]; fits={}
    for degree in DEGREES:
        f=fit(tau,u,d['R'],degree); fits[degree]=f; r,v=state(f)
        own,_=own_rmse(d,tau,f.coefficients,degree); trans,_=transplant_rmse(d,r,v,t0)
        eig=np.linalg.eigvalsh(2*f.A.T@f.A)
        rows.append({'degree':int(degree),'own_rmse_arcsec_axis':own,'transplant_rmse_arcsec_axis':trans,'r0_au':r.tolist(),'v0_au_day':v.tolist(),'elements':elements(r,v),'solve_seconds':f.seconds,'cond_AtA':float(np.linalg.cond(f.A.T@f.A)),'min_hessian_eigenvalue':float(eig[0]),'rank':int(f.rank)})
    # Compute true original hard fit in-memory, using its prescribed Gauss initialization.
    import f_method1 as hard
    y0=min(hard.gauss(0,n//2,n-1),key=lambda y:np.sum(miss(y)**2))
    tic=time.perf_counter(); sol=least_squares(miss,y0,x_scale='jac'); hard_seconds=time.perf_counter()-tic
    hr,hv=sol.x[:3],sol.x[3:]
    hown=float(np.sqrt(np.mean(sol.fun**2)))
    hjson={'r0_au':hr.tolist(),'v0_au_day':hv.tolist(),'rmse_arcsec_axis':hown,'seconds':hard_seconds,'elements':elements(hr,hv)}
    for row in rows:
        row['delta_r0_easy_minus_hard_au']=(np.array(row['r0_au'])-hr).tolist()
        row['delta_v0_easy_minus_hard_au_day']=(np.array(row['v0_au_day'])-hv).tolist()
    out={'n':n,'t0_jd_tt':t0,'hard':hjson,'degrees':rows}
    os.makedirs('easy/results',exist_ok=True)
    json.dump(out,open('easy/results/experiment.json','w'),indent=2)
    lines=['# Easy convex fit results','',f"Fit observations: {n}",'', '| Degree | Own RMSE (arcsec/axis) | Transplant RMSE | cond(AᵀA) | min eig(2AᵀA) | solve s |','|---:|---:|---:|---:|---:|---:|']
    for x in rows: lines.append(f"| {x['degree']} | {x['own_rmse_arcsec_axis']:.6g} | {x['transplant_rmse_arcsec_axis']:.6g} | {x['cond_AtA']:.4g} | {x['min_hessian_eigenvalue']:.4g} | {x['solve_seconds']:.4g} |")
    lines+=['',f"Hard fit RMSE: {hown:.8g} arcsec/axis; runtime {hard_seconds:.4g} s."]
    open('easy/results/summary.md','w').write('\n'.join(lines)+'\n')
    return out,d,tau,u,fits,sol
if __name__=='__main__': run()
