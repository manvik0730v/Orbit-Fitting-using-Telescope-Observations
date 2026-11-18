"""Generate clean, white-background figures for the analysis."""
import os,json
import numpy as np
import matplotlib.pyplot as plt
from .experiments import run
from .geometry import evaluate
from .metrics import angular_residual
def make_plots():
    """Write all requested diagnostic and teaching figures."""
    os.makedirs('easy/figures',exist_ok=True); out,d,tau,u,fits,sol=run(); p='easy/figures/'
    plt.style.use('seaborn-v0_8-whitegrid')
    fig,ax=plt.subplots(figsize=(7,4)); ax.plot([0,4],[0,2],lw=3,color='#3076b8',label='measured ray'); ax.scatter([2.3],[3.5],s=80,color='#d95f02',label='predicted point'); ax.scatter([2.3],[1.15],s=75,color='#2a9d65',label='closest point on ray'); ax.annotate('',xy=(2.3,3.4),xytext=(2.3,1.25),arrowprops=dict(arrowstyle='<->',color='#b32d2e',lw=2)); ax.text(2.45,2.2,'sideways miss'); ax.text(3.05,1.8,'along-ray part',rotation=26); ax.legend(); ax.set(xlim=(-.2,4.3),ylim=(-.2,4),xlabel='x',ylabel='y',title='A telescope measures a ray, not a distance'); fig.tight_layout(); fig.savefig(p+'ray_geometry.png',dpi=180); plt.close(fig)
    # Convex quadratic along a segment, with its chord below the curve.
    f=fits[2]; c0=f.coefficients; rng=np.random.default_rng(2026); c1=c0+rng.normal(size=len(c0))*np.maximum(np.abs(c0),.02)*.1
    A=f.A; b=A@c0-f.residual; ss=np.linspace(0,1,200); vals=np.array([np.linalg.norm(A@(c0+s*(c1-c0))-b)**2 for s in ss]); chord=np.linspace(vals[0],vals[-1],len(ss))
    fig,ax=plt.subplots(figsize=(7,4)); ax.plot(ss,vals,lw=2.5,label='objective along segment'); ax.plot(ss,chord,'--',lw=2,label='chord'); ax.set(xlabel='fraction along coefficient segment',ylabel='sum of squared misses (au²)',title='A convex objective stays below its chord'); ax.legend(); fig.tight_layout(); fig.savefig(p+'convexity_line.png',dpi=180); plt.close(fig)
    deg=[x['degree'] for x in out['degrees']]
    fig,ax=plt.subplots(figsize=(7,4)); ax.plot(deg,[x['own_rmse_arcsec_axis'] for x in out['degrees']],'o-',label='own ray residual'); ax.plot(deg,[x['transplant_rmse_arcsec_axis'] for x in out['degrees']],'o-',label='full-engine transplant'); ax.axhline(out['hard']['rmse_arcsec_axis'],ls='--',color='black',label='hard fit'); ax.set_yscale('log'); ax.set(xlabel='polynomial degree',ylabel='RMSE (arcsec per axis)',title='Degree sweep'); ax.legend(); fig.tight_layout(); fig.savefig(p+'degree_sweep.png',dpi=180); plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(9,4)); names=['Hard','Easy A','Easy B','Easy C']; hs=[out['hard']['r0_au'],* [out['degrees'][i]['r0_au'] for i in range(3)]]; vs=[out['hard']['v0_au_day'],* [out['degrees'][i]['v0_au_day'] for i in range(3)]]
    for ax,vals,title,unit in [(axs[0],hs,'Position at t0 (au)','au'),(axs[1],vs,'Velocity at t0 (au/day)','au/day')]:
        x=np.arange(4); ax.bar(x,[np.linalg.norm(np.array(v)-hs[0] if title.startswith('Position') else np.array(v)-vs[0]) for v in vals],color=['#414a6b','#55a6a6','#e7ae4a','#ca6f62']); ax.set_xticks(x,names); ax.set_title(title+'; bar = distance from hard state'); ax.set_ylabel(unit)
    fig.tight_layout(); fig.savefig(p+'price_of_realism_2.png',dpi=180); plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4)); idx=np.arange(len(tau)); easy=evaluate(fits[2].coefficients,tau,2)-d['R']; easyres=angular_residual(d['ra'],d['dec'],easy); hardres=sol.fun
    ax.plot(idx,np.linalg.norm(easyres,axis=1),label='Easy B'); ax.plot(idx,np.linalg.norm(hardres.reshape(2,-1).T,axis=1),label='Hard fit'); ax.set(xlabel='observation index (time order)',ylabel='2-axis residual norm (arcsec)',title='Residuals across the fit window'); ax.legend(); fig.tight_layout(); fig.savefig(p+'residuals_vs_time.png',dpi=180); plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,4)); ax.axis('off'); keys=[('a (au)','a_au'),('e','e'),('i (deg)','i_deg'),('T − t0 (days)','T_days_from_epoch'),('q (au)','q_au'),('Q (au)','Q_au')]; h=out['hard']['elements']; b=out['degrees'][1]['elements']; cells=[[f'{h[k]:.6g}',f'{b[k]:.6g}'] for label,k in keys]; table=ax.table(cellText=cells,rowLabels=[a for a,k in keys],colLabels=['Hard fit','Easy B'],loc='center',cellLoc='center'); table.auto_set_font_size(False); table.set_fontsize(11); table.scale(1,1.65); ax.set_title('Osculating orbital elements at t0',pad=18); fig.tight_layout(); fig.savefig(p+'elements_comparison.png',dpi=180); plt.close(fig)
if __name__=='__main__': make_plots()
