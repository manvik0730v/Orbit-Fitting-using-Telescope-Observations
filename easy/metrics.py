"""Residual metrics and independent two-body orbital elements."""
import numpy as np
from b_constants import GM_SUN, ARCSEC
from d_relativePosition import result as hard_engine
from c_Position import simulate
from .geometry import evaluate, coefficients_to_state
def angular_residual(ra,dec,vec):
    """Return the original engine's wrapped RA*cos(dec), Dec residuals in arcsec."""
    pred_ra=np.remainder(np.arctan2(vec[:,1],vec[:,0]),2*np.pi)
    pred_dec=np.arctan2(vec[:,2],np.hypot(vec[:,0],vec[:,1]))
    dra=np.remainder(ra-pred_ra+np.pi,2*np.pi)-np.pi
    return np.stack((dra*np.cos(dec),dec-pred_dec),axis=1)*ARCSEC
def own_rmse(data,tau,c,degree):
    """Compute angular residuals from the fitted polynomial ray vectors."""
    p=evaluate(c,tau,degree); vec=p-data['R']
    return float(np.sqrt(np.mean(angular_residual(data['ra'],data['dec'],vec)**2))), angular_residual(data['ra'],data['dec'],vec)
def transplant_rmse(data,r0,v0,t0):
    """Score an easy state with the original full light-time engine."""
    ra,dec=hard_engine(data['jd_utc'],r0,v0,t0,data['d'],data['h'],data['L'])
    dra=np.remainder(data['ra']-ra+np.pi,2*np.pi)-np.pi
    e=np.r_[dra*np.cos(data['dec']),data['dec']-dec]*ARCSEC
    return float(np.sqrt(np.mean(e**2))),e.reshape(2,-1).T
def elements(r,v):
    """Compute standard osculating elements in the project's equatorial frame."""
    r=np.asarray(r); v=np.asarray(v); rn=np.linalg.norm(r); vn=np.linalg.norm(v); mu=GM_SUN
    h=np.cross(r,v); hn=np.linalg.norm(h); n=np.cross([0.,0.,1.],h); nn=np.linalg.norm(n)
    ev=np.cross(v,h)/mu-r/rn; e=np.linalg.norm(ev); energy=vn**2/2-mu/rn; a=-mu/(2*energy)
    inc=np.arccos(np.clip(h[2]/hn,-1,1))
    Om=np.arctan2(n[1],n[0])%(2*np.pi)
    argp=np.arctan2(np.dot(np.cross(n,ev),h)/(nn*e*hn),np.dot(n,ev)/(nn*e))%(2*np.pi) if e>1e-12 and nn>1e-12 else 0.
    nu=np.arctan2(np.dot(np.cross(ev,r),h)/(e*rn*hn),np.dot(ev,r)/(e*rn)) if e>1e-12 else 0.
    E=2*np.arctan2(np.sqrt(1-e)*np.sin(nu/2),np.sqrt(1+e)*np.cos(nu/2))
    M=E-e*np.sin(E); T_days=-M/np.sqrt(mu/a**3); q=a*(1-e); Q=a*(1+e)
    return {'a_au':float(a),'e':float(e),'i_deg':float(np.degrees(inc)),'T_days_from_epoch':float(T_days),'q_au':float(q),'Q_au':float(Q),'raan_deg':float(np.degrees(Om)),'argp_deg':float(np.degrees(argp))}
