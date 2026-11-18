"""Create seven dark-theme slides and merge them into the supplied deck."""
import os
import numpy as np
import textwrap
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pypdf import PdfReader,PdfWriter
from .experiments import run
BG='#303234'; WHITE='#f5f4f1'; ORANGE='#ff7900'
def _base(title,figpath=None):
    fig=plt.figure(figsize=(16,9),facecolor=BG); ax=fig.add_axes([0,0,1,1]); ax.set_axis_off()
    ax.add_patch(Rectangle((.015,.025),.97,.95,fill=False,edgecolor='#c9c8c3',lw=2,transform=ax.transAxes))
    
    ax.plot([.08,.92],[.84,.84],color='#bdbdb8',lw=1.2,transform=ax.transAxes)
    if figpath:
        img=plt.imread(figpath); ia=fig.add_axes([.52,.16,.40,.59]); ia.imshow(img); ia.axis('off')
    fig.text(.5,.91,title,ha='center',va='center',color=WHITE,fontsize=30,fontweight='bold',fontfamily='DejaVu Serif')
    return fig,ax
def _bullets(ax,lines,dy=.10,size=17):
    y=.73
    for line in lines:
        ax.text(.09,y,'□',color='#ff7900',fontsize=size+1,va='center',transform=ax.transAxes)
        wrapped='\n'.join(textwrap.wrap(line,width=42))
        ax.text(.12,y,wrapped,color=WHITE,fontsize=size,fontweight='bold',va='center',fontfamily='DejaVu Sans',transform=ax.transAxes)
        y-=dy * max(1,wrapped.count('\n')+1)
def create():
    """Write slide_text.md and the merged presentation PDF."""
    out,d,tau,u,fits,hard=run(); os.makedirs('easy/slides',exist_ok=True); texts=[]
    specs=[
      ('The idea: a telescope measures a ray',['A telescope reports the direction to Eros.','It does not report how far away Eros is.','We compare the predicted position to the whole ray.'],'ray_geometry.png'),
      ('Measure the miss as a distance',['uᵢ = (cos δ cos α, cos δ sin α, sin δ)','Pᵢ = I − uᵢuᵢᵀ removes the along-ray part.','||Pᵢ(pᵢ − Rᵢ)|| is the sideways distance in au.'],'ray_geometry.png'),
      ('Model Eros’s path as a polynomial',['A: degree 1, 6 unknowns (position and velocity)','B: degree 2, 9 unknowns (headline model)','C: degree 3, 12 unknowns','Degrees 1–5 were also evaluated.'],'degree_sweep.png'),
      ('Easy problem formulation',['minimize Σᵢ ||Pᵢ(pᵢ − Rᵢ)||²','pᵢ = r₀ + v₀τᵢ + ½aτᵢ² + ⅙jτᵢ³ + …','Dropped: light-time, gravity ODE, bound-orbit constraint.','Kept: pyerfa Earth position and telescope offset Rᵢ=r₁+r₂.'],'ray_geometry.png'),
      ('Why it is convex, with no initial guess',['The unknown coefficients enter linearly.','The objective is a sum of squared linear residuals.','A linear least-squares solve gives the global minimum.','There is no starting guess to choose.'],'convexity_line.png'),
      ('Easy results: the full-engine check matters',['Easy B own RMSE: %.3f arcsec/axis'%out['degrees'][1]['own_rmse_arcsec_axis'],'Easy B transplanted RMSE: %.1f arcsec/axis'%out['degrees'][1]['transplant_rmse_arcsec_axis'],'Hard fit RMSE: %.3f arcsec/axis'%out['hard']['rmse_arcsec_axis'],'Polynomial accuracy does not guarantee a physical orbit.'],'degree_sweep.png')]
    for idx,(title,bullets,img) in enumerate(specs,1):
        fig,ax=_base(title,'easy/figures/'+img); _bullets(ax,bullets,dy=.095,size=16)
        fig.savefig('/private/tmp/easy-slide-%02d.png'%idx,dpi=120,facecolor=BG); plt.close(fig)
        texts.append('## Slide %d: %s\n\n'%(idx,title)+'\n'.join('- '+b for b in bullets)+'\n\nFigure: easy/figures/'+img)
    e=out['degrees']; he=out['hard']['elements']; be=e[1]['elements']
    bullets=['Easy B assumes a local polynomial; hard fit enforces gravity and light-time.','State offset: |Δr₀| = %.3f au, |Δv₀| = %.4f au/day'%(np.linalg.norm(np.array(e[1]['r0_au'])-out['hard']['r0_au']),np.linalg.norm(np.array(e[1]['v0_au_day'])-out['hard']['v0_au_day'])),'Elements: Easy B a=%.3f au, e=%.3f; hard a=%.3f au, e=%.3f.'%(be['a_au'],be['e'],he['a_au'],he['e']),'Transplant RMSE: Easy B %.1f vs hard %.3f arcsec/axis.'%(e[1]['transplant_rmse_arcsec_axis'],out['hard']['rmse_arcsec_axis'])]
    fig,ax=_base('Price of Realism II: easy vs hard','easy/figures/price_of_realism_2.png'); _bullets(ax,bullets,dy=.08,size=14)
    fig.savefig('/private/tmp/easy-slide-07.png',dpi=120,facecolor=BG);plt.close(fig)
    texts.append('## Slide 7: Price of Realism II\n\n'+'\n'.join('- '+b for b in bullets)+'\n\nFigure: easy/figures/price_of_realism_2.png')
    open('easy/slides/slide_text.md','w').write('# Slide text\n\n'+'\n\n'.join(texts)+'\n')
    from reportlab.pdfgen import canvas
    temp='/private/tmp/easy-new-slides.pdf'; c=canvas.Canvas(temp,pagesize=(1920,1080))
    for i in range(1,8): c.drawImage('/private/tmp/easy-slide-%02d.png'%i,0,0,width=1920,height=1080);c.showPage()
    c.save()
    original=PdfReader('easy/slides/original_deck.pdf'); new=PdfReader(temp); writer=PdfWriter()
    for i,p in enumerate(original.pages):
        writer.add_page(p)
        if i==7:
            for slide in new.pages[:6]: writer.add_page(slide)
    writer.add_page(new.pages[6])
    with open('easy/slides/Orbit_Fitting_with_easy.pdf','wb') as f: writer.write(f)
if __name__=='__main__': create()
