from pathlib import Path
import importlib.util,json,math,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/mpl_v9')
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
HERE=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('L',HERE/'shared/lfv_letter.py');L=importlib.util.module_from_spec(spec);spec.loader.exec_module(L)
import sys
sys.path.insert(0,str(HERE/'shared/physical_lfv'))
import physical_threebody as P
HH=np.load(HERE/'shared/physical_lfv/gram_320.npy')
plt.rcParams.update({'font.size':10,'font.family':'serif','mathtext.fontset':'cm','axes.spines.top':False,'axes.spines.right':False})
y=L.yukawa(L.YE);fig,ax=plt.subplots(1,2,figsize=(7.2,3.1))
M=np.linspace(2.5,12,80);B=[];C=[]
for m in M:
 cc=P.coefficients(P.fields(m*1000),np.full(8,m*1000));B.append(P.radiative(cc));C.append(P.rate(cc,HH))
ax[0].semilogy(M,B,label=r'$e\gamma$ dipole',color='#15577b');ax[0].semilogy(M,C,label=r'$3e$ tree + finite photon',color='#b54b25')
for lev,col,ls in [(1.5e-13,'#15577b','--'),(6e-14,'#15577b',':'),(4.60517e-15,'#b54b25','--'),(2e-15,'#b54b25',':'),(1e-16,'#555555','-.')]:ax[0].axhline(lev,color=col,lw=.7,ls=ls)
ax[0].set(xlabel=r'$M=\Lambda_D$ [TeV]',ylabel='branching ratio',ylim=(1e-19,1e-11));ax[0].legend(frameon=False,fontsize=8)
_,cp=L.couplings(y,5000);ks=np.geomspace(.01,1,65);bg=[];bc=[];bp=[]
for k in ks:
 cp2={r:(s*(k if r==1 else 1),a*(k if r==1 else 1)) for r,(s,a) in cp.items()};_,_,b=L.dipole(cp2,5000);_,c=L.contact(cp2,5000);bg.append(b);bc.append(c);bp.append(b*L.ALPHA/(3*math.pi)*(math.log(L.MMU**2/L.ME**2)-11/4))
ax[1].loglog(ks,bg,label=r'$e\gamma$',color='#15577b');ax[1].loglog(ks,bc,label=r'$3e$ contact',color='#b54b25');ax[1].loglog(ks,bp,label=r'$3e$ dipole only',color='#555555',ls=':');ax[1].set(xlabel=r'$\kappa_1$ at $M=5$ TeV',ylabel='branching ratio');ax[1].legend(frameon=False,fontsize=8)
for a,t in zip(ax,['(a)','(b)']):a.text(.04,.96,t,transform=a.transAxes,va='top')
fig.tight_layout();fig.savefig(HERE/'figures/lfv_5TeV_alternative.pdf');plt.close(fig)
ll=np.linspace(2,50,100);b=[];c=[]
for lam in ll:
 _,cp=L.couplings(y,lam*1000);_,_,bb=L.dipole(cp,1000);_,cc=L.contact(cp,1000);b.append(bb);c.append(cc)
fig,ax=plt.subplots(figsize=(6.1,3.1));ax.semilogy(ll,b,label=r'$e\gamma$ dipole',color='#15577b');ax.semilogy(ll,c,label=r'$3e$ contact',color='#b54b25');ax.axhline(1.5e-13,color='#15577b',ls='--');ax.axhline(1e-12,color='#b54b25',ls='--');ax.set(xlabel=r'$\Lambda_\ell$ [TeV], fixed 1-TeV poles',ylabel='branching ratio',ylim=(1e-17,1e-8));ax.legend(frameon=False);fig.tight_layout();fig.savefig(HERE/'figures/lfv_scale_suppression.pdf');plt.close(fig)
# Parent threshold running; no obsolete chain-site curves or inferred walking band.
J=json.loads((HERE/'shared/joint_consistency_v9.json').read_text());d=J['formal_links'][2]
fig,ax=plt.subplots(1,2,figsize=(7.2,3.1));ds=np.linspace(1.4,3,160);X=3/(4*math.pi**2)/7.18e-10
ax[0].semilogy(ds,5000*X**(1/(ds-1)),color='#15577b',label='Conditional link')
ax[0].axhline(d['LambdaY_formal_GeV'],color='#b54b25',ls='--',label=r'Formal $\Lambda_Y$')
for dm,col,lab in [(2,'#666666','Walking example'),(3,'#15577b','Canonical choice')]:ax[0].plot(dm,5000*X**(1/(dm-1)),'o',color=col,label=lab)
ax[0].set(xlabel=r'$\Delta_\nu$',ylabel='scale [GeV]',ylim=(1e6,1e25),xlim=(1.4,3.05));ax[0].legend(frameon=False,fontsize=8)
mu=np.geomspace(91.1876,3e16,300);aY=(1/127.951)/(1-.23122)
for N,ND,col in [(2,3,'#15577b'),(3,4,'#b54b25')]:
 bs=4*N/3+20*ND/9+2*N*ND/27
 inv=np.array([1/aY-41/6/(2*math.pi)*math.log(m/91.1876)-sum(db/(2*math.pi)*math.log(m/t) for t,db in [(1749.5,16/9),(6398.6,2),(6398.6,2),(5000,bs)] if m>t) for m in mu])
 ax[1].semilogx(mu,inv,color=col,label=f'Ranks ({N},{N},{ND})')
ax[1].axhline(0,color='#555555',ls=':');ax[1].set(xlabel=r'$\mu$ [GeV]',ylabel=r'$\alpha_Y^{-1}$',ylim=(-15,110));ax[1].legend(frameon=False,fontsize=8)
for a,t in zip(ax,['(a)','(b)']):a.text(.04,.96,t,transform=a.transAxes,va='top')
fig.tight_layout();fig.savefig(HERE/'figures/neutrino_link_consistency.pdf');plt.close(fig)
print('Three flavour vector figures regenerated.')
