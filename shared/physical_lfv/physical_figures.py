from pathlib import Path
import numpy as np,json,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import physical_threebody as p
R=Path(__file__).parent;H=np.load(R/'gram_320.npy');out=json.loads((R/'physical_results_v10.json').read_text());fig=R.parent/'figures'
plt.rcParams.update({'font.size':9,'axes.labelsize':10,'legend.fontsize':8,'savefig.bbox':'tight'})
d=np.linspace(-.012,.012,801);Q=p.fields();c=[p.coefficients(Q,5000*np.tile([1+x,1-x],4)) for x in d];B=np.array([p.radiative(v) for v in c]);r=np.array([p.rate(v,H)/p.radiative(v) for v in c]);rd=out['RD']
f,ax=plt.subplots(2,1,figsize=(3.45,3.9),sharex=True)
ax[0].semilogy(100*d,r,color='#124e78');ax[0].axhline(rd,color='gray',ls='--',label='Dipole only');ax[0].axhline(10*rd,color='#bc6b23',ls=':',label=r'$10R_D$');ax[0].set_ylabel(r'$B_{3e}/B_{e\gamma}$');ax[0].legend(loc='lower left')
ax[1].semilogy(100*d,B,color='#124e78');ax[1].axhline(1.5e-13,color='#a63c33',ls='--',label='MEG II 90% limit');ax[1].set_xlabel(r'Common pole deformation $100\delta$ (%)');ax[1].set_ylabel(r'$B(\mu\to e\gamma)$');ax[1].legend(loc='lower right');f.subplots_adjust(hspace=.13);f.savefig(fig/'fig_pole_alignment.pdf');plt.close(f)
# Common-scale curves: baseline versus split counterexample.
m=np.linspace(3,11,180);rate=[];split=[]
for M in m:
 cc=p.coefficients(p.fields(M*1000),np.full(8,M*1000));rate.append([p.radiative(cc),p.rate(cc,H)])
 cc=p.coefficients(p.fields(M*1000),M*1000*np.tile([1.05,.95],4));split.append([p.radiative(cc),p.rate(cc,H)])
rate=np.array(rate);split=np.array(split)
f,ax=plt.subplots(figsize=(5.5,3.1))
for j,color,label in [(0,'#124e78',r'$\mu\to e\gamma$'),(1,'#bc6b23',r'$\mu\to3e$')]:
 ax.semilogy(m,rate[:,j],color=color,label=label+', aligned');ax.semilogy(m,split[:,j],color=color,ls='--',label=label+r', $\delta=0.05$')
ax.axhline(1.5e-13,color='#124e78',ls=':',lw=.8);ax.axhline(4.60517e-15,color='#bc6b23',ls=':',lw=.8)
ax.set(xlabel=r'$M=\Lambda_D$ (TeV)',ylabel='Branching ratio',ylim=(1e-17,2e-10));ax.legend(ncol=2);f.savefig(fig/'fig_scale_physical.pdf');plt.close(f)
# Mix/pole scatter diagnostic, matching prior independent and not an uncertainty probability.
a=np.load(R/'pole_scan_10.npy');f,ax=plt.subplots(figsize=(5.5,3.1))
for data,color,label in zip(a,['#bc6b23','#124e78'],['Unmixed independent poles','Arbitrary real mixing']):ax.scatter(data[:,1],data[:,2],s=3,alpha=.25,color=color,label=label,rasterized=True)
ax.axvline(1.5e-13,color='gray',ls='--');ax.axhline(rd,color='gray',ls='--');ax.set(xscale='log',yscale='log',xlabel=r'$B(\mu\to e\gamma)$ at 5 TeV',ylabel=r'$B_{3e}/B_{e\gamma}$');ax.legend();f.savefig(fig/'fig_general_mixing.pdf');plt.close(f)
