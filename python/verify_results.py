import numpy as np, json
from scipy.optimize import minimize
from scipy.integrate import solve_ivp
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
m1=200;k1=2e6;c1=2*0.02*np.sqrt(k1*m1);w_op=955*2*np.pi/60
U=50*6.3e-3/w_op;F0=U*w_op**2;wn=np.sqrt(k1/m1)
w=np.linspace(0.5*wn,1.5*wn,3001);band=(w>=0.9*w_op)&(w<=1.1*w_op);KS=np.linspace(0.85,1.15,13)
iop=np.argmin(abs(w-w_op));vr=lambda X,ww:X*ww/np.sqrt(2)*1e3
def frf(m2,k2,c2,k1=k1):
    Z11=k1+k2-m1*w**2+1j*w*(c1+c2);Z12=-(k2+1j*w*c2);Z22=k2-m2*w**2+1j*w*c2
    D=Z11*Z22-Z12**2;return np.abs(F0*Z22/D),np.abs(F0*(-Z12-Z22)/D)
def dh(m2,mu):
    f=1/(1+mu);z=np.sqrt(3*mu/(8*(1+mu)**3));wa=f*wn;return m2*wa**2,2*z*m2*wa
def wbv(m2,k2,c2):return max(vr(frf(m2,k2,c2,k1*s)[0],w)[band].max() for s in KS)
X0=np.abs(F0/(k1-m1*w**2+1j*w*c1))
out={'F0':F0,'U':U*1e6,'v0':vr(X0[iop],w_op),'X0':X0[iop]*1e3,'FT0':abs(k1+1j*w_op*c1)*X0[iop]}
mus=np.round(np.arange(0.02,0.301,0.01),2);vb=[]
for m in mus:
    k2,c2=dh(m*m1,m);vb.append(vr(frf(m*m1,k2,c2)[0],w)[band].max())
vb=np.array(vb);out['muMin']=float(mus[np.argmax(vb<=4.5)]);out['vb']={str(m):float(v) for m,v in zip(mus,vb)}
mu=0.15;m2=30
D=[['Undamped',m2*wn**2,0.0],['Den Hartog',*dh(m2,mu)]]
r=minimize(lambda p:frf(m2,*np.exp(p))[0].max(),np.log(D[1][1:]),method='Nelder-Mead',options={'xatol':1e-8,'fatol':1e-12,'maxiter':10000,'maxfev':10000});D.append(['Minimax',*np.exp(r.x)])
r=minimize(lambda p:wbv(m2,*np.exp(p)),np.log(D[1][1:]),method='Nelder-Mead',options={'xatol':1e-8,'fatol':1e-12,'maxiter':10000,'maxfev':10000});D.append(['Robust',*np.exp(r.x)])
fig,ax=plt.subplots(figsize=(7.5,4.2));ax.axvspan(0.9*w_op/2/np.pi,1.1*w_op/2/np.pi,color='orange',alpha=.15,label='Speed band ±10%')
ax.plot(w/2/np.pi,X0*1e3,'k',lw=1.6,label='No absorber')
sty=['b--','g-.','m:','r-'];out['designs']={}
for (n,k2,c2),s in zip(D,sty):
    X,Xr=frf(m2,k2,c2);ax.plot(w/2/np.pi,X*1e3,s,lw=1.6,label=n)
    out['designs'][n]=dict(k2=k2,c2=c2,zeta=c2/(2*np.sqrt(k2*m2)),f=np.sqrt(k2/m2)/wn,peak=X.max()*1e3,vop=vr(X[iop],w_op),vband=vr(X,w)[band].max(),vwc=wbv(m2,k2,c2),rel=Xr.max()*1e3,FT=abs(k1+1j*w_op*c1)*X[iop])
ax.set_ylim(0,0.5);ax.set_xlabel('Frequency (Hz)');ax.set_ylabel('Machine amplitude X1 (mm)');ax.grid(alpha=.3);ax.legend(fontsize=8);ax.set_title('Machine response with and without absorber')
fig.tight_layout();fig.savefig('fig1_frf.png',dpi=150)
R=out['designs']['Robust'];k2r,c2r=R['k2'],R['c2']
out['tol_k2']={s:wbv(m2,k2r*s,c2r) for s in (0.9,0.95,1.05,1.1)};out['tol_c2']={s:wbv(m2,k2r,c2r*s) for s in (0.5,0.75,1.5,2)}
def run(k2,c2,T=8,te=12):
    def rhs(t,y):
        wt,th=(w_op*t/T,w_op*t*t/(2*T)) if t<T else (w_op,w_op*T/2+w_op*(t-T))
        F=U*wt**2*np.sin(th);x1,v1,x2,v2=y
        return [v1,(F-c1*v1-k1*x1-c2*(v1-v2)-k2*(x1-x2))/m1,v2,(-c2*(v2-v1)-k2*(x2-x1))/m2 if k2 else 0]
    return solve_ivp(rhs,[0,te],[0]*4,max_step=1e-3,rtol=1e-7,atol=1e-11)
fig,ax=plt.subplots(figsize=(7.5,3.6));out['runup']={}
for n,k2,c2,c in [('No absorber',0,0,'k'),('Undamped',m2*wn**2,0,'b'),('Robust',k2r,c2r,'r')]:
    s=run(k2,c2);ax.plot(s.t,s.y[0]*1e3,c,lw=.5,label=n)
    out['runup'][n]=dict(mx=np.abs(s.y[0]).max()*1e3,ss=np.abs(s.y[0][s.t>10]).max()*1e3,rel=np.abs(s.y[0]-s.y[2]).max()*1e3 if k2 else 0)
ax.set_xlabel('Time (s)');ax.set_ylabel('x1 (mm)');ax.grid(alpha=.3);ax.legend(fontsize=8);ax.set_title('Start-up 0 → 955 rpm in 8 s')
fig.tight_layout();fig.savefig('fig2_startup.png',dpi=150)
fig,ax=plt.subplots(figsize=(7.5,3.6));ax.plot(mus,vb,'o-',ms=3);ax.axhline(4.5,color='r',ls='--',label='Zone B/C 4.5 mm/s');ax.axhline(2.3,color='g',ls=':',label='Zone A/B 2.3 mm/s')
ax.set_xlabel('Mass ratio μ = m2/m1');ax.set_ylabel('Worst v_rms in band (mm/s)');ax.grid(alpha=.3);ax.legend(fontsize=8);ax.set_title('Mass-ratio trade-off (Den Hartog tuning)')
fig.tight_layout();fig.savefig('fig3_mass_ratio.png',dpi=150)
G=79.3e9;rho=7850;d=8e-3;Dm=48e-3;C=Dm/d;ke=k2r/4;na=G*d**4/(8*Dm**3*ke);nt=na+2;Ls=nt*d
Fst=m2*9.81/4;yst=Fst/ke;Fa=ke*R['rel']*1e-3;Ks=1+0.5/C;Kw=(4*C-1)/(4*C-4)+0.615/C
tm=Ks*8*Fst*Dm/(np.pi*d**3)/1e6;ta=Kw*8*Fa*Dm/(np.pi*d**3)/1e6;Sut=1610/(8**0.193);Ssu=0.67*Sut;Sse=241/(1-379/Ssu)
L0=Ls+yst+R['rel']*1e-3+0.15*yst+0.01;ms=rho*np.pi*d**2/4*np.pi*Dm*na
out['spring']=dict(ke=ke,na=na,nt=nt,Ls=Ls*1e3,yst=yst*1e3,L0=L0*1e3,L0D=L0/Dm,tm=tm,ta=ta,Sut=Sut,nF=1/(ta/Sse+tm/Ssu),fsurge=0.5*np.sqrt(ke/ms),block_t=m2/rho/0.04*1e3)
out['damper']=dict(eta=2*R['zeta'],fa=R['f']*wn/2/np.pi,Fd=c2r*R['f']*wn*R['rel']*1e-3)
json.dump(out,open('final.json','w'),indent=1,default=float)
d=dict(out);d.pop('vb');print(json.dumps(d,indent=1,default=float))
