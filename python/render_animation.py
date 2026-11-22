# Python twin of DVA_Animation.m (same equations, same layout) -> GIF + MP4 for the report
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter
from matplotlib.patches import Rectangle
m1,k1,c1,m2=200,2e6,800.0,30.0
w_op=955*2*np.pi/60; F0=50*6.3e-3/w_op*w_op**2
k2,c2=288588.818,485.1848
mag,fps,fvis=150,20,1.5; ISO=[2.3,4.5,7.1]
def cplx(w):
    X0=F0/(k1-m1*w**2+1j*w*c1)
    Z11=k1+k2-m1*w**2+1j*w*(c1+c2);Z12=-(k2+1j*w*c2);Z22=k2-m2*w**2+1j*w*c2
    D=Z11*Z22-Z12**2;return X0,F0*Z22/D,-F0*Z12/D
fH=np.linspace(8,24,800);A=np.array([cplx(2*np.pi*f) for f in fH]);A0=abs(A[:,0])*1e3;A1=abs(A[:,1])*1e3
def spring(xc,y0,y1,nc,wd):
    L=y1-y0;lead=.12*L;yz=np.linspace(y0+lead,y1-lead,2*nc+1);xz=xc+wd*np.r_[0,np.tile([1,-1],nc)];xz[-1]=xc
    return np.r_[xc,xz,xc],np.r_[y0,yz,y1]
def dash(xc,y0,y1,L0,w):
    yc=y0+.1*L0;yt=y0+.85*L0;yp=y1-.55*L0;n=np.nan
    return [xc,xc,n,xc-w,xc-w,xc+w,xc+w,n,xc,xc,n,xc-.8*w,xc+.8*w],[y0,yc,n,yt,yc,yc,yt,n,y1,yp,n,yp,yp]
fig=plt.figure(figsize=(12.8,7.2),dpi=75,facecolor='w')
axL=fig.add_axes([.03,.34,.30,.60]);axR=fig.add_axes([.35,.34,.30,.60])
G={}
for ax,key,dv in [(axL,'L',False),(axR,'R',True)]:
    ax.set_xlim(-210,210);ax.set_ylim(-30,640);ax.set_aspect('equal');ax.axis('off')
    ax.plot([-170,170],[0,0],'k',lw=2)
    for xh in range(-165,156,20): ax.plot([xh,xh+14],[0,-14],'k',lw=.8)
    P={'s1':ax.plot([],[],'k',lw=1.8)[0],'d1':ax.plot([],[],'k',lw=1.8)[0],
       'm':ax.add_patch(Rectangle((0,0),230,120,fc=(.27,.47,.72),ec='k',lw=1.5)),
       'disc':ax.plot([],[],'w',lw=2)[0],'ecc':ax.plot([],[],'o',ms=9,mfc=(.85,.1,.1),mec='w')[0],
       'lab':ax.text(0,0,'Pump 200 kg',color='w',weight='bold',ha='center',va='center',fontsize=10),
       'txt':ax.text(0,595,'',ha='center',va='center',fontsize=12,weight='bold')}
    if dv:
        P.update(s2=ax.plot([],[],'k',lw=1.8)[0],d2=ax.plot([],[],'k',lw=1.8)[0],
                 a=ax.add_patch(Rectangle((0,0),130,60,fc=(.93,.55,.15),ec='k',lw=1.5)),
                 alab=ax.text(0,0,'DVA',weight='bold',ha='center',va='center',fontsize=10))
    G[key]=P
axL.set_title('WITHOUT absorber',color=(.75,.1,.1),fontsize=15,weight='bold')
axR.set_title('WITH robust DVA (30 kg)',color=(.1,.5,.15),fontsize=15,weight='bold')
axF=fig.add_axes([.06,.07,.58,.22]);axF.axvspan(.9*w_op/2/np.pi,1.1*w_op/2/np.pi,color=(1,.88,.65),label='Speed band ±10%')
axF.plot(fH,A0,'k',lw=1.6,label='No absorber');axF.plot(fH,A1,'r',lw=2,label='Robust DVA')
axF.set_xlim(8,24);axF.set_ylim(0,.45);axF.grid(alpha=.4);axF.set_xlabel('Frequency (Hz)');axF.set_ylabel('X$_1$ (mm)');axF.legend(loc='upper right',fontsize=8)
cur=axF.axvline(w_op/2/np.pi,color='b',lw=1.5);mk0,=axF.plot([],[],'ko',ms=7);mk1,=axF.plot([],[],'ro',ms=7)
axI=fig.add_axes([.67,.07,.31,.87]);axI.axis('off')
axI.text(0,.97,'Robust Dynamic Vibration Absorber',fontsize=15,weight='bold',va='top')
for j,s in enumerate(['Pump + frame: m$_1$ = 200 kg, k$_1$ = 2.0 MN/m, ζ$_1$ = 2%',f'Unbalance: G6.3, F$_0$ = {F0:.1f} N at 955 rpm',
    f'Absorber: m$_2$ = 30 kg, k$_2$ = {k2/1e3:.1f} kN/m, c$_2$ = {c2:.0f} N·s/m','Limit: ISO 10816-3 Zone B/C = 4.5 mm/s']):
    axI.text(0,.89-.055*(j+1),s,fontsize=11,va='center')
hS=axI.text(0,.58,'',fontsize=13,weight='bold',color=(.1,.2,.6));hF=axI.text(0,.51,'',fontsize=13)
axI.text(0,.38,'Motion magnified ×150 and slowed down.\nAmplitudes and phases are the exact\nsteady-state 2-DOF results at each speed.',fontsize=10,color='.35',va='top')
axI.text(0,.12,'Python twin of DVA_Animation.m\n(same equations; MATLAB original in MATLAB Drive)',fontsize=9,color='.35',va='top')
cols=[(.1,.55,.15),(.45,.6,.1),(.9,.5,0),(.8,.1,.1)]
def readout(h,X,w):
    v=abs(X)*w/np.sqrt(2)*1e3;iz=next(i for i,l in enumerate(ISO+[1e9]) if v<=l)
    h.set_text(f'X$_1$ = {abs(X)*1e3:.3f} mm\nv$_{{rms}}$ = {v:.1f} mm/s  →  Zone {"ABCD"[iz]}');h.set_color(cols[iz])
def up(P,d1,d2,th,dv):
    yb=200+d1;x,y=spring(-55,0,yb,7,22);P['s1'].set_data(x,y);x,y=dash(55,0,yb,200,20);P['d1'].set_data(x,y)
    P['m'].set_xy((-115,yb));yc=yb+52;a=np.linspace(0,2*np.pi,60);P['disc'].set_data(30*np.cos(a),yc+30*np.sin(a))
    P['ecc'].set_data([30*np.sin(th)],[yc+30*np.cos(th)]);P['lab'].set_position((0,yb+106))
    if dv:
        ya=430+d2;yt=yb+120;x,y=spring(-30,yt,ya,5,16);P['s2'].set_data(x,y);x,y=dash(30,yt,ya,110,12);P['d2'].set_data(x,y)
        P['a'].set_xy((-65,ya));P['alab'].set_position((0,ya+30))
N1,N2=6*fps,10*fps;fs=np.linspace(8,24,N2)
def frame(n):
    th=2*np.pi*fvis*n/fps
    if n<N1: f=w_op/2/np.pi;hS.set_text('Scene 1: running at 955 rpm')
    else: f=fs[n-N1];hS.set_text('Scene 2: speed sweep 480 → 1440 rpm')
    w=2*np.pi*f;X0,X1,X2=cplx(w);e=np.exp(1j*th)
    up(G['L'],mag*(X0*e).real*1e3,0,th,False);up(G['R'],mag*(X1*e).real*1e3,mag*(X2*e).real*1e3,th,True)
    readout(G['L']['txt'],X0,w);readout(G['R']['txt'],X1,w);hF.set_text(f'Speed: {f*60:4.0f} rpm  ({f:.1f} Hz)')
    cur.set_xdata([f,f]);mk0.set_data([f],[abs(X0)*1e3]);mk1.set_data([f],[abs(X1)*1e3])
an=FuncAnimation(fig,frame,frames=N1+N2)
an.save('DVA_animation.mp4',writer=FFMpegWriter(fps=fps,bitrate=2500))
fig.set_dpi(50); an.save('DVA_animation.gif',writer=PillowWriter(fps=10),dpi=50) if False else None
