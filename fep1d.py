#!/usr/bin/env python
"""
fep1d.py - a script for the analysis of reaction coordinates.

*  Copyright (C) 2014  Sergei Krivov, Polina Banushkina <krivov@yahoo.com>
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""
import math
import Gnuplot

#----the-gnuplot-py-1.7-mousesupport.patch-----------------------------------
import Gnuplot.gp as gp
def test_mouse():
    import os,tempfile,commands
    
    tmpname = tempfile.mktemp()
    tfile = open(tmpname,"w")
    tfile.write("set mouse")
    tfile.close()
    msg = commands.getoutput(gp.GnuplotOpts.gnuplot_command + " " +
                             tmpname)
    os.unlink(tmpname)
    if msg:  # Gnuplot won"t print anything if it has mouse support
        has_mouse = 0
    else:
        has_mouse = 1
    return has_mouse

def new_init(self, filename=None, persist=None, debug=0, mouse=None):
        if mouse is None:
            mouse = test_mouse()
        if mouse:
            gp.GnuplotOpts.prefer_inline_data = 0
            gp.GnuplotOpts.prefer_fifo_data = 0

        if filename is None:
            self.gnuplot = gp.GnuplotProcess(persist=persist)
        else:
            if persist is not None:
                raise Errors.OptionError(
                    'Gnuplot with output to file does not allow '
                    'persist option.')
            self.gnuplot = _GnuplotFile(filename)

        self._clear_queue()
        self.debug = debug
        self.plotcmd = 'plot'
        if mouse:
            self("set mouse")
        self("set terminal %s" % gp.GnuplotOpts.default_term)

Gnuplot.Gnuplot.__init__=new_init
#----------------------------------------------------------------------


def readRC(name,ind,tmin=None,tmax=None):
  f=open(name)
  lx=[]
  t=0
  for s in f:
    t=t+1
    if tmin!=None and tmin>t:continue
    if tmax!=None and t>tmax:break
    s=s.split()
    lx.append(float(s[ind]))
  return lx  

def writeRC(lx,name):
  f=open(name,'w')
  for x in lx:f.write('%g\n' %x)
  f.close()
  
def writeXY(lxy,name):
  f=open(name,'w')
  for x,y in lxy:f.write('%g %g\n' %(x,y))
  f.close()

def compZh(lx,dx,dt=1):
  """ computes Zh[x], returns dictionary Zh[x]"""
  zh={}
  for x in lx:    
    x=math.floor(x/dx)*dx
    zh[x]=zh.get(x,0)+1
  for x in zh:zh[x]=float(zh[x])/dt/dx
  return zh

def compZhadapt(lx,dx,dt=1,zhmin=100):
  """ computes Zh[x], returns dictionary Zh[x]"""
  zh={}
  for x in lx:    
    x=math.floor(x/dx)*dx
    zh[x]=zh.get(x,0)+1
  for x in zh:zh[x]=float(zh[x])/dt/dx
  keys=zh.keys()
  keys.sort()
  lzh=[]
  sx=0
  szh=0
  sn=0
  zhmin=float(zhmin)/dx
  for x in keys:
    szh=szh+zh[x]
    sx=sx+x*zh[x]
    sn=sn+1
    if szh>zhmin:
      lzh.append((sx/szh,szh/sn))
      sx,szh,sn=0,0,0
  if sn>0:lzh.append((sx/szh,szh/sn))
  return lzh

def compZhadapt2(lx,dx,dt=1,zhmin=100):
  """ computes Zh[x], returns dictionary Zh[x]"""
  zh={}
  for x in lx:    
    x=math.floor(x/dx)*dx
    zh[x]=zh.get(x,0)+1
  for x in zh:zh[x]=float(zh[x])/dt/dx
  keys=zh.keys()
  keys.sort()
  lzh=[]
  sx=0
  szh=0
  sn=0
  zhmin=float(zhmin)/dx
  i0=0
  for i in range(len(keys)):
    x=keys[i]
    szh=szh+zh[x]
    sx=sx+x*zh[x]
    sn=sn+1
    x0=keys[i0]
    while szh-zh[x0]>zhmin: # subtract first point
      szh=szh-zh[x0]
      sx=sx-x0*zh[x0]
      sn=sn-1
      i0=i0+1
      x0=keys[i0]
    if szh>zhmin:
      lzh.append((sx/szh,szh/sn))
  lzh.append((sx/szh,szh/sn))
  return lzh

def compZhkde(lx,dx,dt=1,bw=0.01):
  """ computes Zh[x], returns dictionary Zh[x] using kernel density estimator"""
  from scipy import stats
  import numpy 
  z=len(lx)
  gkde=stats.gaussian_kde(lx,bw)
  xmin=min(lx)
  xmax=max(lx)
  np=int((xmax-xmin)/dx)
  ind=numpy.linspace(xmin,xmax,np)
  pdf=gkde.evaluate(ind)
  lzh=[(ind[i],pdf[i]*z) for i in range(np)]
  return lzh


def compZc(lx,dx=None,dt=1,eps=None):
  """computes Zc(x), returns dictionary Zc[x]"""
  dzc={}
  tmax=len(lx)
  for i in range(0,tmax-dt):
    x=lx[i+dt]
    lastx=lx[i]
    if dx!=None: 
      x=math.floor(x/dx)*dx
      lastx=math.floor(lastx/dx)*dx
    if lastx<x:
      dzc[lastx]=dzc.get(lastx,0)+1
      dzc[x]=dzc.get(x,0)-1
    else:  
      dzc[x]=dzc.get(x,0)+1
      dzc[lastx]=dzc.get(lastx,0)-1
  keys=dzc.keys()
  keys.sort()
  zc={}
  z=0
  for x in keys:
    zc[x]=float(z)/dt/2
    z=z+dzc[x]
    if eps!=None: 
      zc[x+eps]=float(z)/dt/2  
      if zc[x+eps]<1e-8:zc[x+eps]=0
  return zc        

def compZcr(lx,r,dx=None,dt=1,eps=None):
  """computes Zc,r(x), returns dictionary Zc,r[x]"""
  dzc={}
  tmax=len(lx)
  for i in range(0,tmax-dt):
    x=lx[i+dt]
    lastx=lx[i]
    d=abs(x-lastx)
    if r<0 and d<dx: continue 
    else: d=d**r
    if dx!=None: 
      x=math.floor(x/dx)*dx
      lastx=math.floor(lastx/dx)*dx
    if lastx<x:
      dzc[lastx]=dzc.get(lastx,0)+d
      dzc[x]=dzc.get(x,0)-d
    else:  
      dzc[x]=dzc.get(x,0)+d
      dzc[lastx]=dzc.get(lastx,0)-d
  keys=dzc.keys()
  keys.sort()
  zc={}
  z=0
  for x in keys:
    zc[x]=float(z)/dt/2
    z=z+dzc[x]
    if eps!=None:
      zc[x+eps]=float(z)/dt/2  
      if zc[x+eps]<1e-8:zc[x+eps]=0
  return zc        

def compZcrMSM(ekn,r,dx=None,eps=None):
  """computes Zc,r(x) from an MSM, returns dictionary Zc,r[x]"""
  dzc={}
  for x,y in ekn:
    d=abs(y-x)**r*ekn[(x,y)]
    if dx!=None: 
      x=math.floor(x/dx)*dx
      y=math.floor(y/dx)*dx
    if y<x:
      dzc[y]=dzc.get(y,0)+d
      dzc[x]=dzc.get(x,0)-d
    else:  
      dzc[x]=dzc.get(x,0)+d
      dzc[y]=dzc.get(y,0)-d
  keys=dzc.keys()
  keys.sort()
  zc={}
  z=0
  for x in keys:
    zc[x]=float(z)/2
    z=z+dzc[x]
    if eps!=None:
      zc[x+eps]=float(z)/2
      if zc[x+eps]<1e-8:zc[x+eps]=0
  return zc        


def tonatural(lx,dx,mdydx=100):
#  Zc=Zh*(D*dt/pi)**0.5; dt=1
#  D**0.5=pi**0.5*Zc/Zh
  zh=compZh(lx,dx)
  zc=compZc(lx,dx)
  keys=zh.keys()
  keys.sort()
  spi=1./math.sqrt(math.pi) # pi**(-0.5)
  y=0
  x2y={}
  dydx={}
  for x in keys:
    dydx[x]=0
    if zc[x]>0:
      dydx[x]=float(zh[x])/zc[x]*spi
      if zc[x]<20:dydx[x]=min(dydx[x],mdydx) # to avoid gaps at small statistics
    y=y+dydx[x]*dx
    x2y[x]=y
  a=dydx.values()
  a.sort()
  ly=[]
  for x in lx:
    x0=math.floor(x/dx)*dx
    ly.append(x2y[x0]+(x-x0)*dydx[x0])
  return ly  

def tonatural2(lx,dx,zhmin=1000, returnx2y=0):
#  Zc1=dt*D*Zh; dt=1
#  D=Zc1/Zh
#  dy/dx=D^-1/2
  zh=compZh(lx,dx)
  zc1=compZcr(lx,1,dx)
  keys=zh.keys()
  keys.sort()
  y=0
  x2y={}
  dydx={}
  s=0
  s2=0
  s3=0
  lxt=[]
  for i in range(len(keys)):
    x=keys[i]
    lxt.append(x)
    s+=zh[x]
    s2+=zh[x]*zc1[x]
    s3+=1.
    if s*dx>zhmin:
      s2=float(s2)/s
      if s2<1:cdydx=0
      else: cdydx=(float(s)/s3/s2)**0.5
      for x in lxt:
        dydx[x]=cdydx
        y=y+dydx[x]*dx
        x2y[x]=y
      s,s2,s3,lxt=0,0,0,[]
  for x in lxt: # leftovers
    dydx[x]=cdydx
    y=y+dydx[x]*dx
    x2y[x]=y
  a=dydx.values()
  a.sort()
  ly=[]
  for x in lx:
    x0=math.floor(x/dx)*dx
    ly.append(x2y[x0]+(x-x0)*dydx[x0])
  if returnx2y==1:return ly,x2y 
  return ly

def tonatural3(lx,dx, returnx2y=0):
#  Zc1=dt*D*Zh; Zh=2Zc-1; dt=1
#  D=Zc1/(2Zc-1)
  zcm1=compZcr(lx,-1,dx)
  zc1=compZcr(lx,1,dx)
  keys=zcm1.keys()
  keys.sort()
  y=0
  x2y={}
  dydx={}
  xl=0
  for x in keys:
    dydx[x]=0
    if zc1[x]>0:
      dydx[x]=(float(2*zcm1[x])/zc1[x])**0.5
    y=y+dydx[x]*(x-xl)
    x2y[x]=y
    xl=x
  a=dydx.values()
  a.sort()
  ly=[]
  for x in lx:
    x0=math.floor(x/dx)*dx
    ly.append(x2y[x0]+(x-x0)*dydx[x0])
  if returnx2y==1:return ly,x2y 
  return ly 
  
def toZa(lx,dx):
  zh=compZh(lx,dx)
  keys=zh.keys()
  keys.sort()
  zt=sum(zh.values()) # total Z
  x2y={}
  y=0
  for x in keys:
    y=y+float(zh[x])/zt
    x2y[x]=y
  ly=[x2y[math.floor(x/dx)*dx] for x in lx]
  return ly  

def comppfold(lx,dx,x0,x1):
  zh=compZh(lx,dx)
  zc=compZc(lx,dx)
  x2y={}
  keys=zh.keys()
  keys.sort()
  if x0>x1:x0,x1,rev=x1,x0,1
  else: rev=0
  lpfold=[]
  y=0
  for x in keys:
    if x>x0 and x<x1:
      if zc[x]>0:y=y+float(zh[x])/zc[x]/zc[x]
    lpfold.append((x,y))
  s=lpfold[-1][1]-lpfold[0][1]
  if rev==0:lpfold=[(x,y/s) for x,y in lpfold]  
  else:lpfold=[(x,1-y/s) for x,y in lpfold]  
  return lpfold  

def topfold(lx,dx,x0,x1):
  zh=compZh(lx,dx)
  zc=compZc(lx,dx)
  x2y={}
  keys=zh.keys()
  keys.sort()
  if x0>x1:x0,x1,rev=x1,x0,1
  else: rev=0
  y=0
  for x in keys:
    if x>x0 and x<x1:
      if zc[x]>0:y=y+float(zh[x])/zc[x]/zc[x]
    x2y[x]=y
  s=x2y[math.floor(x1/dx)*dx]-x2y[math.floor(x0/dx)*dx]
  ly=[]
  for x in lx:
    if x<x0:x=x0
    if x>x1:x=x1
    x=math.floor(x/dx)*dx
    if rev==0: ly.append(x2y[x]/s)
    else:ly.append(1-x2y[x]/s)
  return ly  

  
def comppfoldMSM(lx,dx,x0,x1):
  """pfold coordinate from an MSM"""
  import scipy
  import scipy.sparse
  from scipy.sparse.linalg import spsolve
  
  if x0>x1:x0,x1,rev=x1,x0,1
  else: rev=0
  ekn={}
  iso={}
  ni={}
  lasti=None
  for x in lx:
    if x<x0:x=x0
    if x>x1:x=x1
    x=math.floor(x/dx)*dx
    if iso.has_key(x): i=iso[x]
    else:
      i=len(iso)
      iso[x]=i
    ni[i]=ni.get(i,0)+1
    if lasti!=None: ekn[(lasti,i)]=ekn.get((lasti,i),0)+1
    lasti=i
  size=len(iso)
  a=scipy.sparse.lil_matrix((size,size),dtype='d')
  b=scipy.zeros(size,dtype='d')
  i1p=iso[math.floor(x1/dx)*dx]
  i0p=iso[math.floor(x0/dx)*dx]
  for i,j in ekn:
    if i==i1p or i==i0p:continue
    a[i,j]=float(ekn[(i,j)])/ni[i]
  for i in range(size):a[i,i]-=1
  a[i0p,i0p]=1
  a[i1p,i1p]=1
  b[i1p]=1
  a=a.tocsr()
  pf=spsolve(a,b)
  pfx={}
  if rev==0:
    for x in iso:pfx[x]=pf[iso[x]]
  else:  
    for x in iso:pfx[x]=1-pf[iso[x]]
  return pfx

def topfoldMSM(lx,dx,x0,x1):
  """transform to pfold coordinate from an MSM"""
  pf=comppfoldMSM(lx,dx,x0,x1)
  lpf=[]
  if x0>x1:x0,x1=x1,x0
  for x in lx:
    if x<x0:x=x0
    if x>x1:x=x1
    x=math.floor(x/dx)*dx
    lpf.append(pf[x])
  return lpf  

def compD(lx,dx,dt=1):
  """ compute D(x), returns list of pairs [x,D(x)]"""
  zh=compZh(lx,dx,dt)
  zc=compZc(lx,dx,dt)
  ld=[]
  keys=zh.keys()
  a=math.pi/dt
  for x in keys:
    d=(float(zc[x])/zh[x])**2*a
    ld.append((x,d))
  return ld

def compD2(lx,dx,dt=1):
  """ compute D(x), returns list of pairs [x,D(x)]"""
  zc1=compZcr(lx,1,dx,dt)
  zc2=compZcr(lx,-1,dx,dt)
  ld=[]
  keys=zc1.keys()
  for x in keys:
    if zc2[x]>0:
      d=float(zc1[x])/zc2[x]/2/dt
      ld.append((x,d))
  return ld
 
  
def compalpha(lx,dx,dt1,dt2):  
  zc1=compZc(lx,dx,dt1)
  zc2=compZc(lx,dx,dt2)
  keys=zc1.keys()
  lalpha=[]
  ldt=math.log(dt1)-math.log(dt2)
  for x in keys:
    if zc1[x]>0 and zc2[x]>0:
      al=1+(math.log(zc1[x])-math.log(zc2[x]))/ldt
      lalpha.append((x,al))
  return lalpha

def compptpx(lx,dx,x0,x1):
  """compute p(TP|x)"""
  zheq=compZh(lx,dx)
  zhtp={}
  def addZh(zh,ltpx):
    for x in ltpx:    
      x=math.floor(x/dx)*dx
      zh[x]=zh.get(x,0)+1
  b=3
  ltpx=[]
  if x0>x1:x0,x1=x1,x0
  for x in lx:
    ltpx.append(x)
    if b==3:
      if x<=x0:b=0
      if x>=x1:b=1
    if x<=x0:
      if b==1: addZh(zhtp,ltpx)
      b=0
      ltpx=[]  
    if x>=x1:
      if b==0: addZh(zhtp,ltpx)
      b=1
      ltpx=[]  
  ptpx=[]
  for x in zheq:
    ptpx.append((x,float(zhtp.get(x,0))/dx/zheq[x]))
  ptpx.sort()  
  return ptpx  

def comp_ekn_tp(traj,x0,x1,dt=1,dx=None,react=0):
  """computes MSM by using transition paths""" 
  def process(traj):
    if traj[0]==traj[-1] and react==1:return
    n=len(traj)
    if n<2:return
    for i in range(1,n): # from i-dt to i
      j=i-dt
      if j<0:j=0
      key=traj[j],traj[i]
      ekn[key]=ekn.get(key,0)+1
    for i in range(max(n-dt,1),n-1):
      key=traj[i],traj[-1]
      ekn[key]=ekn.get(key,0)+1
    if dt>n-1:
      key=traj[0],traj[-1]
      ekn[key]=ekn.get(key,0)+dt-n+1
    
  ekn={}
  ok=False
  lx=[]
  if dx!=None: 
    x0=math.floor(x0/dx)*dx
    x1=math.floor(x1/dx)*dx
  for x in traj:
    if dx!=None:x=math.floor(x/dx)*dx
    if x<=x0:
      lx.append(x0)
      if ok:process(lx)
      lx=[x0]
      ok=True
      continue
    if x>=x1:  
      lx.append(x1)
      if ok:process(lx)
      lx=[x1]
      ok=True
      continue
    lx.append(x)
  process(lx)
  for ij in ekn:ekn[ij]=float(ekn[ij])/dt   
  return ekn


def comp_mfpte(lx,x0,x1,dt=1):
  if x0>x1:x0,x1=x1,x0
  b=3
  nij1={}
  nij2={}
  t=0
  for i in xrange(0,len(lx),dt):
    x=lx[i]
    t+=dt
    if b==3:
      if x<=x0:b=0
      if x>=x1:b=1
    if x<=x0:
      key=b,0
      nij1[key]=nij1.get(key,0)+t
      nij2[key]=nij2.get(key,0)+1
      b=0
      t=0
    if x>=x1:
      key=b,1
      nij1[key]=nij1.get(key,0)+t
      nij2[key]=nij2.get(key,0)+1
      b=1
      t=0
  t01e=float(nij1[(0,0)]+nij1[(0,1)])/nij2[(0,1)]   
  t10e=float(nij1[(1,0)]+nij1[(1,1)])/nij2[(1,0)]
  tp01e=float(nij1[(0,1)])/nij2[(0,1)]   
  tp10e=float(nij1[(1,0)])/nij2[(1,0)]
  n10e=nij2[(1,0)]
  n01e=nij2[(0,1)]
  print 'between %g and %g estimated from trajectory with dt %g' %(x0,x1,dt) 
  print 'N < %g   > %g' %(n10e,n01e)  
  print 'mfpt < %g   > %g' %(t10e,t01e)  
  print 'mtpt < %g   > %g' %(tp10e,tp01e)
  return n10e,n01e,t10e,t01e,tp10e,tp01e  
  

def comp_mfpt(lx,dx,ldt):
  """ compute mfpt using the Kramers equation"""
  zh=compZh(lx,dx)
  zc1=compZcr(lx,1,dx,1)
  keys=zh.keys()
  while 1:
    s=raw_input('enter x0,x1\n')
    s=s.split()
    try: x0,x1=float(s[0]),float(s[1])
    except : break 
    if x0>x1:x0,x1=x1,x0
    t=0
    za=0
    keys.sort()
    s=0
    xl=None
    for x in keys:
      if x>x1:break
      if x<x0:continue
      if xl!=None: s=s+(1/zc1[x]+1/zc1[xl])/2*(x-xl)
      xl=x
    NAB=1./s
    s2=0
    mfpt01=0
    mfpt10=0
    mtpt=0
    xl=None
    for x in keys:
      if x0<x and x<x1 and xl!=None: 
        s2=s2+(1/zc1[x]+1/zc1[xl])/2*(x-xl)
      xl=x
      q=s2/s
      if q<0:q=0
      if q>1:q=1
      mfpt01=mfpt01+zh[x]*(1-q)*dx
      mfpt10=mfpt10+zh[x]*q*dx
      mtpt=mtpt+zh[x]*q*(1-q)*dx
    mfpt01=mfpt01/NAB
    mfpt10=mfpt10/NAB
    mtpt=mtpt/NAB
      
    print 'between %g and %g estimated analytically with Z_C1 and q' %(x0,x1) 
    print 'N %g ' %(NAB)  
    print 'mfpt < %g   > %g ' %(mfpt10,mfpt01)  
    print 'mtpt %g  ' %(mtpt)  
    print 
    for dt in ldt:comp_mfpte(lx,x0,x1,dt)
    print


def comp_mfpt_q(lq,dx,ldt=[1,],xq=3):
  """ compute properties using the diffusive model along q and directly from trajectories"""
  
  print xq
  if xq==2 or xq==3:
    if xq==2:lx,q2x=natural2(lq,dx,returnx2y=1)
    if xq==3:lx,q2x=natural3(lq,dx,returnx2y=1)
    zc=compZcr(lx,-1,0.1,1,1e-5)
    lxy=[(x,-math.log(2*float(zc[x]))) for x in zc if zc[x]>0]
    lxy.sort()
    g=Gnuplot.Gnuplot()
    g.plot(Gnuplot.Data(lxy, with_='lines lw 1'))
  while 1:
    if xq==2 or xq==3: s=raw_input('enter x0,x1\n')
    else: s=raw_input('enter q0,q1\n')
    s=s.split()
    try: x0,x1=float(s[0]),float(s[1])
    except : break 
    if x0>x1:x0,x1=x1,x0
    
    if xq==2 or xq==3:  
      q0=min([(abs(q2x[q]-x0),q) for q in q2x])[1]
      q1=min([(abs(q2x[q]-x1),q) for q in q2x])[1]
    else:
      q0=x0
      q1=x1
    
    dq2=0
    for i in xrange(len(lq)-1):dq2+=(lq[i+1]-lq[i])**2
    N10=dq2/2/(q1-q0)
    print 'N %g<->%g estimated as NAB/(q1-q0) %g' %(q0,q1,N10)  
    def comp_mtpt(lq,q0,q1,N10):
      s1=0
      s2=0
      s3=0
      for q in lq:
        q=(q-q0)/(q1-q0)
        if q<0:q=0
        if q>1:q=1
        s1=s1+q*(1-q)
        s2=s2+q
        s3=s3+1-q
      return s1/N10,s2/N10,s3/N10
    mtpt,mfpt1,mfpt2=comp_mtpt(lq,q0,q1,N10)
    print 'mfpt %g<->%g estimated as <q>/J and <1-q>/J %g %g' %(q0,q1,mfpt1,mfpt2)
    print 'mtpt %g<->%g estimated as <q(1-q)>/J %g' %(q0,q1,mtpt)
    print 
    if xq==2 or xq==3:
      for dt in ldt: 
        n10e,n01e,t10e,t01e,tp10e,tp01e=comp_mfpte(lx,x0,x1,dt)
        print N10/(n10e+n01e)*2-1,mfpt1/t10e-1,mfpt2/t01e-1, mtpt/(tp10e+tp01e)*2-1
    else:  
      for dt in ldt: 
        n10e,n01e,t10e,t01e,tp10e,tp01e=comp_mfpte(lq,q0,q1,dt)
        print N10/(n10e+n01e)*2-1,mfpt1/t10e-1,mfpt2/t01e-1, mtpt/(tp10e+tp01e)*2-1
    print


def fep1d(ldat, cfep=1, cfep1=0, cfepr=0, hfep=0, D=0, alpha=None, pfold=0, pfoldMSM=0, ptpx=0, ldt=[1], x0=None, x1=None, dx=0.1, transformto=None, writeps=0, writexy=0, writerc=0, testoptimality=0, tmin=None, tmax=None, mfpt=0, epszc=None, glines=None, nameout=None, mfptq=None):
  """ 
  fep1d.py - a script for the analysis of reaction coordinates
  Authors: Polina Banushkina, Sergei Krivov
  License: GPL
  
  *  fep1d.py  rc.dat  -- mandatory argument
  *  optional arguments are given in the form 
  fep1d.py  rc.dat  --argument1=value1 --argument2=value2 
  
  Function fep1d() has the following values by default
  fep1d(ldat,cfep=1,cfep1=0,hfep=0,D=0,alpha=None,pfold=0,pfoldMSM=0,ptpx=0,
  ldt=[1],x0=None,x1=None,dx=0.1,transformto=None,writeps=0,writexy=0,writerc=0,
  testoptimality=0,tmin=None,tmax=None,mfpt=0,glines=None,nameout=None):

  --cfep=1  - computes cut based (CFEP) profile
  --cfep1=1 - computes cut based (CFEP_1) profile
  5 different methods have been impleentd to estimate the density of states  
  --hfep=1  - computes histogram based (HFEP) conventional profile. 
  --hfep=2  - computes histogram based (HFEP) conventional profile with adpative choice of dx.
  --hfep=3  - estimates the denisty of states as 2Z_{C,-1}
  --hfep=4  - computes histogram based (HFEP) conventional profile with another adpative choice of dx.
  --hfep=5  - estimates the density of states using kernel density estimator gkde from scipy.stats module
  --hfep=[1,3] - e.g., can be used to plot a few estimates simultaneosly. 
  --dx='size' of histogram bins
  --ldt=[dt] -- the time step.
  
  
  Reaction coordinate transformations:
  --transformto=natural, natural2, natural3
  --transformto=Za 
  --transformto=pfold 
  --transformto=pfoldMSM
    For pfold and pfoldMSM transformation arguments  --x0=value-of-x0 
    --x1=value-of-x1 and --dx=     should be specified
 
  Writing to files:
  --nameout  -- the prefix for output files, default value - the name
      of the first reaction coordinate file
  --writeps=1 - prints the plot to postscript file 'nameout.fep1d.eps'
  --writexy=1 - writes the x,y coordinates of the plots to file 'nameout.*.xy'
  --writerc=1 -- writes transformed RC to file 'nameout.transformto.dat'


  Analysis:
  --D=1 -- computes diffusion coefficient 
  --alpha=[dt1,dt2] - list of time steps dt1 and dt2 to build the profiles. 
  --ptpx=1 -- computes the probability of being on the transition path. 
    Arguments --x0= , --x1=  and  --dx= should be specified
  --testoptimality=1 -- tests the optimality of the reaction coordinate. 
    Arguments --x0= and --x1= should be specified (x0<x1)
  --testoptimality=2 -- tests the optimality of the RC, without transforming it first to q(RC). Usefull to confirm that the RC is the committor.

  
  Computing pfold:
  --pfold=1 
  --pfoldMSM=1 
    Arguments --x0= , --x1=  and  --dx= should be specified.
 
  --mfpt=1 - computes the mean first passage time, 
    mean transitions path time and  number of transitions from the free energy profile and compares with those computed directly from the trajectory.

  --mfptq=1 - the same as mfpt, but the RC is assumed to be the committor.
  
  --glines="['set xrange [5:13]','set yrange [60:180]']" - list of commands
    to gnuplot
  """ 
  if len(ldat)==0:
    print fep1d.__doc__
    return
  g=Gnuplot.Gnuplot()
  lg=[]
  if nameout==None:nameout=ldat[0]
  if epszc!=None and epszc>dx:
    print 'epszc>dx, redefined to dx/2'
    epszc=dx/2.
  for dat in ldat:
    if ':' in dat:
      dat,ind=dat.split(':')
      ind=int(ind)-1
    else:
      ind=0  
    lx=readRC(dat,ind,tmin,tmax) 
    nameout=dat
#-------- optimality test of rc -----------------------------
    if testoptimality!=0:
      xlabel='x'
      if testoptimality==1: xlabel='pfoldMSM'
      ylabel='F_{C,1}/kT'  
      g('set xlabel \'%s\' ' %xlabel)
      g('set ylabel \'%s\' ' %ylabel)
      if len(ldt)==1:ldt=[2**i for i in range(17)]
      if testoptimality==1: lx=topfoldMSM(lx,dx,x0,x1)  
      for dt in ldt:
        ekn=comp_ekn_tp(lx,0.0,1.0,dt,dx=0.0001)
        zc=compZcrMSM(ekn,1,eps=epszc)
        lxy=[(x,-math.log(float(zc[x]))) for x in zc if zc[x]>0]
        lxy.sort()
        if writexy==1:writeXY(lxy,'%s.F_C,1.dt%s.xy' %(nameout,dt))
        lg.append(Gnuplot.Data(lxy, with_='lines lw 5', title='%s F_{C,1} dt=%g' %(nameout,dt)))
        g.plot(*lg)
    else:
      ylabel='F/kT'  
#-------- rc transformations ---------------------------------
      xlabel='x'
      if transformto=='natural':
        lx=tonatural(lx,dx)
        xlabel='natural'
      elif transformto=='natural2':
        lx=tonatural2(lx,dx)
        xlabel='natural2'
      elif transformto=='natural3':
        lx=tonatural3(lx,dx)
        xlabel='natural3'
      elif transformto=='Za':
        lx=toZa(lx,dx)
        xlabel='Za'
      elif transformto=='pfold':
        lx=topfold(lx,dx,x0,x1)
        xlabel='pfold'
      elif transformto=='pfoldMSM':
        lx=topfoldMSM(lx,dx,x0,x1)
        xlabel='pfold'
      if writerc==1:writeRC(lx,'%s.%s.dat' %(nameout,transformto))  
    
#------------- profiles ----------------------------------------  
      for dt in ldt:
        pref='dt=%i' %dt
        prefw='.dt=%i' %dt
        if len(ldt)==1: pref,prefw='',''
        if not isinstance(hfep,list):hfep=[hfep,]
        if hfep!=0:
          if 1 in hfep:
            zh=compZh(lx,dx,dt)
            lxy=[(x,-math.log(float(zh[x]))) for x in zh if zh[x]>0]
            lxy.sort()
            if writexy==1:writeXY(lxy,'%s.F_H%s.xy' %(nameout,prefw))
            lg.append(Gnuplot.Data(lxy, with_='lines lw 5', title='%s F{H1} %s' %(nameout,pref)))
          if 2 in hfep:
            lzh=compZhadapt(lx,dx,dt,zhmin=1000)
            lxy=[(x,-math.log(float(zh))) for x,zh in lzh if zh>0]
            lxy.sort()
            if writexy==1:writeXY(lxy,'%s.F_H%s.xy' %(nameout,prefw))
            lg.append(Gnuplot.Data(lxy, with_='lines lw 5', title='%s F{H2} %s' %(nameout,pref)))
          if 3 in hfep:
            zc=compZcr(lx,-1,dx,dt,eps=epszc)
            lxy=[(x,-math.log(2*float(zc[x]))) for x in zc if zc[x]>0]
            lxy.sort()
            if writexy==1:writeXY(lxy,'%s.F_H%s.xy' %(nameout,prefw))
            lg.append(Gnuplot.Data(lxy, with_='lines lw 5', title='%s F{H3} %s' %(nameout,pref)))
          if 4 in hfep:
            lzh=compZhadapt2(lx,dx,dt,zhmin=1000)
            lxy=[(x,-math.log(float(zh))) for x,zh in lzh if zh>0]
            lxy.sort()
            if writexy==1:writeXY(lxy,'%s.F_H%s.xy' %(nameout,prefw))
            lg.append(Gnuplot.Data(lxy, with_='lines lw 5', title='%s F{H4} %s' %(nameout,pref)))
          if 5 in hfep:
            lzh=compZhkde(lx,dx,dt,bw=0.003)
            lxy=[(x,-math.log(float(zh))) for x,zh in lzh if zh>0]
            lxy.sort()
            if writexy==1:writeXY(lxy,'%s.F_H%s.xy' %(nameout,prefw))
            lg.append(Gnuplot.Data(lxy, with_='lines lw 5', title='%s F{H5} %s' %(nameout,pref)))
        if cfep==1:
          zc=compZc(lx,dx,dt,eps=epszc)
          lxy=[(x,-math.log(zc[x])) for x in zc if zc[x]>0]
          lxy.sort()
          if writexy==1:writeXY(lxy,'%s.F_C%s.xy' %(nameout,prefw))
          lg.append(Gnuplot.Data(lxy, with_='lines lw 5', title='%s F_C %s' %(nameout,pref)))
        if cfepr!=0:
          zc=compZcr(lx,cfepr,dx,dt,eps=epszc)
          lxy=[(x,-math.log(zc[x]*2.)) for x in zc if zc[x]>0]
          lxy.sort()
          if writexy==1:writeXY(lxy,'%s.F_{C,%g}%s.xy' %(nameout,cfepr,prefw))
          lg.append(Gnuplot.Data(lxy, with_='lines lw 5', title='%s F_{C,1} %s' %(nameout,pref)))
        if cfep1==1:
          zc=compZcr(lx,1,dx,dt,eps=epszc)
          lxy=[(x,-math.log(zc[x])) for x in zc if zc[x]>0]
          lxy.sort()
          if writexy==1:writeXY(lxy,'%s.F_{C,1}%s.xy' %(nameout,prefw))
          lg.append(Gnuplot.Data(lxy, with_='lines lw 5', title='%s F_{C,1} %s' %(nameout,pref)))
        g.plot(*lg)
      
#-------------- other quantities ----------------------------------
      if D!=0:
        if D==1:lxy=compD(lx,dx,dt)
        if D==2:lxy=compD2(lx,dx,dt)
        lxy.sort()
        if writexy==1:writeXY(lxy,'%s.D.xy' %(nameout))
        lg.append(Gnuplot.Data(lxy, with_='lines lw 5', title='%s D' %nameout, axes= 'x1y2'))
        g('set y2tics \n set ytics nomirror')
        g('set y2label \'D\' ')
  
      if alpha!=None:
        dt1=int(alpha[0])
        dt2=int(alpha[1])
        lxy=compalpha(lx,dx,dt1,dt2)
        lxy.sort()
        if writexy==1:writeXY(lxy,'%s.alpha.xy' %(nameout))
        lg.append(Gnuplot.Data(lxy, with_='lines lw 5', title='%s alpha' %nameout, axes= 'x1y2'))
        g('set y2tics \n set ytics nomirror')
        g('set y2label \'alpha\' ')

      if pfold==1:
        lxy=comppfold(lx,dx,x0,x1)
        lxy.sort()
        if writexy==1:writeXY(lxy,'%s.pfold.xy' %(nameout))
        lg.append(Gnuplot.Data(lxy, with_='lines lw 5', title='%s pfold' %nameout, axes= 'x1y2'))
        g('set y2tics \n set ytics nomirror')
        g('set y2label \'pfold\' ')

      if pfoldMSM==1:
        dpf=comppfoldMSM(lx,dx,x0,x1)
        lxy=dpf.items()
        lxy.sort()
        if writexy==1:writeXY(lxy,'%s.pfoldMSM.xy' %(nameout))
        lg.append(Gnuplot.Data(lxy, with_='lines lw 5', title='%s pfoldMSM' %nameout, axes= 'x1y2'))
        g('set y2tics \n set ytics nomirror')
        g('set y2label \'pfold\' ')
     
      if ptpx==1:
        lxy=compptpx(lx,dx,x0,x1)
        lxy.sort()
        if writexy==1:writeXY(lxy,'%s.ptpx.xy' %(nameout))
        lg.append(Gnuplot.Data(lxy, with_='lines lw 5', title='%s p(TP|x)' %nameout, axes= 'x1y2'))
        lxy2=[(x,x*(1-x)*2) for x in [0.01*i for i in range(100)]]
        lg.append(Gnuplot.Data(lxy2, with_='lines lw 5', title='%s p(TP|x)_th' %nameout, axes= 'x1y2'))
        g('set y2tics \n set ytics nomirror')
        g('set y2label \'ptpx\' ')
  
  g('set xlabel \'%s\' ' %xlabel)
  g('set ylabel \'%s\' ' %ylabel)
  if glines:
    for l in glines:g(l)
  g.plot(*lg)
  if writeps==1: g.hardcopy(filename='%s.fep1d.eps' %nameout,color=1,fontsize=20)

#-------------- interactive analysis ----------------------------------
  if mfpt==1:comp_mfpt(lx,dx,ldt)
  if mfptq!=None:comp_mfpt_q(lx,dx,ldt,mfptq)

  raw_input("press key")  
#---END----------------------------------------------------------------

def parsecommandline(func,dbg=0):
  import sys
  def coerce(arg_value):
    try:
      return int(arg_value)
    except (TypeError,ValueError):
      pass
    try:
      return float(arg_value)
    except (TypeError,ValueError):
      pass
    if arg_value[0] in ('[', '{', '('):
      return eval(arg_value)
    return arg_value

  def parse_args(args):
    kw = {}
    pos = []
    for arg in args:
      if arg.startswith('--') and '=' in arg:
        name, value = arg.split('=', 1)
        kw[name[2:]] = coerce(value)
      else:
        pos.append(coerce(arg))
    return pos, kw
    
  argv=sys.argv
  pos,kw=parse_args(sys.argv[1:])
  print pos,kw
  if dbg==1:func(pos,**kw)  
  else:
    try:func(pos,**kw)  
    except TypeError , e:
      print e
      print func.__doc__

if __name__=="__main__":
  parsecommandline(fep1d,dbg=1)

