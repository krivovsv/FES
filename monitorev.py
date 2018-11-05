#!/usr/bin/env python
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
import os,time

plt.ion()
fig=plt.figure(figsize=(10,6))

def plot(name):
  plt.clf()
  ax1=fig.add_subplot(111)
  ax1.set_ylabel('$\mathrm{F/kT}$')
  ax1.set_xlabel('r')
  ax2=ax1.twinx()
  ax2.set_ylabel('$-\ln \mathrm{Z_C/Z_C}$')
  f=open(name)
  ls=[]
  for l in f:
    l=l.split()
    l=[float(s) for s in l]
    ls.append(l)
  f.close()
  lx=[l[0] for l in ls]
  ax1.plot(lx,[l[1] for l in ls], 'xb',label='$-\ln \mathrm{Z_H}(r)$')
  ax1.plot(lx,[l[2] for l in ls], '+r',label='$-\ln \mathrm{2Z_{C,-1}}(r)$')
  ax2.plot(lx,[l[3] for l in ls], '-k',label='$-\ln \mathrm{Z_C/Z_C}(r,1)$')
  ax2.plot(lx,[l[4] for l in ls], '-', color= '0.5',label='$-\ln \mathrm{Z_C/Z_C}(r,dt>1)$')
  for i in range(5,len(ls[0])):
    ax2.plot(lx,[l[i] for l in ls], '-', color='0.5')
#  l1,lab1=ax1.get_legend_handels_labels()
#  l2,lab2=ax2.get_legend_handels_labels()
  fig.legend(loc='upper center',mode='expand',ncol=4)
  fig.canvas.draw()
  fig.canvas.flush_events()
  plt.draw()

name='ev.zc'
ts=os.path.getmtime(name)
plot(name)
while True:
  time.sleep(0.1)
  tsn=os.path.getmtime(name)
  if tsn!=ts:
    ts=tsn
    time.sleep(0.1)
    plot(name)
  
