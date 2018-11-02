#!/usr/bin/env python
import Gnuplot
import os,time

g=Gnuplot.Gnuplot()
g('set style data lines')
g('set xrange [0:1]')
g('set y2tics nomirror')
g('set xlabel \'R\'')
g('set ylabel \'F/kT\'')
g('set y2label \'-ln Z_{C,1}/kT\'')
g('set key outside bottom spacing 1.5')

def plot(name):
  f=open(name)
  l=f.readline()
  f.close()
  l=l.split()
  str='plot \'%s\' using 1:2 title \'Z_H(r)\' with points axes x1y1, ' % name
  str=str+ '\'%s\' using 1:3 title \'Z_{C,-1}/2(r)\' with points  axes x1y1, ' % name 
  str=str+ '\'%s\' using 1:4 title \'Z_{C,1}(r,1)\' axes  x1y2 lt -1 lw 3, ' % name
  for i in range(5,len(l)+1): 
    str=str+' \'%s\' using 1:%i title \'Z_{C,1}(r,%i)\' axes x1y2, ' %(name,i,2**(i-4)) 
#  print str
  g(str)


name='rc.zc'
ts=os.path.getmtime(name)
plot(name)
while True:
  time.sleep(0.1)
  tsn=os.path.getmtime(name)
  if tsn!=ts:
    ts=tsn
    time.sleep(0.1)
    plot(name)
  
