# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 15:38:59 2017

@author: jamesleech
"""

#3d plot of koma2 data points

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import csv,json,math
c=csv.reader(open("table3q.csv"),delimiter="&")

cs=[]
for i in c:
  cs.append({"i":eval(i[1]),"j":eval(i[2]),"k":eval(i[3]),
             "x1":float(i[4]),"x2":float(i[5]),"x3":float(i[6]),
             "v":float(i[8].replace("\\",""))})

#print json.dumps(cs, indent=4)

def perm(i,a,b,c):
  perm={0:(a,b,c),1:(a,c,b),2:(b,c,a),3:(b,a,c),4:(c,a,b),5:(c,b,a)}
  return perm[i]
  
xxlist = [] 
yylist = [] 
zzlist = [] 
rlist = []   
  
for i in cs:
  a1=i["x1"]
  b1=i["x2"]
  c1=i["x3"]
  for ix in range(6):
    a,b,c=perm(ix,a1,b1,c1)
    if a==0:
      x=0
      y=b
    else:
      x=(c*c+a*a-b*b)/a/2
      y=math.sqrt(c*c-x*x)
    rx=a/math.sqrt(2)
    ry=0
    lx=(a-2*x)/math.sqrt(6)
    ly=(-2*y)/math.sqrt(6)
  
    r2=rx*rx+ry*ry+lx*lx+ly*ly
    xx=2*(rx*lx+ry*ly)/r2
    yy=(rx*rx+ry*ry-lx*lx-ly*ly)/r2
    zz=2*(rx*ly-ry*lx)/r2  
  
    #print r2,a,b,c,xx,yy,zz,i["v"]
  
#add xx, yy, hyper-radius data to be plotted to lists: 
    xxlist.append(-yy)
    yylist.append(xx)
    zzlist.append(zz)
    rlist.append(math.sqrt(r2))


fig = plt.figure() 

ax = fig.add_subplot(111, projection = '3d')
ax.axis([-1.0,1.0,-1.0,1.0])

ax.scatter(xxlist, yylist, rlist, zdir = 'z')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('Potential, V')

plt.show() 