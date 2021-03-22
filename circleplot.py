# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 10:47:44 2017

@author: jamesleech
"""

#***Adapted from Milovan Suvakov's code:***
import csv,json,math, numpy as np, matplotlib.pylab as plt, matplotlib.colors as colors
c=csv.reader(open("table3q.csv"),delimiter="&")

#x,y and hyperradius^2 r2 to be used for plot
xxlist = []
yylist = []
r2list = [] 

cs=[]
for i in c:
  cs.append({"i":eval(i[1]),"j":eval(i[2]),"k":eval(i[3]),"x1":float(i[4]),"x2":float(i[5]),"x3":float(i[6]),"v":float(i[8].replace("\\",""))})

#print json.dumps(cs, indent=4)

#for all 6 permutations 
def perm(i,a,b,c):
  perm={0:(a,b,c),1:(a,c,b),2:(b,c,a),3:(b,a,c),4:(c,a,b),5:(c,b,a)}
  return perm[i]
  
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
    
    #add xx, yy, hyper-radius data to be plotted to lists: 
    yylist.append(xx)
    xxlist.append(-yy)
    r2list.append(r2)    

#Plot of xx, yy colorcoded by hyper-radius squared:
plt.scatter(xxlist,yylist,c = r2list, cmap = 'jet', marker = 'x')

#formatting:
plt.axis([-1.0,1.0,-1.0,1.0])
plt.xlabel('x')
plt.ylabel('y')
plt.title()
plt.colorbar(mappable=None, cax=None, ax=None, label = 'Hyper-radius squared')

plt.show() 


