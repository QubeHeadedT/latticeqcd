# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 15:41:06 2017

@author: jamesleech
"""

import matplotlib.pylab as plt, math
from scipy.optimize import curve_fit
import numpy

#Script to calculate the constant K of coloumb term for pillars of results
#using a fit curve. 

#import file for hyper-radius, hyper-angles and side lengths 
f = open('r2_a_b_c_x_y_z_v', 'r') 
#split into rows - reading only equilateral triangles (0:48)
lines = f.readlines()
f.close() 

#lists for values to filled from lines
rlist = [] 
alist = [] 
blist = [] 
clist = [] 
xlist = [] 
ylist = [] 
zlist = [] 
vlist = [] 
AnaA = []   #Analytical A value: (1/a + 1/b + 1/c)*R

#Select which measurements to consider based on X,Y hyperangles
X = 0.0
Y = 0.
err = 1e-15   #to account for computing precision errors on zeros

#separate lines to split into columns (2D array)
for line in lines: 
    data = line.split()
    
    #Filter to take only results which match X,Y
    if (abs(float(data[4]) - X) < err) and (abs(float(data[5]) - Y) < err):
        R = math.sqrt(float(data[0]))   #take square root of r2
        rlist.append(R)   
        alist.append(float(data[1]))
        blist.append(float(data[2]))
        clist.append(float(data[3]))
        xlist.append(-float(data[5]))
        ylist.append(float(data[4]))
        zlist.append(float(data[6]))
        vlist.append(float(data[7]))
        abc = (1/float(data[1])) + (1/float(data[2])) + (1/float(data[3]))
        AnaA.append(abc * R)

#plot R against V:
plt.scatter(rlist,vlist, c = 'red', marker = 'x')
plt.xlabel('Hyper-radius')
plt.ylabel('Potential, V')

#Take the Average analytical A
AveanaA = sum(AnaA)/float(len(AnaA))
print 'Average A analytical = ', AveanaA

#fit potential model to curve - using least squares fit - define V function: 
def potential(R, A, B, C):
    return (A/R) + B*R + C 
    
#optimise parameters
popt, pcov = curve_fit(potential, rlist, vlist)
print 'A fitted =', popt[0]
print 'B fitted =', popt[1]
 

#produce list for fit curve using parameters calculated
#fitresult = [] 
#for r in rlist: 
 #   fitresult.append((popt[0]/r) + popt[1]*r + popt[2])
#plt.plot(rlist, fitresult)

#produce curve using calculated parameters
x = numpy.linspace(0,12,100)
y = potential(x, popt[0],popt[1],popt[2])
plt.plot(x,y)
plt.axis([0, 12, 0,2])
plt.show() 

#calculate K = Aanalytical/Afit
print 'Ratio =', abs(AveanaA/popt[0])

k = AveanaA/popt[0]
print k

secondres = [] 

def Veff(x):
    return (x[4] - (1./k)*((1/x[1]) + (1/x[2]) + (1/x[3])) - popt[2])/x[0]
    

secondres = map(Veff, list(zip(rlist, alist, blist, clist, vlist)))

plt.plot(rlist, secondres)
#plt.show() 


    
    





