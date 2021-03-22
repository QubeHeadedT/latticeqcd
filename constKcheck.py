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

#-------------Select one from below:-----------------------------
dataset = 'data58'
dataset = 'data60'
dataset = 'r2_a_b_c_x_y_z_v'
#----------------------------------------------------------------

#-----------Modify title based on datafile selected:-------------

if dataset == 'data58':
  datastring = 'Takahashi, beta = 5.8'
if dataset == 'data60':
  datastring = 'Takahashi, beta = 6.0'
if dataset == 'r2_a_b_c_x_y_z_v':
  datastring = 'Koma data, beta = 6.0'


#----------------------------------------------------------------

#---------------Read data from file:-----------------------------
f = open(dataset, 'r')
#split into rows 
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
Y = 0.0
err = 1e-3 #to account for computing precision errors on zeros

#separate lines to split into columns (2D array)
for line in lines: 
    data = line.split()

    x = -float(data[5]) #coordinate flipping
    y = float(data[4]) #coordinate flipping
    
    #Filter to take only results which match X,Y
    if abs(x - X) < err and abs(float(y - Y)) < err:
        R = math.sqrt(float(data[0]))   #take square root of r2
        rlist.append(R)   
        alist.append(float(data[1]))
        blist.append(float(data[2]))
        clist.append(float(data[3]))
        xlist.append(x) #coordinate flipping
        ylist.append(y)  #coordinate flipping
        zlist.append(float(data[6]))
        vlist.append(float(data[7]))
        abc = (1/float(data[1])) + (1/float(data[2])) + (1/float(data[3]))
        AnaA.append(abc * R)      
#--------------------------------------------------------------------------

#plot R against V:
plt.scatter(rlist,vlist, c = 'red', marker = 'x')
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
print 'C fitted =', popt[2] 
 
#produce curve using calculated parameters
x = numpy.linspace(0,12,100)
y = potential(x, popt[0],popt[1],popt[2])
plt.plot(x,y)
plt.axis([0, 12, 0,4])

#calculate and print K = Aanalytical/Afit
k = AveanaA/popt[0]
print 'k =', abs(k)

#Formatting:
plt.xlabel('Hyper-radius')
plt.ylabel('Potential, V')
plt.title('Hyper-radius vs. Potential V at (x,y) = ' + '(' + str(X) + ', ' + str(Y) + ')' + '\n' + datastring)
plt.text(6.0, 0.5, 'Analytical A = ' + str(AveanaA) + '\n' + 'Ratio K = ' + str(k) + '\n' + 'Error in (x,y) values = ' + str(err))
plt.text(2.0, 3.0, 'Fitted Values:' + '\n' + 'A = ' + str(popt[0]) + '\n' + 'B = ' + str(popt[1]) + '\n' + 'C = ' + str(popt[2]))

plt.show() 


  
    





