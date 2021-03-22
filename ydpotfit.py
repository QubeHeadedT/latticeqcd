# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 09:52:59 2017

@author: jamesleech
"""

import math 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import numpy as np 

#Plots polar angle vs potential for all 220 results.
#Compares a weighted sum of delta and y-string potential to results

#------define functions:-------------

#Normalisation constants to make results match constantcheck.py fit values
Ak = 8.0311785558
Bk = 0.07637747
#^^^^ Maybe not needed? ^^^^

#repeated terms for r components:
def r1(alph, phi):
    return math.sqrt(1 + math.sin(float(alph))*math.sin(((math.pi)/6.0) - phi))
    
def r2(alph, phi): 
    return math.sqrt(1+math.sin(float(alph))*math.sin(((math.pi)/6.0) + phi)) 

def r3(alph, phi): 
    return math.sqrt(1-math.sin(float(alph))*math.cos(phi)) 
    
#Potentials: 
def Vy(alph): 
    return math.sqrt((3.0/2.0)*(1+abs(math.cos(alph)))) 

def Vdelta(alph, phi):
    return (r1(alph,phi) + r2(alph,phi) + r3(alph,phi)) 
    
def Vcoulomb(alph, phi): 
    return Ak*(1.0/r1(alph,phi)) + (1.0/r2(alph,phi)) + (1.0/r3(alph,phi)) 

def B(alph,phi,x):
    return (Bk/(x*(math.sqrt(3)-3)+3))*(x*Vy(phi) + (1.0-x)*Vdelta(alph,phi)) 
        
#------------------------------------
        
#import file for hyper-radius, hyper-angles and side lengths 
f = open('r2_a_b_c_x_y_z_v', 'r') 
#split into rows - reading only equilateral triangles (0:48)
lines = f.readlines()
f.close() 


alphlist = [] #polar angle
philist = [] #azimuthal angle 
rlist = [] #hyper-radius 
vlist = [] #potential
fitresults = [] #potential as predicted by model 
y = [] 

    
#separate lines to split into columns (2D array)
for line in lines: 
    data = line.split()
    if float(data[5]) == 0.0 and float(data[4]) == 0.0:
        phi = math.pi/2.0
        print 'hello'
    else:
        phi = math.atan2(float(data[5]),float(data[4]))
    #take subset of phi values  
    if phi == (math.pi/2.0): 
        #collect data for plot: 
        R = math.sqrt(float(data[0]))   #take square root of r2
        rlist.append(R) 
        alphlist.append(math.acos(-float(data[6])))
        philist.append(phi)
        vlist.append(float(data[7]))
        y.append(data[5])
        
        
#lists of constants to be put into potential model 
phis = [] 
x = [] 
c = [] 
rs = [] 
alphas = np.arange(0.0, math.pi, 0.01)

i = 0 
xvalue = 0.5 #Delta-Y string split. x = 1 corresponds to fully Y string. 
cvalue = 1.08764152 #constant in potential model 
rvalue = 2.0    #arbitrary r value 
phivalue = 0.0   

#----------Calculate a fit to a mix of delta and y string potential------------
#constants into lists of correct length:
while i < len(alphas): 
    phis.append(phivalue)
    x.append(xvalue)
    c.append(cvalue) 
    rs.append(rvalue)
    i = i + 1

#weighted sum of potentials (weighting parameter x):
def Vsum(t): #t = (R,alph,phi,x,c)
    return (t[0]*(B(t[1],t[2],t[3])) 
        - (Vcoulomb(t[1],t[2]))/t[0] + t[4])

#Calculate results
fitresults = map(Vsum, list(zip(rs,alphas,phis,x,c)))


#plt.plot(map(math.cos,alphas), fitresults)
#---------------2D PLOT --------------------------------
#fig = plt.figure() 

plt.scatter(alphlist, vlist)
plt.plot(alphas, fitresults)
plt.xlabel('alpha')
plt.ylabel('Potential, V')

#--------------- 3D PLOT -------------------------------
#ax = fig.add_subplot(111, projection = '3d')
#ax.scatter(philist, map(math.cos,alphlist), vlist, zdir = 'z')
#ax.plot(phis, alphlist, fitresults, rstride = 4, cstride = 4, linewidth = 0) 
#ax.set_xlabel('phi')
#ax.set_ylabel('alpha')
#ax.set_zlabel('Potential, V')

plt.show() 
