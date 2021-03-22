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

#repeated terms for r components:
def r1(alph, phi):
    return math.sqrt(1 + math.sin(float(alph))*math.sin(((math.pi)/6.0) - phi))
    
def r2(alph, phi): 
    return math.sqrt(1+math.sin(float(alph))*math.sin(((math.pi)/6.0) - phi)) 

def r3(alph, phi): 
    return math.sqrt(1-math.sin(float(alph))*math.cos(phi)) 
    
#Potentials: 
def Vy(alph): 
    return math.sqrt((3.0/2.0)*(1+abs(math.cos(alph)))) 

def Vdelta(alph, phi):
    return r1(alph,phi) + r2(alph,phi) + r3(alph,phi) 
    
def Vcoulomb(alph, phi): 
    return (1.0/r1(alph,phi)) + (1.0/r2(alph,phi)) + (1.0/r3(alph,phi)) 
        
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
    
#separate lines to split into columns (2D array)
for line in lines: 
    data = line.split()
    
    phi = math.atan2(float(data[5]),float(data[4]))
    #take subset of phi values  
    if phi == (math.pi/2.0): 
        #collect data for plot: 
        R = math.sqrt(float(data[0]))   #take square root of r2
        rlist.append(R) 
        alphlist.append(math.acos(float(data[6])/R))
        philist.append(phi)
        vlist.append(float(data[7]))
        
#lists of constants to be put into potential model 
phis = [] 
x = [] 
c = [] 
alphas = np.arange(0.01, math.pi, 0.01)

i = 0 
xvalue = 0.5 #Delta-Y string split. x = 1 corresponds to fully Y string. 
cvalue = 1.08764152 #constant in potential model 
phivalue = (math.pi/2.0)

#constants into lists of correct length
while i < len(alphas): 
    phis.append(phivalue)
    x.append(xvalue)
    c.append(cvalue) 
    i = i + 1

#weighted sum of potentials (weighting parameter x):
def Vsum(t): #t = (R,alph,phi,x,c)
    return t[0]*(t[2]*Vy(t[1]) + (1.0-t[3])*Vdelta(t[1],t[2])) - (Vcoulomb(t[1],t[2]))/t[0] + t[4]

fitresults = map(Vsum, list(zip(rlist,alphas,phis,x,c)))

fig = plt.figure() 

plt.scatter(alphlist, vlist)
plt.plot(alphas, fitresults)
plt.xlabel('cos(alpha)')
plt.ylabel('Potential, V') 

#******3D PLOT******
#ax = fig.add_subplot(111, projection = '3d')
#ax.scatter(philist, alphlist, vlist, zdir = 'z')
#ax.plot(phis, alphlist, fitresults, rstride = 4, cstride = 4, linewidth = 0) 
#ax.set_xlabel('phi')
#ax.set_ylabel('alpha')
#ax.set_zlabel('Potential, V')

plt.show() 