# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 13:27:16 2017

@author: jamesleech
"""

#Script for fitting a surface to Vstar potential with a linear split
#between delta and y string potential. 

import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit 
import numpy as np 
from pylab import meshgrid

#-------------Select the dataset from below:----------------
#dataset = 'data58' 
#dataset = 'data60'
dataset = 'r2_a_b_c_x_y_z_v'
#-----------------------------------------------------------

#-------------Function Definitions:------------------------------------
#Calculates Vstar:
def Vstar(v,k,a,b,c,D,R):   #v potential, R hyper-radius, k & D fitting consts.
    return (v - (1/k)*((1/a)+(1/b)+(1/c)) - D)/R 

#calculates 2D radius from centre of circle at origin
def Radius(x,y): 
    return x*x+y*y
#Definitions for delta and y string potentials in terms of (x,y): 
def sinacos(x,y):   #a repeated term
    return np.sin(np.arccos(np.sqrt(1.0-x*x-y*y))) 

#functions to calculate distance components:
def r1(x,y): 
    return np.sqrt(1.0+sinacos(x,y)*np.sin((np.pi/6.0)-np.arctan2(y,x))) 
    
def r2(x,y): 
    return np.sqrt(1.0+sinacos(x,y)*np.sin((np.pi/6.0)+np.arctan2(y,x)))

def r3(x,y): 
    return np.sqrt(1.0 - sinacos(x,y)*np.cos(np.arctan2(y,x))) 
    
#String and coulomb potentials: 
def Vy(x,y): 
    return np.sqrt((3.0/2.0)*(1.0+abs(np.sqrt(1-x*x-y*y))))
    
def Vdelta(x,y): 
    return r1(x,y) + r2(x,y) + r3(x,y) 
    
def Vcoul(x,y): 
    return k*((1.0/r1(x,y))+(1.0/r2(x,y))+(1.0/r3(x,y))) 
    
#Potential w/ Linear combination of delta and Vy:  ***To be fitted!***
def Vlin(r, A, B):   #t balance parameter. 0 => all delta, 1 => all y string
    return A*Vy(r[0],r[1])+(B)*Vdelta(r[0],r[1]) 

#----------------------------------------------------------------------

#------------Constant Definitions:-------------------------------------
#Constants from least-squares fit, change for each dataset:
if dataset == 'data60':
  k = -8.5154214306
  D = 0.922169639534
if dataset == 'data58':
  k = -14.7780882602
  D = 0.787366529133
if dataset == 'r2_a_b_c_x_y_z_v':
  k = -8.0311785558
  D = 1.0876415217

#Select string to modify graph title to fit data:
if dataset == 'data58':
  datastring = 'Takahashi data, beta = 5.8'
if dataset == 'data60':
  datastring = 'Takahashi data, beta = 6.0'
if dataset == 'r2_a_b_c_x_y_z_v':
  datastring = 'Koma & Koma data, beta = 6.0'

#--------------------------------------------------------------------


#-------------Sub-routine to read x,y from file from maxr radius-------
#--------------and return best fit params and covariances:-------------
def rmaxAB(maxr): 

    
    f = open(dataset, 'r')

    if dataset == 'r2_a_b_c_x_y_z_v':
    #Don't read quark-diqark collions in Koma data - avoid division by zero [0:1187]
      lines = f.readlines()[0:1187]
    else:
      lines = f.readlines()
      
    f.close() 
        
    #lists for values to filled from lines
    rlist = []  
    xlist = [] 
    ylist = [] 
    vstarlist = [] 
        
        #separate lines to split into columns (2D array)
    for line in lines: 
        data = line.split()
        R = math.sqrt(float(data[0]))   #take square root of r2 for hyper-rad.
        #Filter results within rmax: 
        if (Radius(float(data[4]),float(data[5])) < maxr) and R > HRmin: 
            rlist.append(R)   
            ylist.append(float(data[4]))   #coordinate flipping
            xlist.append(-float(data[5]))    #coordinate flipping 
            vstarlist.append(Vstar(float(data[7]),k,float(data[1]),
                                   float(data[2]),
                                    float(data[3]),D,R))
        
    #calculate best fit params.        
    popt, pocv = curve_fit(Vlin,[xlist,ylist], vstarlist)  
        
    return popt, pocv         
#-------------------------------------------------------------------------
        
#------------Subroutine to write A/B for range of rmax to file:-------------
def write_rmaxAB(minr,maxr, numr): #start, stop and # of steps of rmax
#NOTE: avoid minr = 0.0 
    
    rs = np.linspace(minr, maxr, numr) 
    
    f = open('rmaxAB.txt','w') 
    
    for r in rs: 
        popt, pocv = rmaxAB(r)  #calculate A,B parameters
        AdivB = popt[0]/popt[1]    #A/B
        f.write(str(r) + ' ' + str(AdivB) + '\n') 
        
    f.close() 
    return; 
#-------------------------------------------------------------------------
    
#---------Subroutine to read data from rmaxAB and produce a plot:---------
def plot_rmaxAB(minr, maxr, numr): #start, stop and #of steps of rmax 
#NOTE: avoid minr = 0.0 
    
    #produce file to be plotted: 
    write_rmaxAB(minr, maxr, numr)    
    
    f = open('rmaxAB.txt', 'r') 
    
    lines = f.readlines()   
    f.close() 
    
    #lists to be plotted: 
    r = [] 
    AdivB = [] 
    
    #read data from file and fill out lists
    for line in lines: 
            data = line.split()
            r.append(data[0])
            AdivB.append(float(data[1]))
            
    #plot lists:
    plt.plot(r, AdivB) 
    plt.xlabel('maximum fitting radius')
    plt.ylabel('sigY/sigD')
    plt.title('Ratio sigY/sigD vs. maximum fitting radius \n' + datastring + ', #steps = ' + str(numr))
    plt.show()
    
    return;    
#-------------------------------------------------------------------------
    
#--------------Surface fit plot subroutine:-------------------------------
    
def surfaceplot(maxr):
    
    #------------READ DATA FROM FILE:--------------
    f = open(dataset, 'r') 

    lines = f.readlines()   
    f.close() 
        
    #lists for values to filled from lines
    rlist = []  
    xlist = [] 
    ylist = [] 
    vstarlist = [] 
        
        #separate lines to split into columns (2D array)
    for line in lines: 
        data = line.split()
        R = math.sqrt(float(data[0]))   #take square root of r2 for hyper-rad.
        #Filter results within rmax: 
        if (Radius(float(data[4]),float(data[5])) < maxr) and R > HRmin: 
            rlist.append(R)   
            ylist.append(float(data[4]))   #coordinate flipping
            xlist.append(-float(data[5]))    #coordinate flipping 
            vstarlist.append(Vstar(float(data[7]),k,float(data[1]),
                                   float(data[2]),
                                    float(data[3]),D,R))
        
        
    #----------Calculate best fit params:---------------     
    popt, pocv = curve_fit(Vlin,[xlist,ylist], vstarlist)
    
    #-----Calculate Curve Fit Surface:-----
    print 'A =', popt[0]
    print 'B =', popt[1]
    print pocv 

    fitdata = [] 

    #range of (x,y) surface - ***Make into a circle?***:
    xrang = np.arange(-rmax, rmax, 0.005)
    yrang = np.arange(-rmax, rmax, 0.005)
    X, Y = meshgrid(xrang,yrang)
    
    #apply fitted function to surface: 
    fitdata = Vlin([X,Y],popt[0], popt[1])
    
    #-----Produce 3D plot (x,y,V*):------
    #Formatting: 
    fig = plt.figure() 
    ax = fig.add_subplot(111, projection = '3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('Vstar')
    plt.title('x, y, V*' + '\n' + datastring + '\n' + 'rmax = ' + str(maxr) + ', sigY = ' + str(popt[0]) + ', sigD = ' + str(popt[1]))
    ax.axis([-1.0,1.0,-1.0,1.0])


    #Datapoints (x,y,Vstar)
    ax.scatter(xlist, ylist, vstarlist, zdir = 'z')

    #Fitted surface: 
    #ax.plot_surface(X, Y, fitdata, color = 'red')

    plt.show() 
      
    return;
    
#-------------------------------------------------------------------------

#--------------------------Code to be run:--------------------------------

#---------Constants:-------
#parameter to limit radius of data analysed (0.0 -> 1.0): 
rmax = 1.0

#Minimum hyper-radius to be considered (coulomb consideration):  
HRmin = 0.0
#--------------------------

#-----rmax vs A/B:------
plot_rmaxAB(0.01, rmax, 500)

#-----Fitted surface to data plot:-----
#surfaceplot(rmax) 

#-------------------------------------------------------------------------
