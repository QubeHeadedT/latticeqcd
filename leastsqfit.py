# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 14:02:46 2017

@author: jamesleech
"""

#Performs a 3D fit of (A, sigY, C) using a function other than curve_fit 
#It uses leastsq from scipy library. 

from scipy.optimize import leastsq 
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#--------------------Potential Function Definitions:--------------------------
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
    
def Vcoul(A,x,y): 
    return A*((1.0/r1(x,y))+(1.0/r2(x,y))+(1.0/r3(x,y))) 
#-----------------------------------------------------------------------------

#--------------------Pull relevant data from data file:-----------------------
def readdata():
  f = open('r2_a_b_c_x_y_z_v', 'r')     #open file
  #Format - r2, a, b, c, x, y, z, v
  lines = f.readlines()[0:995] #Remove quark-diquark and linear ones
  f.close() 

  #make lists to be filled: 
  rlist = [] 
  xlist = [] 
  ylist = [] 
  vlist = [] 

  for line in lines: 
    data = line.split()
    R = np.sqrt(float(data[0])) 
    rlist.append(R) 
    xlist.append(-float(data[5]))     #coordinate flipping
    ylist.append(float(data[4]))      #coordinate flipping 
    vlist.append(float(data[7])) 
    
  return rlist, xlist, ylist, vlist 
#-----------------------------------------------------------------------------
    
#-----------------Produce 3D plot of (x, y, V):-------------------------------
def plot3D(): 

  #lists to be plotted:
  rlist = [] 
  xlist = [] 
  ylist = [] 
  vlist = [] 
  
  rlist, xlist, ylist, vlist = readdata()
  
  #-----Produce 3D plot (x,y,V*):------
  #Formatting: 
  fig = plt.figure() 
  ax = fig.add_subplot(111, projection = '3d')
  ax.set_xlabel('x')
  ax.set_ylabel('y')
  ax.set_zlabel('V')
  ax.axis([-1.0,1.0,-1.0,1.0])


  #Datapoints (x,y,Vstar)
  ax.scatter(xlist, ylist, vlist, zdir = 'z')
  
  plt.show() 
  
  return 
#-----------------------------------------------------------------------------
  
#-----------------Perform fitting for (x, y, V):------------------------------
def fitAYC(): 
  
  #Set up lists
  rlist = [] 
  xlist = [] 
  ylist = [] 
  vlist = [] 
  
  #Set list values from datafile 
  rlist, xlist, ylist, vlist = readdata() 
  
  #------------Perform fitting:----------------
  
  #just Y potential like in paper: 
  def Vfit1(r, coeff):   #r = [x, y] , coeff = [A, sigY, C]
    return -Vcoul(coeff[0], r[0], r[1])/r[2] + coeff[1]*Vy(r[0],r[1])*r[2] + coeff[2]  
  
  #define difference between fitted and actual values: 
  def residuals(coeff, y, r):
    res=[]
    for i in range(len(r)): 
      res.append(y[i]-Vfit1(r[i],coeff)) 
    return res
    
  #initial guess for parameters:
  coeff0 = np.array([0.1230, 0.1027, 0.9085], dtype = float) 
  
  #calculated fitted parameters: 
  x, flag = leastsq(residuals, coeff0, args = (vlist, list(zip(xlist, ylist,rlist)))) 
  
  #Print results: 
  print x
  print flag

  return; 

#-----------------------------------------------------------------------------
  
def AYCerror(): 
  
  #set up lists:
  rlist = [] 
  xlist = [] 
  ylist = [] 
  vlist = [] 
  
  #set list values from datafile:
  rlist, xlist, ylist, vlist = readdata() 
  
    #just Y potential like in paper: 
  def Vfit1(r, coeff):   #r = [x, y] , coeff = [A, sigY, C]
    return -Vcoul(coeff[0], r[0], r[1])/r[2] + coeff[1]*Vy(r[0],r[1])*r[2] + coeff[2]  
  
  #define difference between fitted and actual values: 
  def residuals(coeff, y, r):
    res=[]
    for i in range(len(r)): 
      res.append(y[i]-Vfit1(r[i],coeff))
    return res
    
  #initial guess for parameters:
  coeff0 = np.array([0.1230, 0.1027, 0.9085], dtype = float) 
  
  x, flag = leastsq(residuals, coeff0, args = (vlist, list(zip(xlist, ylist,rlist))))
  
  #caluclate error for each (x,y,V) value: 
  fiterrlist = map(abs, residuals(x, vlist, list(zip(xlist, ylist,rlist))))   
  
  #------caluculate average error value:-----

  print '[A, sigY, C] = ' + str(x) 
  
  av_err = sum(fiterrlist)/len(fiterrlist) 
  print 'Average fit error = ', av_err 
  
  return 
  
#-----------------Calculate the average AYC fitting error:--------------------
  
#-----------------------------------------------------------------------------



  
  
#=========================Run-time code:======================================
  
AYCerror() 
