# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 14:57:15 2017

@author: jamesleech
"""

#Script to produce a 3D plot of (x, y, V) from an old koma datafile

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#------------READ DATA FROM FILE:--------------

#-------Select from below:----------
#f = open('data57', 'r') 
#f = open('data58', 'r')
f = open('data60', 'r')
#f = open('r2_a_b_c_x_y_z_v', 'r') 
#-----------------------------------


#------Read x, y, v from the data file in plottable lists:---------------
lines = f.readlines()   
f.close() 
        
#lists for values to filled from lines
xlist = [] 
ylist = [] 
vlist = [] 
        
#separate lines to split into columns (2D array)
for line in lines: 
  data = line.split() 
  ylist.append(float(data[4]))   #coordinate flipping
  xlist.append(-float(data[5]))    #coordinate flipping 
  vlist.append(float(data[7]))
#----------------------------------------------------------------------
  
#--------------Plot results:-------------------------------------------
  
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


