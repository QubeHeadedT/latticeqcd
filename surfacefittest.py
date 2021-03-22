# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 15:55:43 2017

@author: jamesleech
"""

import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit 
import random
from pylab import meshgrid 

#surface fitting practise with an arbitrary data set 

#------Function to be fitted:
def func(r, a, b):  
    return a*r[0] + b*r[1]


#------Random data set with noise:
r = np.array([[1.0, 2.0, 3.0, 4.0, 5.0],[2.0, 3.0, 1.0, 5.0, 4.0]])
zfunc = func(r,1,2)
znoise = [random.random(),random.random(),random.random(),
          random.random(),random.random()] 
z = zfunc + znoise 

#------Find curve fit: 
popt, pocv = curve_fit(func,r,z) 

fitdata = [] 

xrang = np.arange(0,5.0,0.05)
yrang = np.arange(0,5.0,0.05)

X, Y = meshgrid(xrang,yrang)

fitdata = func([X,Y],popt[0],popt[1])

#-----Plot:
fig = plt.figure() 
ax = fig.add_subplot(111, projection = '3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z') 

#Plot data: 
ax.scatter(r[0], r[1], z, zdir = 'z')

#Plot curve fit:
ax.plot_surface(X, Y, fitdata, color = 'red')

plt.show() 