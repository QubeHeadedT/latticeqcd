# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:06:07 2017

@author: jamesleech
"""

import matplotlib.pylab as plt

#Produces a circle plot in x,y plane of data points taken, 
#color-coded by hyper-radius squared

#Open file (select from comments) and read lines:
#------------------------
f = open('r2_a_b_c_x_y_z_v', 'r')

#Format of datafiles - r2, a, b, c, x, y, z, v
#------------------------

lines = f.readlines()[0:995] 
f.close() 

#lists of data to be plotted:
xlist = [] 
ylist = [] 
r2list = [] 

#strip useful data from f as floats
for line in lines: 
    data = line.split()
    xlist.append(float(data[4]))
    ylist.append(-float(data[5]))
    #note x, y coordinate flipping above for r3 term -> 0
    r2list.append(float(data[0])) 
     
#----------Plot the data:------------------------
    
plt.scatter(xlist, ylist, c = r2list, cmap = 'jet', marker = 'x') 

#formatting:
plt.axis([-1.0,1.0,-1.0,1.0])
plt.xlabel('x')
plt.ylabel('y')
plt.title('Koma data')
plt.colorbar(mappable=None, cax=None, ax=None, label = 'Potential, V')

plt.show() 
#------------------------------------------------
