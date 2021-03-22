# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 11:27:38 2017

@author: jamesleech
"""
import math

#import file for hyper-radius, hyper-angles and side lengths 
f = open('r2_a_b_c_x_y_z_v', 'r') 
#split into rows - reading only equilateral triangles (0:48)
lines = f.readlines()[47:407]
f.close() 

#lists for values to filled from lines
rlist = [] 
vlist = [] 


#separate lines to split into columns (2D array)
for line in lines: 
    data = line.split()

    #Filter to take only results which match X,Y
    R = math.sqrt(float(data[0]))   #take square root of r2
    rlist.append(R)   
    vlist.append(float(data[7]))
    
    
print list(zip(rlist, vlist))
    
        