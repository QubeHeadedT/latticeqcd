# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 13:18:17 2017

@author: jamesleech
"""

import math 
import matplotlib.pylab as plt 
import numpy as np 
from scipy.optimize import curve_fit


#normalisation constant : 

T = (3.8/75.0) #for y0 line
#T = 2.55/10.0 #for x5 line

#------------------Potential Function Definitions:--------------------
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
def Vy(x,y,R): 
    return (T)*np.sqrt(3)*(1.0/R)*(np.sqrt((3.0/2.0)*
            (1.0+abs(np.sqrt(1-x*x-y*y)))))
    
def Vdelta(x,y,R): 
    return (T)*(1.0/R)*(r1(x,y) + r2(x,y) + r3(x,y))
    
def Vstar(v,k,a,b,c,D,R):   #v potential, R hyper-radius, k & D fitting consts.
    return (v - (1/k)*((1/a)+(1/b)+(1/c)) - D)/R 
    
#Potential w/ Linear combination of delta and Vy:  ***To be fitted!***
def Vlin(r, A, B):   #t balance parameter. 0 => all delta, 1 => all y string
    return A*Vy(r[0],r[1], r[2])+(B)*Vdelta(r[0],r[1], r[2]) 
    #r = [x, y, R] 
    
#--------------------------------------------------------------------

#---------------------------Constants:-------------------------------
#Constants from least-squares fit:
k = -8.0311785558   #Aana/Afit 
D = 1.08764152      #Const.
#--------------------------------------------------------------------

#-----------------Subroutine to produce a circle plot:---------------
def circleplot():
    
    #import file for hyper-radius, hyper-angles and side lengths 
    f = open('r2_a_b_c_x_y_z_v', 'r') 

    lines = f.readlines()
    f.close() 

    #lists for values to filled from lines
    rlist = []  
    xlist = [] 
    ylist = []
    alist = [] 
    blist = [] 
    clist = [] 
    vlist = [] 

    for line in lines: 
        data = line.split()
        R = math.sqrt(float(data[0]))   #take square root of r2
        x = -float(data[5])         #coordinate flipping
        y = float(data[4])      #coordinate flipping 
        rlist.append(R)   
        xlist.append(x)
        ylist.append(y)
        alist.append(float(data[1]))
        blist.append(float(data[2]))
        clist.append(float(data[3]))
        vlist.append(float(data[7]))
        
    plotscatter(xlist,ylist, 'g')
    return 

#------------------------------------------------------------------

#-----------Sub-routine to plot scatter given xlist and ylist:-----
def plotscatter(xlist,ylist, color):
    plt.scatter(xlist,ylist, c = color)
    plt.xlabel('x')
    plt.ylabel('y') 
    return

#--------Subroutine to isolate y = 0 line--------------------------
def y0line():
    
    #import file for hyper-radius, hyper-angles and side lengths 
    f = open('r2_a_b_c_x_y_z_v', 'r') 

    lines = f.readlines()
    f.close() 
    
    err = 1e-15 #accounts for numerical error in zeroes

    #lists for values to filled from lines
    rlist = []  
    xlist = [] 
    ylist = []
    alist = [] 
    blist = [] 
    clist = [] 
    vlist = [] 

    for line in lines: 
        data = line.split()
        
        R = math.sqrt(float(data[0]))   #take square root of r2
        x = -float(data[5])             #coordinate flipping
        y = float(data[4])              #coordinate flipping 

        #Filter only x = 0 line: 
        if abs(y) < err and x > -0.5:        
            rlist.append(R)   
            xlist.append(x)
            ylist.append(y)
            alist.append(float(data[1]))
            blist.append(float(data[2]))
            clist.append(float(data[3]))
            vlist.append(float(data[7]))
    
    return xlist, ylist, alist, blist, clist, rlist, vlist 

#-------------------------------------------------------------------

#---------Subroutine to isolate x = -0.5 line:----------------------
def x5line():
    
    #import file for hyper-radius, hyper-angles and side lengths 
    f = open('r2_a_b_c_x_y_z_v', 'r') 

    lines = f.readlines()
    f.close() 
    
    err = 1e-3 #accounts for numerical precision error 

    #lists for values to filled from lines
    rlist = []  
    xlist = [] 
    ylist = []
    alist = [] 
    blist = [] 
    clist = [] 
    vlist = [] 

    for line in lines: 
        data = line.split()
        
        R = math.sqrt(float(data[0]))   #take square root of r2
        x = -float(data[5])             #coordinate flipping
        y = float(data[4])              #coordinate flipping 

        #Filter only x = 0 line: 
        if abs(x + 0.5) < err:        
            rlist.append(R)   
            xlist.append(x)
            ylist.append(y)
            alist.append(float(data[1]))
            blist.append(float(data[2]))
            clist.append(float(data[3]))
            vlist.append(float(data[7]))
    
    return xlist, ylist, alist, blist, clist, rlist, vlist 

#-------------------------------------------------------------------

#--------Subroutine filter lists with range of Hyper-radius:--------
def binsplit(xlist, ylist, alist, blist, clist, rlist, vlist, HRmin, HRmax) :
        
    data = list(zip(xlist, ylist, alist, blist, clist, rlist, vlist)) 
    
    xbinlist = [] 
    ybinlist = [] 
    abinlist = [] 
    bbinlist = [] 
    cbinlist = [] 
    rbinlist = []
    vbinlist = [] 
    
    for line in data: 
        if HRmin < line[5] < HRmax:
            xbinlist.append(line[0])
            ybinlist.append(line[1])
            abinlist.append(line[2])
            bbinlist.append(line[3])
            cbinlist.append(line[4])
            rbinlist.append(line[5])
            vbinlist.append(line[6]) 
            
    return xbinlist, ybinlist, abinlist, bbinlist, cbinlist, rbinlist, vbinlist 
#-------------------------------------------------------------------

#---------Sub-routine to calculate the potentials for a point-------
def pointpot(x, y, a, b, c, v, R): 
    
    valVstar = Vstar(v,k,a,b,c,D,R)   
    valV = v/R
    
    return valVstar, valV
#-------------------------------------------------------------------
    
#---------Sub-routine to calculate Vd and Vy for given (x,y):-------
def dypot(x,y,R):
    
    valVd = Vdelta(x,y,R)
    valVy = Vy(x,y,R) 
    
    return valVd, valVy
#-------------------------------------------------------------------
    
#-----------Function to plot y0 line graph--------------------------
def ploty0(minHR, maxHR):
    
    #isolate line
    xlist, ylist, alist, blist, clist, rlist, vlist = y0line() 
    
    #isolate particular HR bin 
    xlist, ylist, alist, blist, clist, rlist, vlist = binsplit(xlist, 
                        ylist, alist, blist, clist, rlist, vlist, minHR, maxHR)
            
    #Scatter plot of vstar and v/HR for data:          
    params = list(zip(xlist, ylist, alist, blist, clist, vlist, rlist)) 
    vstarlist = [] 
    v_rlist = [] 
    
    #Calculate the potentials for the points
    for line in params: 
        if line[2] != 0 and line[3] != 0 and line[4] != 0:  
            vstar, v_r = pointpot(line[0], line[1], line[2], line[3],
                               line[4], line[5], line[6])
            vstarlist.append(vstar)
            v_rlist.append(v_r)
        else: 
            #When a, b, or c == 0: 
            vstarlist.append(0.0)
            v_rlist.append(0.0)
        
    
    #vary x for this line - equivalent to polar coordinate:
    plt.scatter(xlist, vstarlist, c = 'g')
    #plt.scatter(xlist, v_rlist, c = 'r')
    #------------------------------------------
    
    #add a curve fit line for linear combination of Y and Delta: 
    #popt, pocv = popt, pcov = curve_fit(Vlin, [xlist,np.zeros(len(xlist))],
                                               #vstarlist)
    #print coefficients: 
    #print 'A = ', popt[0]
   # print 'B = ', popt[1]     
    
    #Formatting: 
    plt.xlabel('y')
    plt.ylabel('Potential, V')
    plt.text(0.0,0.06,'T ratio = ' + str(T) + '\n' + 'T*Vdelta, T*Sqrt(3)*Vy')
    plt.title('Vstar Isosceles Triangle, y = 0 line' + '\n' + 'Blue: Delta string, Green: Y string (Koma data)')
    
    return 

              
#---------Same as above but for x =-0.5 line:--------------------------
    
def plotx5(minHR, maxHR):
        
         #isolate line
    xlist, ylist, alist, blist, clist, rlist, vlist = x5line() 
    
    #isolate particular HR bin 
    xlist, ylist, alist, blist, clist, rlist, vlist = binsplit(xlist, 
                        ylist, alist, blist, clist, rlist, vlist, minHR, maxHR)
                        
    #Scatter plot of vstar and v/HR for data:          
    params = list(zip(xlist, ylist, alist, blist, clist, vlist, rlist)) 
    vstarlist = [] 
    v_rlist = [] 
    
    #Calculate the potentials for the points
    for line in params: 
        if line[2] != 0 and line[3] != 0 and line[4] != 0:  
            vstar, v_r = pointpot(line[0], line[1], line[2], line[3],
                               line[4], line[5], line[6])
            vstarlist.append(vstar)
            v_rlist.append(v_r)
        else: 
            #When a, b, or c == 0: 
            vstarlist.append(0.0)
            v_rlist.append(0.0)
        
    
    #vary x for this line - equivalent to polar coordinate:
    plt.scatter(ylist, vstarlist, c = 'g')
    #plt.scatter(ylist, v_rlist, c = 'r')
    
    #Formatting: 
    plt.xlabel('y')
    plt.ylabel('Potential, V')
    plt.text(0.0,0.06,'T ratio = ' + str(T) + '\n' + 'T*Vdelta, T*Sqrt(3)*Vy')
    plt.title('Vstar Right Angle Triangles - x = -0.5 line' + '\n' + 'Blue: Delta string, Green: Y string (Koma data)')

    return 
#-------------------------------------------------------------------
        
    
#--Subroutine to plot lines for delta and y string potentials along y0 line:---
    
def dylinesy0(R):  
    
    x = np.linspace(-0.5, 1.0, 200) 
    
    dlist = []  #delta potentials 
    ylist = []  #y potentials 
    
    dlist = Vdelta(x, 0.0, R) 
    ylist = Vy(x, 0.0, R) 
    
    plt.plot(x, dlist)
    plt.plot(x, ylist) 
    
    return 
#-------------------------------------------------------------------
    
#------Same as above but for x = -0.5 line:-------------------------
    
def dylinesx5(R): 
    
    y = np.linspace(-1.0, 1.0, 200)
    
    dlist = [] 
    ylist = [] 
    
    dlist = Vdelta(-0.5, y, R)
    ylist = Vy(-0.5, y, R)
    
    plt.plot(y, dlist)
    plt.plot(y, ylist)
    
    return 
#-------------------------------------------------------------------
#*******************************************************************
#___________________________________________________________________
    
    
#------------****RUN-TIME CODE****----------------------------------

#*********DEMONSTRATION OF LINE TAKEN:***********

#xlist, ylist, alist, blist, clist, rlist, vlist = x5line() 
#circleplot()

#---------Bin Splitting:------------
#xlist, ylist, alist, blist, clist, rlist, vlist = binsplit(xlist, 
                        #ylist, alist, blist, clist, rlist, vlist, 0.0, 4.0)

#-----Plot line to be evaluated:------- 
#plotscatter(xlist, ylist, 'r')

#plt.show() 

#********SCATTER VS LINES PLOTS: ****************

#----y = 0 line:----
ploty0(0.0, 15.0) 
dylinesy0(2.0)

#----x = -0.5 line:---- 
#plotx5(0.0,15.0) 
#dylinesx5(10.0) 


plt.show() 
#------------------------------------------------------------------
