import csv, math

#Extracts data from csv files from old Koma data and outputs 
#a data text file in the format r2,a,b,c,x,y,z,v to be used in the 
#data-analysis apparatus built for the new Koma data.



#------Define symmeteries from 6 possible configurations of 3 bodies:------
def perm(i,a,b,c):
  perm={0:(a,b,c),1:(a,c,b),2:(b,c,a),3:(b,a,c),4:(c,a,b),5:(c,b,a)}
  return perm[i]
#--------------------------------------------------------------------------


   
#--------------Files:-------------------------------------------------------
#------select csv data----------:
#c=csv.reader(open('beta57.csv'),delimiter="&") 
#c=csv.reader(open('beta58.csv'),delimiter="&") 
c=csv.reader(open('beta60.csv'),delimiter="&") 


#-------select datafile to write to:---------------------
#f = open('data57', 'w')
#f = open('data58', 'w')
f = open('data60', 'w')

#***NOTE: CHECK THE OPTIONS ABOVE ARE FOR THE SAME BETA VALUE***
#---------------------------------------------------------------------------

#-----Note on (i, j, k) notation: ------------
#for the three quark positions r1, r2, r3: 
#In this data: 

#r1 = (i, 0, 0), r2 = (0, j, 0), r3 = (0, 0, k)
#-----------------------------------------------


#----------Read csv, produce data, write to output file:--------------------
for l in c:
  #split xyz stirng into an array: 
  ijk = l[0].replace("(", "").replace(")", "").replace(" ", "").split(",")
  i = float(ijk[0])
  j = float(ijk[1])
  k = float(ijk[2]) 
         
  #strip out potential value 
  #(split trick to remove bracket value at the end of string): 
  v = l[1].split("(")[0] 

  #calculate a, b, c data from the (x,y,z): 
  a1 = math.sqrt(i*i + j*j)
  b1 = math.sqrt(j*j + k*k)
  c1 = math.sqrt(i*i + k*k) 
  
  #For all six possible permutations: 
  for i in range(6): 
    a,b,c = perm(i, a1, b1, c1) 
    
    #------------#Calculate hyper-radius and (x, y, z):---------------
    x=(c*c+a*a-b*b)/a/2
    y=math.sqrt(c*c-x*x)
  
    #rho and lambda coordinates:
    rx=a/math.sqrt(2)
    ry=0
    lx=(a-2*x)/math.sqrt(6)
    ly=(-2*y)/math.sqrt(6)
  
    #calculate x, y, z and hyper-radius^2 values: 
    r2=rx*rx+ry*ry+lx*lx+ly*ly
    xx=2*(rx*lx+ry*ly)/r2
    yy=(rx*rx+ry*ry-lx*lx-ly*ly)/r2
    zz=2*(rx*ly-ry*lx)/r2 
    #-----------------------------------------------------------------
    
    #-------Write r2, a, b, c, x, y, z, v to data file:---------------
    f.write(str(r2) + ' ' + str(a) + ' ' + str(b) + ' ' + str(c) + ' '
    + str(xx) + ' ' + str(yy) + ' ' + str(zz) + ' ' + v + '\n')    
    
#---------------------------------------------------------------------------
  
f.close() 

