import matplotlib.pyplot as plt
import math
from astropy.io import fits
import numpy as np

def function1(D,n,opt_dep):
    col_den=float(D)*float(n) #Calculation for column density, obtained by integrating volumetric density over a column, in this case, depth D times density n
    col_den=col_den*3.086*10.**18#Parsec to centimeter conversion to put column density in units of cm^-2
    cro_sec=opt_dep/col_den # Cross section is equal to optical depth divided by column density.
    return(col_den,cro_sec)

def function2(D,n,cro_sec,I_0,S):
    step_arr=np.linspace(0,D*3.086*10.**18,10000) #Define array 0 to D stepsize ds
    int_arr=np.zeros((1,len(step_arr))) #Initialize another array using same steps but holding intensity values
    int_arr[0]=I_0 #set first element to I_0
    col_den=float(D)*float(n)*3.086*10.**18
    opt_dep=cro_sec*float(D)*float(n)*3.086*10.**18 #Optical depth is equal to cross section times column density
    Si=np.ones((1,len(step_arr)))
    for i in range(1,len(int_arr)):
        int_arr[i]=int_arr[0]*math.e**(-opt_dep)+Si[i]
   
    print(len(int_arr))
    return(int_arr)

###########################################################################################################################################################

print('Problem 1:')

#Enter cloud parameters
D=100 #Cloud depth in pc
n=1 #Density in cm^-3
opt_dep=10**-3 #Optical depth is 10^-3

col_den,cro_sec=function1(D,n,opt_dep)
print('The column density of the cloud, given D=%3.0f pc and n=%1.0f cm^-3, is: N=%4.3e cm^-2.' % (D,n,col_den))
print('The cross section needed to have a total optical depth of %1.0e is: %4.3e cm^2' % (opt_dep,cro_sec))

opt_dep=1
col_den,cro_sec=function1(D,n,opt_dep)

print('The cross section needed to have a total optical depth of %1.0e is: %4.3e cm^2' % (opt_dep,cro_sec))


opt_dep=10**3
col_den,cro_sec=function1(D,n,opt_dep)

print('The cross section needed to have a total optical depth of %1.0e is: %4.3e cm^2' %(opt_dep,cro_sec))

################################################################################################################################################################

print('Prolem 2:')

S=np.ones((1,100000))
I_D=function2(D,n,cro_sec,0,np.ones((1,10000)))
print(I_D)
