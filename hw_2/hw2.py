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
    step_arr=np.linspace(0,D*3.086*10.**18,1000) #Define array 0 to D stepsize ds
    int_arr=np.zeros((len(step_arr),1)) #Initialize another array using same steps but holding intensity values
    int_arr[0]=I_0 #set first element to I_0
    opt_dep=cro_sec*float(D)*float(n)*3.086*10.**18
    temp=0
    for i in range(0,len(int_arr)):
        temp=temp+S*math.e**(-opt_dep[i])
        if i != 0:
            int_arr[i]=I_0*math.e**(-opt_dep[i])+temp
    return(float(int_arr[len(int_arr)-1]))

def function3(v,cro_sec_0):
    gauss=cro_sec_0*math.e**(-np.square(v-500)/500)
    return gauss

###########################################################################################################################################################

print('Problem 1:')

#Enter cloud parameters
D=100 #Cloud depth in pc
n=1 #Density in cm^-3
opt_dep=np.array([10.**-3,1.,10.**3]) #Optical depth is 10^-3

col_den,cro_sec=function1(D,n,opt_dep[0])
cro_sec=np.zeros((3,1))
print('The column density of the cloud, given D=%3.0f pc and n=%1.0f cm^-3, is: N=%4.3e cm^-2.' % (D,n,col_den))

for i in range(0,3):

    col_den,cro_sec[i]=function1(D,n,opt_dep[i])

    print('The cross section needed to have a total optical depth of %1.0e is: %4.3e cm^2' % (opt_dep[i],cro_sec[i]))

################################################################################################################################################################

print('Problem 2:')

print('The function has been made!')

################################################################################################################################################################

print('Problem 3:')

print('The function has been made!')

v=np.linspace(1,1000,1000)
cro_secs=[[]]*3
for i in range(0,3):

    cro_secs[i]=function3(v,cro_sec[i])
    plt.plot(v,cro_secs[i])
    plt.xlabel('frequency (arbitrary units)')
    plt.ylabel('cross section (arbitrary units)')
    plt.title('Cross Section vs Frequency for Max = %4.3e cm^2' % cro_sec[i])
    plt.show()

#################################################################################################################################################################

print('Problem 4:')

print('(a)')
I_0=0
T=10**-3
I_D=np.zeros((1000,))
x=cro_secs[0]
for i in range(0,1000):
    I_D[i]=function2(D,n,x,I_0,1)

plt.plot(x,I_D)
plt.show()
