import matplotlib.pyplot as plt
import math
from astropy.io import fits
import numpy as np

def function1(D,n,opt_dep):
    #Function for problem 1, takes a depth D, density n, and optical depth opt_dep.  Returns a 2 element list containing the column density and cross section 

    col_den=float(D)*float(n) #Calculation for column density, obtained by integrating volumetric density over a column, in this case, depth D times density n
    col_den=col_den*3.086*10.**18 #Parsec to centimeter conversion to put column density in units of cm^-2
    cro_sec=opt_dep/col_den #Cross section is equal to optical depth divided by column density.
    
    return(col_den,cro_sec)

def function2(D,n,cro_sec,I_0,S):
    #Fucntion for problem 2, cro_sec is a single value

    step_arr=np.linspace(0,D*3.086*10.**18,1000) #Define array 0 to D with stepsize ds
    int_arr=np.zeros((len(step_arr),)) #Initialize another array using same steps but holding intensity values
    int_arr[0]=I_0 #set first element to I_0
    opt_dep=cro_sec*float(D)*float(n)*3.086*10.**18
    opt_deps=np.linspace(0,opt_dep,1000)
    temp=0.
    
    for i in range(0,len(int_arr)):
    
        temp=temp+S*(math.e**(-opt_dep+opt_deps[i]))*opt_dep/1000
        if i != 0:
            int_arr[i]=I_0*(math.e**(-opt_deps[i]))+temp

    return(float(int_arr[len(int_arr)-1]))

def function3(v,cro_sec_0,center,width):
    #Function for problem 3.  returns an array of cross sections in the shape of a gaussian given an input array of frequencies.

    gauss=cro_sec_0*math.e**(-np.square(v-center)/width**2)
    print(len(gauss))
    return gauss

def function4(v,cro_secs,I_0,S):
    #Function for problem 4. Plots intensity versus frequency given conditions
    #Cro_secs is an array of cross sections
    #I_0 is initial intensity
    #S is source function
    #T is optical depth

    D=100.
    n=1.
    I_D=np.zeros((1000,))
    
    for i in range(0,len(cro_secs)):
    
        I_D[i]=function2(D,n,cro_secs[i],I_0,S)

    plt.plot(v,I_D)
    plt.xlabel('Frequency')
    plt.ylabel('Intensity')
    plt.title('Observed Intensity from cloud')
    plt.show()

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

v=np.linspace(10**-7,10**-6,1000)
cro_secs=[[]]*4
center=5*10**-7
width=2.5*10**-8
for i in range(0,3):
    
    cro_secs[i]=function3(v,cro_sec[i],center,width)
    plt.plot(v,cro_secs[i])
    plt.xlabel('frequency')
    plt.ylabel('cross section (cm^-2)')
    plt.title('Cross Section vs Frequency for Max = %4.3e cm^2' % cro_sec[i])
    plt.show()

#################################################################################################################################################################

print('Problem 4:')

print('(a)')

I_0=0
S=1

function4(v,cro_secs[0],I_0,S) #cro_secs[0] gives cross sections for optical depth of 10^-3

print('(b)')

I_0=10
S=1

function4(v,cro_secs[0],I_0,S)

print('(c)')

I_0=1
S=10

function4(v,cro_secs[0],I_0,S)

print('(d)')

I_0=1
S=1
cro_secs[3]=np.ones((1000,))*100

function4(v,cro_secs[3],I_0,S)

print('(e)')

I_0=1
S=1000

function4(v,cro_secs[2],I_0,S)

print('(f)')

I_0=1000
S=1

function4(v,cro_secs[2],I_0,S)

