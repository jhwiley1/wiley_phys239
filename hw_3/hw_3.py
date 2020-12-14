import matplotlib.pyplot as plt
import math
from astropy.io import fits
import numpy as np
from scipy.fft import fft

def function1(y_0,v_0,plott,plottt):

    Ze=1.6*10**-19 #fixed ion charge at the origin, coulombs
    m_e=9.109*10**-31 #electron mass in kg
    a_0=5.29*10**-11 #Bohr radius, meters

    e_e=-1.6*10**-19 #electron charge in coulombs
    x_0=np.float128(-300)
    x_0=x_0*a_0 #electron initial position in x
    y_0=y_0*a_0 #electron initial position in y
    v_e_x=v_0*10**5

    N = 100000

    k=np.float128(8.99*10**9) #Coulomb's constant in SI units

    r_x=np.zeros(N,dtype=np.float128) #initialize position array
    r_y=np.zeros(N,dtype=np.float128) #initialize position array
    a_x=np.zeros(N,dtype=np.float128) #initialize x acceleration array
    a_y=np.zeros(N,dtype=np.float128) #initialize y acceleration array
    v_x=np.zeros(N,dtype=np.float128) #initialize x velocity array
    v_y=np.zeros(N,dtype=np.float128) #initialize y velocity array

    r_x[0]=np.float128(x_0) #set initial position of x position array to x_0
    r_y[0]=np.float128(y_0) #set initial position of y position array to y_0

    v_x[0]=np.float128(v_e_x)
    v_y[0]=np.float128(0)
    d_t=np.float128(10**-17) #time interival
    
    for i in range(1,N):

        r=np.float128(np.sqrt(r_x[i-1]**2+r_y[i-1]**2))
        force=k*Ze*e_e/(r**2)
        force_x=force*np.cos(r_x[i-1]/r)*np.sign(r_x[i-1])
        force_y=force*np.abs(np.sin(r_y[i-1]/r))*np.sign(r_y[i-1])

        a_x[i-1]=force_x/m_e
        a_y[i-1]=force_y/m_e

        v_x[i]=v_x[i-1]+a_x[i-1]*d_t
        v_y[i]=v_y[i-1]+a_y[i-1]*d_t

        r_x[i]=r_x[i-1]+v_x[i-1]*d_t
        r_y[i]=r_y[i-1]+v_y[i-1]*d_t
   

    time=np.linspace(0,N,N)*d_t

    r=np.float128(np.sqrt(r_x[N-1]**2+r_y[N-1]**2))
    force=k*Ze*e_e/(r**2)
    force_x=force*np.cos(r_x[N-1]/r)*np.sign(r_x[N-1])
    force_y=force*np.sin(r_y[N-1]/r)*np.sign(r_y[N-1])
    a_x[N-1]=force_x/m_e
    a_y[N-1]=force_y/m_e

    if plott == 1:

        plt.plot(r_x[0:10000],r_y[0:10000]) #Plot particle path, x and y
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Particle Path')
        plt.ticklabel_format(axis="both",style="sci",scilimits=(0,0))
        plt.scatter([x_0,0],[y_0,0],marker=(5,1),c=['r','k'])
        plt.show()

        plt.plot(time[0:10000],v_x[0:10000]) #plot x velocity
        plt.xlabel('time')
        plt.ylabel('X velocity')
        plt.title('X Velocity vs. Time')
        plt.ticklabel_format(axis="both",style="sci",scilimits=(0,0))
        plt.show()

        plt.plot(time[0:10000],v_y[0:10000]) #plot y velocity
        plt.xlabel('Time')
        plt.ylabel('Y Velocity')
        plt.title('Y Velocity vs. Time')
        plt.ticklabel_format(axis="both",style="sci",scilimits=(0,0))
        plt.show()

        plt.plot(time[0:10000],a_x[0:10000]) #plot X acceleration
        plt.xlabel('Time')
        plt.ylabel('X Acceleration')
        plt.title('X Acceleration vs. Time')
        plt.ticklabel_format(axis="both",style="sci",scilimits=(0,0))
        plt.show()

        plt.plot(time[0:10000],a_y[0:10000]) #plot Y acceleration
        plt.xlabel('Time')
        plt.ylabel('Y Acceleration')
        plt.title('Y Acceleration vs. Time')
        plt.ticklabel_format(axis="both",style="sci",scilimits=(0,0))
        plt.show()

    accel=np.sqrt(a_x**2+a_y**2)

    yf=fft(accel)
    xf=np.linspace(0.0,1.0/(2.0*d_t),N)
    y_plot=2.0/N * np.abs(yf[0:N])

    temp=0
    i=2

    while temp == 0:

        if y_plot[i]>y_plot[i-1] and y_plot[i-1]>y_plot[i-2] and y_plot[i-2]>y_plot[i-3] and y_plot[i-3]>y_plot[i-4] and y_plot[i-4]>y_plot[i-5]:
            temp=1
            peak_y=np.max(y_plot[i-5:i+500])
            peak_pos=np.argmax(y_plot[i-5:i+500])
            peak_x=xf[peak_pos]
            
            if plottt == 1:

                plt.plot(xf[i-5:i+500],y_plot[i-5:i+500])
                plt.xlabel('v')
                plt.ylabel('S')
                plt.title('Power Spectrum of Emitted Radiation for v=%2.1e and b=%2.1e'% (v_e_x, y_0))
                plt.show()
                
        i=i+1
    return peak_x

y_0=input('Enter Y_0 in units of bohr radii (ex. 200): ')
v_0=input('Enter V_0 in units of 10^5 m/s (ex. 10): ')

y_0=np.float128(y_0)
v_0=np.float128(v_0)

function1(y_0,v_0,1,1)
y_0_array=[130,135,140,145,150,155,160,165,170,175,180,185,190,195,200]
v_0=10
peak_array=np.zeros(15)

print('Generating peak frequency vs. impact parameter plot...')

for i in range(0,15):
    
    y_0=np.float128(y_0_array[i])
    v_0=np.float128(v_0)
    peak_array[i]=function1(y_0,v_0,0,0)

plt.plot(y_0_array,peak_array)
plt.xlabel('b (a0)')
plt.ylabel('v')
plt.title('Peak Frequency vs Impact Parameter for V_0=10^6 m/s')
plt.show()


v_0_array=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]
y_0=200

peak_array=np.zeros(20)

print('Generating peak frequency vs. initial velocity plot...')

for i in range(0,20):

    y_0=np.float128(y_0)
    v_0=np.float128(v_0_array[i])
    peak_array[i]=function1(y_0,v_0,0,0)

plt.plot(v_0_array,peak_array)
plt.xlabel('velocity')
plt.ylabel('peak frequency')
plt.title('Peak Frequency vs. Initial Velocity for b = 250 a0')
plt.show()
