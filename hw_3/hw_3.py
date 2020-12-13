import matplotlib.pyplot as plt
import math
from astropy.io import fits
import numpy as np

Ze=1.6*10**-19 #fixed ion charge at the origin, coulombs
m_e=9.109*10**-31 #electron mass in kg
a_0=5.29*10**-11 #Bohr radius, meters

e_e=-1.6*10**-19 #electron charge in coulombs
x_e=-400*a_0 #electron initial position in x
y_e=400*a_0 #electron initial position in y
v_e_x=10**5 #electron x velocity, m/s

k=np.float128(8.99*10**9) #Coulomb's constant in SI units

r_x=np.zeros(100000,dtype=np.float128) #initialize position array
r_y=np.zeros(100000,dtype=np.float128) #initialize position array
a_x=np.zeros(100000,dtype=np.float128)
a_y=np.zeros(100000,dtype=np.float128)
v_x=np.zeros(100000,dtype=np.float128)
v_y=np.zeros(100000,dtype=np.float128)

r_x[0]=np.float128(x_e)
r_y[0]=np.float128(y_e)

v_x[0]=np.float128(v_e_x)
v_y[0]=np.float128(0)
d_t=np.float128(10**-17) #time interival

for i in range(1,100000):

    r=np.float128(np.sqrt(r_x[i-1]**2+r_y[i-1]**2))
    force=k*Ze*e_e/(r**2)
    force_x=force*np.cos(r_x[i-1]/r)*np.sign(r_x[i-1])
    force_y=force*np.sin(r_y[i-1]/r)*np.sign(r_y[i-1])

    a_x[i-1]=force_x/m_e
    a_y[i-1]=force_y/m_e

    v_x[i]=v_x[i-1]+a_x[i-1]*d_t
    v_y[i]=v_y[i-1]+a_y[i-1]*d_t

    r_x[i]=r_x[i-1]+v_x[i-1]*d_t
    r_y[i]=r_y[i-1]+v_y[i-1]*d_t
   

time=np.linspace(0,10000,100000)*d_t

r=np.float128(np.sqrt(r_x[99999]**2+r_y[99999]**2))
force=k*Ze*e_e/(r**2)
force_x=force*np.cos(r_x[99999]/r)*np.sign(r_x[99999])
force_y=force*np.sin(r_y[99999]/r)*np.sign(r_y[99999])
a_x[99999]=force_x/m_e
a_y[99999]=force_y/m_e


plt.plot(r_x,r_y)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Particle Path')
plt.ticklabel_format(axis="both",style="sci",scilimits=(0,0))
plt.scatter([x_e,0],[y_e,0],marker=(5,1),c='y')
plt.show()

plt.plot(time,v_x)
plt.xlabel('time')
plt.ylabel('X velocity')
plt.title('X Velocity vs. Time')
plt.ticklabel_format(axis="both",style="sci",scilimits=(0,0))
plt.show()

plt.plot(time,v_y)
plt.xlabel('Time')
plt.ylabel('Y Velocity')
plt.title('Y Velocity vs. Time')
plt.ticklabel_format(axis="both",style="sci",scilimits=(0,0))
plt.show()

plt.plot(time,a_x)
plt.xlabel('Time')
plt.ylabel('X Acceleration')
plt.title('X Acceleration vs. Time')
plt.ticklabel_format(axis="both",style="sci",scilimits=(0,0))
plt.show()

plt.plot(time,a_y)
plt.xlabel('Time')
plt.ylabel('Y Acceleration')
plt.title('Y Acceleration vs. Time')
plt.ticklabel_format(axis="both",style="sci",scilimits=(0,0))
plt.show()

print(v_x[0])
