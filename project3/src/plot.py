import numpy as np 
import pyarma as pa 
import matplotlib.pyplot as plt 

data = pa.cube()
matrix = pa.mat()
data.load("data/r_Euler.dat")
r_e = np.array(data)
matrix.load("data/r_a.dat")
r_a = np.array(data)
data.load("data/v_Euler.dat")
v_e = np.array(data)
data.load("data/r_RK4.dat")
r_RK4 = np.array(data)
data.load("data/r_RK4.dat")
v_RK4 = np.array(data)

particles = np.size(r_e, axis=1)
print(r_e[1, 1, 0])
for i in range(particles):
    #plt.plot(r_e[1,i,:], label="Euler")

    plt.plot(r_e[0,i,:],r_e[1,i,:], label="RK4")
    #plt.plot(r_a[1,0,:], "--", label="analytic")

    plt.legend()
plt.show()

"""
for i in range(particles):
    plt.plot(r_e[0,i,:], r_e[1,i,:], label="Euler")
    plt.plot(r_RK4[0,i,:], r_RK4[1,i,:], label="RK4")
    plt.legend()
plt.show()
"""


