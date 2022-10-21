import numpy as np 
import pyarma as pa 
import matplotlib.pyplot as plt 

data = pa.cube()

data.load("data/r_Euler.dat")
r_e = np.array(data)
data.load("data/v_Euler.dat")
v_e = np.array(data)

data.load("data/r_RK4.dat")
r_RK4 = np.array(data)
data.load("data/r_RK4.dat")
v_RK4 = np.array(data)

print(np.shape(r_e))
plt.plot(r_e[0,0,:], r_e[1,0,:], label="Euler")
plt.plot(r_RK4[0,0,:], r_RK4[1,0,:], label="RK4")
plt.legend()
plt.show()


