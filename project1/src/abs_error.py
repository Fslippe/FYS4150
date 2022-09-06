import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa
plt.rcParams.update({'font.size': 22})
A = pa.mat()

A.load("n10.dat")
np.loadtxt("x_u.txt")
plt.plot(A[:,0], A[:,1], label="Steps: 10")
A.load("n99.dat")
plt.plot(A[:,0], A[:,1], label="Steps: 100")
A.load("n1000.dat")
plt.plot(A[:,0], A[:,1], label="Steps: 1000")
A.load("n10000.dat")
plt.plot(A[:,0], A[:,1], label="Steps: 10000")
plt.legend()
plt.show()
