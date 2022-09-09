import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa
plt.rcParams.update({'font.size': 22})
A = pa.mat()
B = pa.mat()

def abs_error(file1, file2, steps):
    A.load(file1)
    B.load(file2)
    plt.plot(A[:,0], np.log10(pa.abs(B[:,1] - A[:,1])), label="Steps: %i" %(steps))

def rel_error(file1, file2, steps):
    A.load(file1)
    B.load(file2)
    plt.plot(A[:,0], np.log10(pa.abs((B[:,1] - A[:,1]) / B[:,1])), label="Steps: %i" %(steps))

def max_error(file1, file2):
    A.load(file1)
    B.load(file2)
    return  pa.max(pa.abs((B[:,1] - A[:,1]) / B[:,1]))

a = -1/2

A.load("special_n100.dat")
plt.plot(A[:,0], A[:,1])
plt.show()

for i in range(20):
    a = -1/(2 - 1*a)
    print(a)
abs_error("n10.dat", "x_u10.dat", 10)
abs_error("n100.dat", "x_u100.dat", 100)
abs_error("n1000.dat", "x_u1000.dat", 1000)
abs_error("n10000.dat", "x_u10000.dat", 10000)


plt.title("Logarithm of the absolute error")
plt.legend()
plt.show()

rel_error("n10.dat", "x_u10.dat", 10)
rel_error("n100.dat", "x_u100.dat", 100)
rel_error("n1000.dat", "x_u1000.dat", 1000)
rel_error("n10000.dat", "x_u10000.dat", 10000)
plt.title("Logarithm of the relative error")
plt.legend()
plt.show()

e_10 = max_error("n10.dat", "x_u10.dat")
e_100 = max_error("n100.dat", "x_u100.dat")
e_1000 = max_error("n1000.dat", "x_u1000.dat")
e_10000 = max_error("n10000.dat", "x_u10000.dat")
plt.plot([10, 100, 1000, 10000], [np.log(e_10[0,:]), np.log(e_100[0,:]), np.log(e_1000[0,:]), np.log(e_10000[0,:])])
plt.ylabel("log($max(\epsilon$))")
plt.xlabel("n")
plt.show()
