import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa
#plt.rcParams.update({'font.size': 22})
A = pa.mat()
B = pa.mat()


def abs_error(file1, file2, steps):
    A.load(file1)
    B.load(file2)
    plt.plot(A[:,0], np.log10(pa.abs(B[:,1] - A[:,1])), label="Steps: %i" %(steps))
    plt.xlabel("$x$")
    plt.ylabel("$\log_{10}(\Delta_i)$")

def rel_error(file1, file2, steps):
    A.load(file1)
    B.load(file2)
    plt.plot(A[:,0], np.log10(pa.abs((B[:,1] - A[:,1]) / B[:,1])), label="Steps: %i" %(steps))
    plt.xlabel("$x$")
    plt.ylabel("$\log_{10}(\epsilon_i)$")

def max_error(file1, file2):
    A.load(file1)
    B.load(file2)
    max = pa.max(pa.abs((B[:,1] - A[:,1]) / B[:,1]))
    print(max[0])
    return max

abs_error("data/n10.dat", "data/x_u10.dat", 10)
abs_error("data/n100.dat", "data/x_u100.dat", 100)
abs_error("data/n1000.dat", "data/x_u1000.dat", 1000)
abs_error("data/n10000.dat", "data/x_u10000.dat", 10000)
plt.title("Logarithm of the absolute error")
plt.legend()
plt.savefig("../figures/abs_error.pdf")
plt.show()

rel_error("data/n10.dat", "data/x_u10.dat", 10)
rel_error("data/n100.dat", "data/x_u100.dat", 100)
rel_error("data/n1000.dat", "data/x_u1000.dat", 1000)
rel_error("data/n10000.dat", "data/x_u10000.dat", 10000)
plt.title("Logarithm of the relative error")
plt.legend()
plt.savefig("../figures/rel_error.pdf")

plt.show()

e_10 = max_error("data/n10.dat", "data/x_u10.dat")
e_100 = max_error("data/n100.dat", "data/x_u100.dat")
e_1000 = max_error("data/n1000.dat", "data/x_u1000.dat")
e_10000 = max_error("data/n10000.dat", "data/x_u10000.dat")
e_100000 = max_error("data/n100000.dat", "data/x_u100000.dat")
e_1000000 = max_error("data/n1000000.dat", "data/x_u1000000.dat")
e_10000000 = max_error("data/n10000000.dat", "data/x_u10000000.dat")

plt.plot([1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7], [(e_10[0,:]), (e_100[0,:]), (e_1000[0,:]), (e_10000[0,:]), e_100000[0,:], e_1000000[0,:], e_10000000[0,:]], marker="o", linestyle="--")
plt.ylabel("$max(\epsilon$)")
plt.title("$max(\epsilon_i)$ for differnt n on a logarithmic scale")
plt.xlabel("n")
plt.xscale("log")
plt.yscale("log")
plt.savefig("../figures/max_rel_error.pdf")

plt.show()
