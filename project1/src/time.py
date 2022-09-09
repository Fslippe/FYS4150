import matplotlib.pyplot as plt
import pyarma as pa
import numpy as np

#Using armadillo to import data files
A = pa.mat()
B = pa.mat()
#Creating arrays to fill
mean_g = np.zeros(6)
std_g = np.zeros(6)
mean_s = np.zeros(6)
std_s = np.zeros(6)
x = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6] #list with number of steps for the different datafiles

#Loading arrays for different number of steps and assigning values to arrays
#General algorithm:
A.load("data/general_algorithm_time10.dat")
mean_g[0] = pa.mean(A)[0]
std_g[0] = pa.stddev(A)[0]
A.load("data/general_algorithm_time100.dat")
mean_g[1] = pa.mean(A)[0]
std_g[1] = pa.stddev(A)[0]
A.load("data/general_algorithm_time1000.dat")
mean_g[2] = pa.mean(A)[0]
std_g[2] = pa.stddev(A)[0]
A.load("data/general_algorithm_time10000.dat")
mean_g[3] = pa.mean(A)[0]
std_g[3] = pa.stddev(A)[0]
A.load("data/general_algorithm_time100000.dat")
mean_g[4] = pa.mean(A)[0]
std_g[4] = pa.stddev(A)[0]
A.load("data/general_algorithm_time1000000.dat")
mean_g[5] = pa.mean(A)[0]
std_g[5] = pa.stddev(A)[0]

#Special algorithm:
B.load("data/special_algorithm_time10.dat")
mean_s[0] = pa.mean(B)[0]
std_s[0] = pa.stddev(B)[0]
B.load("data/special_algorithm_time100.dat")
mean_s[1] = pa.mean(B)[0]
std_s[1] = pa.stddev(B)[0]
B.load("data/special_algorithm_time1000.dat")
mean_s[2] = pa.mean(B)[0]
std_s[2] = pa.stddev(B)[0]
B.load("data/special_algorithm_time10000.dat")
mean_s[3] = pa.mean(B)[0]
std_s[3] = pa.stddev(B)[0]
B.load("data/special_algorithm_time100000.dat")
mean_s[4] = pa.mean(B)[0]
std_s[4] = pa.stddev(B)[0]
B.load("data/special_algorithm_time1000000.dat")
mean_s[5] = pa.mean(B)[0]
std_s[5] = pa.stddev(B)[0]

#plots:
plt.title("Mean time over 1000 runs for different number of steps")
plt.plot(x, mean_g, label="General algorithm")
plt.plot(x, mean_s, label="Special algorithm")
plt.xlabel("Number of steps (n)")
plt.ylabel("Time (s)")
plt.xscale("log")
plt.legend()
#plt.savefig("../figures/mean_time.pdf")
plt.show()

plt.title("Special algorithm mean time compared to general algorithm")
plt.plot(x, mean_s/mean_g, label="Special algorithm")
plt.xlabel("Number of steps (n)")
plt.ylabel("$mean_s/mean_g$")
plt.xscale("log")
plt.legend()
plt.savefig("../figures/mean_time_rel.pdf")
plt.show()

plt.title("Standard deviation over 1000 runs for different number of steps")
plt.plot(x, std_g, label="General algorithm")
plt.plot(x, std_s, label="Special algorithm")
plt.xscale("log")
plt.xlabel("Number of steps (n)")
plt.ylabel("$\sigma$ (s)")
plt.legend()
#plt.savefig("../figures/std_time.pdf")
plt.show()
