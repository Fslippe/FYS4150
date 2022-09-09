import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

A = pa.mat()
mean_g = np.zeros(6)
std_g = np.zeros(6)
B = pa.mat()
mean_s = np.zeros(6)
std_s = np.zeros(6)

A.load("general_algorithm_time10.dat")
mean_g[0] = pa.mean(A)[0]
std_g[0] = pa.stddev(A)[0]
A.load("general_algorithm_time100.dat")
mean_g[1] = pa.mean(A)[0]
std_g[1] = pa.stddev(A)[0]
A.load("general_algorithm_time1000.dat")
mean_g[2] = pa.mean(A)[0]
std_g[2] = pa.stddev(A)[0]
A.load("general_algorithm_time10000.dat")
mean_g[3] = pa.mean(A)[0]
std_g[3] = pa.stddev(A)[0]
A.load("general_algorithm_time100000.dat")
mean_g[4] = pa.mean(A)[0]
std_g[4] = pa.stddev(A)[0]
A.load("general_algorithm_time1000000.dat")
mean_g[5] = pa.mean(A)[0]
std_g[5] = pa.stddev(A)[0]


B.load("special_algorithm_time10.dat")
mean_s[0] = pa.mean(B)[0]
std_s[0] = pa.stddev(B)[0]
B.load("special_algorithm_time100.dat")
mean_s[1] = pa.mean(B)[0]
std_s[1] = pa.stddev(B)[0]
B.load("special_algorithm_time1000.dat")
mean_s[2] = pa.mean(B)[0]
std_s[2] = pa.stddev(B)[0]
B.load("special_algorithm_time10000.dat")
mean_s[3] = pa.mean(B)[0]
std_s[3] = pa.stddev(B)[0]
B.load("special_algorithm_time100000.dat")
mean_s[4] = pa.mean(B)[0]
std_s[4] = pa.stddev(B)[0]
B.load("special_algorithm_time1000000.dat")
mean_s[5] = pa.mean(B)[0]
std_s[5] = pa.stddev(B)[0]

x = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6]
plt.plot(x, mean_g)
plt.plot(x, mean_s)
plt.xscale("log")

plt.show()
plt.plot(x, std_g)
plt.plot(x, std_s)
plt.xscale("log")

plt.show()
