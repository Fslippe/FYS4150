import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa
A = pa.mat()
B = pa.mat()
A.load("data/x_u1000.dat")
print(A)
x = A[:,0]
u = A[:,1]

plt.plot(x,u)
plt.title('Plot of $u(x) = 1 - (1 - e^{-10})x - e^{-10x}$')
plt.xlabel("$x$")
plt.ylabel("$u(x)$")
#plt.savefig("../figures/x_u_plot.pdf")
plt.show()

B.load("data/n10.dat")
plt.plot(B[:,0], B[:,1], label="n=10")
B.load("data/n100.dat")
plt.plot(B[:,0], B[:,1], label="n=100")
B.load("data/n1000.dat")
plt.title("plot of analytical and nummerical solution for different number of steps")
plt.plot(B[:,0], B[:,1], label="n=1000")
plt.plot(x,u, label="Analytic")
plt.legend()
plt.savefig("../figures/analytical_compare.png")
plt.show()
