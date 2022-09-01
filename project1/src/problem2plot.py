import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("x_u.txt")
x = data [:,0]
u = data [:,1]

plt.plot(x,u)
plt.xlabel("$x$")
plt.ylabel("$u(x)$")
plt.savefig("../figures/x_u_plot.pdf")
