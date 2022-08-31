import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("x_u.txt")
x = data [:,0]
u = data [:,1]

plt.plot(x,u)
plt.xlabel("$x$")
plt.ylabel("$u(x)$")
plt.savefig("/Users/alessiocanclini/FYS4150/project1/figures/x_u_plot.pdf")
