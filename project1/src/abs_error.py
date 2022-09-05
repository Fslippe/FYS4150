import numpy as np
import matplotlib.pyplot as plt

x_u = np.loadtxt("x_u.txt")
x_v = np.loadtxt("x_v.txt")

plt.plot(np.log10(abs(x_u[:,1] - x_v[:,1])))
plt.title("Absolute error")
plt.show()


rewwkjw
