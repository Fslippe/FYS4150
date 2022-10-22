import numpy as np 
import pyarma as pa 
import matplotlib.pyplot as plt 
import os
plt.rcParams.update({"lines.linewidth": 2})
plt.rcParams.update({"font.size": 14})


def run(N, T, n, interaction, omega_min, omega_max, omega_step, compile=False):
    """

    """
    if compile == True:
        print("Compiling...\n")
        compile_cpp = os.system("g++ -std=c++11 time_dep.cpp penningtrap.cpp particle.cpp -o time_dep -larmadillo -O2"),
    run_cpp = os.system("./time_dep %s %s %s %s %s %s %s" %(N, T, n, interaction, omega_min, omega_max, omega_step))   

  
    data = pa.mat()
    data.load("data/frac_p_left.dat")
    r = np.array(data)
    return r

def plot_p_fraction_frequency():
    r = run(10000, 500, 10, "false", 0.2, 2.5, 0.2, False)
    plt.plot(r[:,0], r[:,1], label = "$f = 0.1$")
    plt.plot(r[:,0], r[:,2], label = "$f = 0.4$")
    plt.plot(r[:,0], r[:,3], label = "$f = 0.7$")
    plt.legend()
    plt.show()


plot_p_fraction_frequency()
