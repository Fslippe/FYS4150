import seaborn as sns
import numpy as np 
import pyarma as pa 
import matplotlib.pyplot as plt 
sns.set_style("darkgrid")

def analytic(T):
    J = 1
    beta = 1/T
    Z = 2 * np.exp(beta*8) + 2*np.exp(-beta*8) + 12
    E = (J / Z) * (16*np.exp(- beta * 8 * J) - 16*np.exp(beta * 8 * J))
    E_2 = (J*J / Z) * (128*np.exp(- beta * 8 * J) + 128*np.exp(beta * 8 * J))
    e = E / 4   
    M = (8*np.exp(beta * 8 * J) + 16) / Z
    M_2 = (32*np.exp(beta * 8 * J) + 32) / Z
    m = M / 4
    C_v = (E_2 - E*E) / (4 * T*T)
    X  = (M_2 - M*M) / (4 * T)
    return E, E_2, e, M, M_2, m, C_v, X

def plot_T(data):

    t = data[0,:]
    E, E_2, e, M, M_2, m, C_v, X = analytic(t) 

    plt.subplot(2, 2, 1)
    sns.lineplot(t, e, color="r")
    plt.xlabel(r"$T [J/k_b]$")
    plt.ylabel(r"$\left<\epsilon\right> [J]$")
    sns.scatterplot(t, data[1,:])

    plt.subplot(2, 2, 2)
    sns.lineplot(t, m, color="r")
    sns.scatterplot(t, data[2,:])
    plt.xlabel(r"$T [J/k_b]$")
    plt.ylabel(r"$\left< m \right> [J]$")

    plt.subplot(2, 2, 3)
    sns.lineplot(t, C_v, color="r")
    sns.scatterplot(t, data[3,:])
    plt.xlabel(r"$T [J/k_b]$")
    plt.ylabel(r"$C_v [k_b]$")
 
    plt.subplot(2, 2, 4)
    sns.lineplot(t, X, color="r")
    sns.scatterplot(t, data[4,:], label="Numeric")
    plt.xlabel(r"$T [J/k_b]$")
    plt.ylabel(r"$\chi[k_b^{-1}]$")

    plt.show()

def plot_diff(data):
    t = 1
    E, E_2, e, M, M_2, m, C_v, X = analytic(t) 
    cycles = data[0,:]
    print(cycles)
    e_diff = np.abs(data[1,:] - e)
    print(e_diff)
    m_diff = np.abs(data[2,:] - m)
    C_v_diff = np.abs(data[3,:] - C_v)
    X_diff = np.abs(data[4,:] - X)

    plt.subplot(2, 2, 1)

    sns.lineplot(cycles, e_diff, markers=True, label="Analytic")
    plt.xscale("log")
    plt.yscale("log")

    plt.subplot(2, 2, 2)

    sns.lineplot(cycles, m_diff, markers=True, label="Analytic")
    plt.xscale("log")
    plt.yscale("log")

    plt.subplot(2, 2, 3)

    sns.lineplot(cycles, C_v_diff, markers=True, label="Analytic")
    plt.xscale("log")
    plt.yscale("log")

    plt.subplot(2, 2, 4)
    sns.lineplot(cycles, X_diff, markers=True, label="Analytic")
    plt.xscale("log")
    plt.yscale("log")
    plt.show()


def main():
    data_1 = pa.mat()
    data_2 = pa.mat()
    data_2.load("data/T_val_1mill.dat")
    data_1.load("data/cycle_val.dat")
    plot_diff(np.array(data_1))

    plot_T(np.array(data_2))

if __name__ == "__main__":
    main()
