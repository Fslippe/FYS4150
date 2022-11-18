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

def plot_T(data, plot_analytic=True):
    t = data[0,:]
    E, E_2, e, M, M_2, m, C_v, X = analytic(t) 

    plt.figure()
    if plot_analytic:
        sns.lineplot(t, e, color="r", label="Analytic")
    plt.xlabel(r"$T$ $[J/k_b]$")
    plt.ylabel(r"$\left<\epsilon\right>$ $[J]$")
    sns.scatterplot(t, data[1,:],  label="Numeric")
    plt.savefig("../figures/numeric_analytic_e_T.pdf", dpi=300, bbox_inches='tight')

    plt.figure()
    if plot_analytic:
        sns.lineplot(t, m, color="r", label="Analytic")
    sns.scatterplot(t, data[2,:],  label="Numeric")
    plt.xlabel(r"$T$ $[J/k_b]$")
    plt.ylabel(r"$\left< m \right>$ $[J]$")
    plt.savefig("../figures/numeric_analytic_m_T.pdf", dpi=300, bbox_inches='tight')

    plt.figure()
    if plot_analytic:
        sns.lineplot(t, C_v, color="r", label="Analytic")
    sns.scatterplot(t, data[3,:],  label="Numeric")
    plt.xlabel(r"$T$ $[J/k_b]$")
    plt.ylabel(r"$C_v$ $[k_b]$")
    plt.savefig("../figures/numeric_analytic_c_v_T.pdf", dpi=300, bbox_inches='tight')

    plt.figure()
    if plot_analytic:
        sns.lineplot(t, X, color="r", label="Analytic")
    sns.scatterplot(t, data[4,:], label="Numeric")
    plt.xlabel(r"$T$ $[J/k_b]$")
    plt.ylabel(r"$\chi$ $[k_b^{-1}]$")
    plt.savefig("../figures/numeric_analytic_X_T.pdf", dpi=300, bbox_inches='tight')

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

    plt.figure()
    sns.lineplot(cycles, e_diff, markers=True, label=r"$|\epsilon_{Analytic}-\epsilon_{Numeric}|$")
    plt.xlabel(r"$MCMC$ $cycles$")
    plt.ylabel(r"$\Delta \left<\epsilon\right>$ $[J]$")
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("../figures/numeric_diff_e.pdf", dpi=300, bbox_inches='tight')

    plt.figure()
    sns.lineplot(cycles, m_diff, markers=True, label=r"$|m_{Analytic}-m_{Numeric}|$")
    plt.xlabel(r"$MCMC$ $cycles$")
    plt.ylabel(r"$\Delta \left< m \right>$ $[J]$")
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("../figures/numeric_diff_m.pdf", dpi=300, bbox_inches='tight')

    plt.figure()
    sns.lineplot(cycles, C_v_diff, markers=True, label=r"$|C_{v,Analytic}-C_{v,Numeric}|$")
    plt.xlabel(r"$MCMC$ $cycles$")
    plt.ylabel(r"$\Delta C_v$ $[k_b]$")
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("../figures/numeric_diff_c_v.pdf", dpi=300, bbox_inches='tight')

    plt.figure()
    sns.lineplot(cycles, X_diff, markers=True, label=r"$|\chi_{Analytic}-\chi_{Numeric}|$")
    plt.xlabel(r"$MCMC$ $cycles$")
    plt.ylabel(r"$\Delta \chi$ $[k_b^{-1}]$")
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("../figures/numeric_diff_X.pdf", dpi=300, bbox_inches='tight')
    plt.show()

def plot_data(data, T, order, savename):
    cycles = data[0,:]
    e = data[1,:]
    m = data[2,:]
    C_v = data[3,:]
    X = data[4,:]


    plt.figure()
    plt.title(r"%s  $T=$%.1f $J/k_b$" %(order, T))
    sns.lineplot(cycles, e, markers=True, label=r"$\epsilon$")
    plt.xlabel(r"$MCMC$ $cycles$")
    plt.ylabel(r"$\left<\epsilon\right>$ $[J]$")
    plt.xscale("log")
    plt.savefig("../figures/%s_e.pdf" %(savename), dpi=300, bbox_inches='tight')

    plt.figure()
    plt.title(r"%s  $T=$%.1f $J/k_b$" %(order, T))
    sns.lineplot(cycles, m, markers=True, label=r"$m$")
    plt.xlabel(r"$MCMC$ $cycles$")
    plt.ylabel(r"$\left< m \right>$ $[J]$")
    plt.xscale("log")
    plt.savefig("../figures/%s_m.pdf" %(savename), dpi=300, bbox_inches='tight')

    plt.figure()
    plt.title(r"%s  $T=$%.1f $J/k_b$" %(order, T))
    sns.lineplot(cycles, C_v, markers=True, label=r"$C_{v}$")
    plt.xlabel(r"$MCMC$ $cycles$")
    plt.ylabel(r"$C_v$ $[k_b]$")
    plt.xscale("log")
    plt.savefig("../figures/%s_c_v.pdf" %(savename), dpi=300, bbox_inches='tight')

    plt.figure()
    plt.title(r"%s  $T=$%.1f $J/k_b$" %(order, T))
    sns.lineplot(cycles, X, markers=True, label=r"$\chi$")
    plt.xlabel(r"$MCMC$ $cycles$")
    plt.ylabel(r"$\chi$ $[k_b^{-1}]$")
    plt.xscale("log")
    plt.savefig("../figures/%s_X.pdf" %(savename), dpi=300, bbox_inches='tight')
    plt.show()

def main():
    data_1 = pa.mat()
    data_2 = pa.mat()
    data_3 = pa.mat()
    data_4 = pa.mat()
    data_5 = pa.mat()
    data_6 = pa.mat()

    data_6.load("cycles_L_20_T_2.4_order.dat")
    data_5.load("cycles_L_20_T_1_order.dat")
    data_4.load("cycles_L_20_T_2.4.dat")
    data_3.load("cycles_L_20_T_1.dat")
    data_2.load("data/T_val_1mill.dat")
    data_1.load("data/cycle_val.dat")
    #plot_data(np.array(data_3), 1, "random", savename="numeric_L_20_T_1")
    #plot_data(np.array(data_4), 2.4, "random", savename="numeric_L_20_T_2.4")
    #plot_data(np.array(data_5), 1, "ordered", savename="numeric_L_20_T_1_order")
    #plot_data(np.array(data_6), 2.4, "ordered", savename="numeric_L_20_T_2.4_order")


    #plot_diff(np.array(data_1))

    plot_T(np.array(data_2))

if __name__ == "__main__":
    main()
