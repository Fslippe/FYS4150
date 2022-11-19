import seaborn as sns
import numpy as np 
import pyarma as pa 
import matplotlib.pyplot as plt 
import os
import time
sns.set_style("darkgrid")

def run(threads, T, lattice_dim, order, T_min=0.5, T_max=4, cycles=False, temp=False, hist=False, compile=False):
    """
    Compile and run c++ code for different parameters
    - threads               Threads to use in parallelization
    - T                     Temperature
    - lattice_dim           Dimension of lattice
    - order                 ordered initial lattice spins (true), random (false)
    - cycles=False          if True run cycle loop for chosen parameters 
    - temp=False            if True run temperature loop for chosen parameters 
    - hist=False            if True run histogram for chosen parameters 
    - compile=False         if True compile c++ file 
    """
    if compile == True:
        print("Compiling...\n")
        compile_cpp = os.system("g++ main.cpp ising_model.cpp -o main -larmadillo -fopenmp -O"),
    
    if cycles:
        save_cycles = "data/cycles_L_%i_T_%.1f_%s.dat" %(lattice_dim, T, order)
    else:
        save_cycles = "none"

    if temp:
        save_temp = "data/temp_L_%i_T_%.1f_%s.dat" %(lattice_dim, T, order)
    else:
        save_temp = "none"

    if hist:
        save_hist = "data/histogram_L_%i_T_%.1f_%s.dat" %(lattice_dim, T, order)
    else:
        save_hist = "none"

    print("Running")
    run_cpp = os.system("./main %s %s %s %s %s %s %s %s %s" %(threads, T, lattice_dim, order, save_cycles, save_temp, save_hist, T_min, T_max))   


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

    if plot_analytic:
        save = "numeric_analytic_"
        
    else:
        save = ""
        

    plt.figure()
    if plot_analytic:
        sns.lineplot(x=t, y=e, color="r", label="Analytic")
    plt.xlabel(r"$T$ [ $J/k_b$ ]")
    plt.ylabel(r"$\left<\epsilon\right>$ [ $J$ ]")
    sns.scatterplot(x=t, y=data[1,:],  label="Numeric")
    plt.savefig("../figures/%se_T.pdf" %(save), dpi=300, bbox_inches='tight')

    plt.figure()
    if plot_analytic:
        sns.lineplot(x=t, y=m, color="r", label="Analytic")
    sns.scatterplot(x=t, y=data[2,:],  label="Numeric")
    plt.xlabel(r"$T$ [ $J/k_b$ ]")
    plt.ylabel(r"$\left< m \right>$ [ $J$ ]")
    plt.savefig("../figures/%sm_T.pdf" %(save), dpi=300, bbox_inches='tight')

    plt.figure()
    if plot_analytic:
        sns.lineplot(x=t, y=C_v, color="r", label="Analytic")
    sns.scatterplot(x=t, y=data[3,:],  label="Numeric")
    plt.xlabel(r"$T$ [ $J/k_b$ ]")
    plt.ylabel(r"$C_v$ [ $k_b$ ]")
    plt.savefig("../figures/%sc_v_T.pdf"%(save), dpi=300, bbox_inches='tight')

    plt.figure()
    if plot_analytic:
        sns.lineplot(x=t, y=X, color="r", label="Analytic")
    sns.scatterplot(x=t, y=data[4,:], label="Numeric")
    plt.xlabel(r"$T$ [ $J/k_b$ ]")
    plt.ylabel(r"$\chi$ $[k_b^{-1}]$")
    plt.savefig("../figures/%sX_T.pdf" %(save), dpi=300, bbox_inches='tight')

    plt.show()

def plot_diff(data, T):
    E, E_2, e, M, M_2, m, C_v, X = analytic(T) 
    print("analytic  ", C_v)
    cycles = data[0,:]
    #print(cycles)
    e_diff = np.abs(data[1,:] - e)
    print(data[4,:])
    m_diff = np.abs(data[2,:] - m)
    C_v_diff = np.abs(data[3,:] - C_v)
    X_diff = np.abs(data[4,:] - X)

    plt.figure()
    sns.lineplot(x=cycles, y=e_diff, markers=True, label=r"$|\epsilon_{Analytic}-\epsilon_{Numeric}|$")
    #plt.xlabel(r"$MCMC$ $cycles$")
    #plt.ylabel(r"$\Delta \left<\epsilon\right>$ [ $J$ ]")
    #plt.xscale("log")
    #plt.yscale("log")
    #plt.savefig("../figures/numeric_diff_e.pdf", dpi=300, bbox_inches='tight')
    #plt.figure()

    sns.lineplot(x=cycles, y=m_diff, markers=True, label=r"$|m_{Analytic}-m_{Numeric}|$")
    #plt.xlabel(r"$MCMC$ $cycles$")
    #plt.ylabel(r"$\Delta \left< m \right>$ [ $J$ ]")
    #plt.xscale("log")
    #plt.yscale("log")
    #plt.savefig("../figures/numeric_diff_m.pdf", dpi=300, bbox_inches='tight')
    #plt.figure()

    sns.lineplot(x=cycles, y=C_v_diff, markers=True, label=r"$|C_{v,Analytic}-C_{v,Numeric}|$")
    #plt.xlabel(r"$MCMC$ $cycles$")
    #plt.ylabel(r"$\Delta C_v$ [ $k_b$ ]")
    #plt.xscale("log")
    #plt.yscale("log")
    #plt.savefig("../figures/numeric_diff_c_v.pdf", dpi=300, bbox_inches='tight')
    #plt.figure()

    sns.lineplot(x=cycles, y=X_diff, markers=True, label=r"$|\chi_{Analytic}-\chi_{Numeric}|$")
    #plt.xlabel(r"$MCMC$ $cycles$")
    plt.ylabel(r"$\Delta (numerical analytical)$")
    plt.xscale("log")
    #plt.yscale("log")
    plt.savefig("../figures/numeric_analytic.pdf", dpi=300, bbox_inches='tight')
    plt.show()

def plot_data(data, data_order, T, savename):
    cycles = data_order[0,:]
    e_order = data_order[1,:]
    m_order = data_order[2,:]
    C_v_order = data_order[3,:]
    X_order = data_order[4,:]
    e = data[1,:]
    m = data[2,:]
    C_v = data[3,:]
    X = data[4,:]


    plt.figure()
    plt.title(r"$T=$ %.1f $J/k_b$" %(T))
    sns.lineplot(x=cycles, y=e, markers=True, label=r"$\epsilon$")
    sns.lineplot(x=cycles, y=e_order, linestyle="--", markers=True, label=r"$\epsilon_{order}$")
    plt.xlabel(r"$MCMC$ $cycles$")
    plt.ylabel(r"$\left<\epsilon\right>$ [ $J$ ]")
    plt.xscale("log")
    plt.savefig("../figures/%s_e.pdf" %(savename), dpi=300, bbox_inches='tight')

    plt.figure()
    plt.title(r"$T=$%.1f $J/k_b$" %(T))
    sns.lineplot(x=cycles, y=m, markers=True, label=r"$m$")
    sns.lineplot(x=cycles, y=m_order, linestyle="--", markers=True, label=r"$m_{order}$")
    plt.xlabel(r"$MCMC$ $cycles$")
    plt.ylabel(r"$\left< m \right>$ [ $J$ ]")
    plt.xscale("log")
    plt.savefig("../figures/%s_m.pdf" %(savename), dpi=300, bbox_inches='tight')

    plt.show()

def plot_hist(data, savename):
    print("min: ", np.min(data[:,0]))
    print("max: ", np.max(data[:,0]))
    print("mean: ", np.mean(data[:,0]))
    sns.histplot(data[:,0], stat="probability", bins=33)
    plt.xlabel(r"$\epsilon$ [ $J$ ]")
    plt.savefig("../figures/%s_m.pdf" %(savename), dpi=300, bbox_inches='tight')

    plt.show()

def timing_test():
    print("\n\nTime used with 1 thread:\n")
    thread1 = time.time()
    os.system("./main %s %s %s %s %s %s %s" %(1, 1, 5, True, "none", "test", "none"))
    thread1_total = time.time() - thread1

    time.sleep(10) # to cool down computer for fair compairison
    print("\n\nTime used with 8 threads:\n")
    thread8 = time.time()
    os.system("./main %s %s %s %s %s %s %s" %(8, 1, 5, True, "none", "test", "none"))
    thread8_total = time.time() - thread8
    time.sleep(10) # to cool down computer for fair compairison
    print("\n\nTime used with 8 threads:\n")
    thread4 = time.time()
    os.system("./main %s %s %s %s %s %s %s" %(4, 1, 5, True, "none", "test", "none"))
    thread4_total = time.time() - thread4
    print("\n Speed up factor from 1 to 8 threads: ", thread8_total/thread1_total)
    print("\n Speed up factor from 1 to 4 threads: ", thread4_total/thread1_total)
    print("\n Speed up factor from 4 to 8 threads: ", thread8_total/thread4_total)



def main():
    cycle_L_2_1 = pa.mat()
    temp_L_2 = pa.mat()
    cycle_L20_1 = pa.mat()
    cycle_L20_2_4 = pa.mat()
    cycle_L20_1_order = pa.mat()
    cycle_L20_2_4_order = pa.mat()
    histogram_T_1 = pa.mat()
    histogram_T_2_4 = pa.mat()
    test = pa.mat()

    run_all = False
    run_long = False

    if run_all:
        run(threads=8, T=2.4, lattice_dim=2, order="true", cycles=True, compile=True)
        run(threads=8, T=2.4, lattice_dim=2, order="false", cycles=True, temp=True)
        run(threads=8, T=1.0, lattice_dim=2, order="true", cycles=True)
        run(threads=8, T=1.0, lattice_dim=2, order="false", cycles=True)
        run(threads=8, T=2.4, lattice_dim=20, order="false", cycles=True, hist=True)
        run(threads=8, T=1.0, lattice_dim=20, order="false", cycles=True, hist=True)
        run(threads=8, T=2.4, lattice_dim=20, order="true", cycles=True)
        run(threads=8, T=1.0, lattice_dim=20, order="true", cycles=True)
    #run(threads=8, T=2.4, lattice_dim=20, order="false", cycles=True, hist=True, compile=True)
    #run(threads=8, T=1.0, lattice_dim=20, order="false", cycles=True, hist=True)
    
    if run_long:
        print("Running lattice 40")
        test = time.time()
        run(threads=8, T=1.0, lattice_dim=40, order="false", temp=True, T_min=2.1, T_max = 2.4)
        print(time.time()-test)

        print("Running lattice 60")
        test = time.time()
        run(threads=8, T=1.0, lattice_dim=60, order="false", temp=True, T_min=2.1, T_max = 2.4)
        print(time.time()-test)

        print("Running lattice 80")
        test = time.time()
        run(threads=8, T=1.0, lattice_dim=80, order="false", temp=True, T_min=2.1, T_max = 2.4)
        print(time.time()-test)

        print("Running lattice 100")
        test = time.time()
        run(threads=8, T=1.0, lattice_dim=100, order="false", temp=True, T_min=2.1, T_max = 2.4)
        print(time.time()-test)



    # cycle loop T = 2.4 
    cycle_L20_2_4_order.load("data/cycles_L_20_T_2.4_true.dat")
    cycle_L20_2_4.load("data/cycles_L_20_T_2.4_false.dat")
    
    # cycle loop T = 1.0 
    cycle_L20_2_4_order.load("data/cycles_L_20_T_2.4_true.dat")
    cycle_L20_2_4.load("data/cycles_L_20_T_2.4_false.dat")
    histogram_T_1.load("data/histogram_L_20_T_1.0_false.dat")
    histogram_T_2_4.load("data/histogram_L_20_T_2.4_false.dat")
    cycle_L20_1_order.load("data/cycles_L_20_T_1.0_true.dat")
    cycle_L20_1.load("data/cycles_L_20_T_1.0_false.dat")
    temp_L_2.load("data/temp_L_2_T_2.4_false.dat")
    cycle_L_2_1.load("data/cycles_L_2_T_1.0_false.dat")

    test.load("data/temp_L_3_T_1.0_false.dat")

    plot_data(np.array(cycle_L20_1), np.array(cycle_L20_1_order), 1, savename="numeric_L_20_T_1")
    plot_data(np.array(cycle_L20_2_4), np.array(cycle_L20_2_4_order), 2.4, savename="numeric_L_20_T_2_4")
    # Histogram T

    plot_hist(np.array(histogram_T_1), "histogram_T_1")
    plot_hist(np.array(histogram_T_2_4),"histogram_T_2_4")

    plot_data(np.array(cycle_L20_1), np.array(cycle_L20_1_order), 1, savename="numeric_L_20_T_1")

    plot_diff(np.array(cycle_L_2_1), 1)
    plot_T(np.array(temp_L_2))
    #timing_test()


if __name__ == "__main__":
    main()
