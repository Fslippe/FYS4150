import pyarma as pa 
import numpy as np 
import matplotlib.pyplot as plt 
import os 
import time 
plt.rcParams.update({"font.size": 12})

def run(h_in,
    dt_in, T_in, xc_in, sigma_x_in, px_in, yc_in, sigma_y_in, py_in, v0_in, compile=False):
    """
    Compile and run c++ code for different parameters
    - h_in                  
    - dt_in
    - T_in
    - xc_in
    - sigma_x_in
    - px_in
    - yc_in
    - sigma_y_in
    - py_in
    - v0_in         
    - compile=False         if True compile c++ file 
    """

    if compile == True:
        print("Compiling...\n")
        compile_cpp = os.system("g++ main.cpp double_slit_box.cpp -o main -larmadillo -O2"),
    start = time.time()
    run_cpp = os.system("./main %s %s %s %s %s %s %s %s %s %s" %(h_in, dt_in, T_in, xc_in, sigma_x_in, px_in, yc_in, sigma_y_in, py_in, v0_in))   
    print("time: ",time.time()-start)
    
def plot(t, data, save):
    p_sum  = np.zeros(len(t))
    for i in range(len(t)):
        p_sum[i] = np.sum(data[:,i,:])

    plt.plot(t, p_sum-1)
    plt.ylabel(r"$\sum_{i,j} u_{ij}^{t*} u_{ij}^t -1$")
    plt.xlabel("t")
    plt.savefig("../figures/%s.pdf" %(save), dpi=300, bbox_inches="tight")
    plt.show()

def plot_time(t, dt, data, V, save):
    xy = np.linspace(0, 1, len(data[:,0,0]))
    x, y = np.meshgrid(xy, xy)
    for t_i in t:
        plt.figure()
        plt.title(r"$t=$%.3f" %(t_i))
        plt.contourf(x,y,data[:,int(t_i/dt),:]) 
        plt.colorbar()
        plt.contourf(x,y,V.T, levels=np.linspace(1e9, 1e10,2), cmap="gray") 
        plt.savefig("../figures/%s_%.3f.pdf" %(save, t_i), dpi=300, bbox_inches="tight")
        plt.xlabel("x")
        plt.ylabel("y")

    
    plt.show()

def main():
    runcpp = False
    dt = 2.5*10**(-5)

    if runcpp:
        run(h_in=0.005,
            dt_in=dt,
            T_in=0.008,
            #T_in=0.0002,
            xc_in=0.25,
            sigma_x_in=0.05,
            px_in=200,
            yc_in=0.5,
            sigma_y_in=0.2,
            py_in=0,
            v0_in=1e10,
            compile=True)


    t_points = np.arange(0, 1+dt, dt)
    DS = pa.cube()
    NS = pa.cube()
    DS2 = pa.cube()
    DS2_real = pa.cube()
    DS2_imag = pa.cube()
    V = pa.mat()

    t = pa.mat()
    t.load("data/t8.dat")
    DS.load("data/double_slit_sigmay_01.dat")
    V.load("data/V.dat")

    NS.load("data/no_slit_test.dat")
    DS2.load("data/double_slit_sigmay_02.dat")
    DS2_real.load("data/double_slit_sigmay_02_real.dat")
    DS2_imag.load("data/double_slit_sigmay_02_imag.dat")


    plot_time([0, 0.001, 0.002], dt, np.array(DS2), np.array(V), "aas")

    plot(np.array(t), np.array(DS), "double_slit_p_diff")
    plot(np.array(t), np.array(NS), "no_slit_p_diff")

    
if __name__ == "__main__":
    main()