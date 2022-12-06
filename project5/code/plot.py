import pyarma as pa 
import numpy as np 
import matplotlib.pyplot as plt 
import os 

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
        compile_cpp = os.system("g++ main.cpp double_slit_box.cpp -o main -larmadillo "),

    run_cpp = os.system("./main %s %s %s %s %s %s %s %s %s %s" %(h_in, dt_in, T_in, xc_in, sigma_x_in, px_in, yc_in, sigma_y_in, py_in, v0_in))   


def main():
    runcpp = True
    if runcpp:
        run(h_in=0.005,
            dt_in=2.5*10**(-5),
            #T_in=0.008,
            T_in=0.0001,
            xc_in=0.25,
            sigma_x_in=0.05,
            px_in=200,
            yc_in=0.5,
            sigma_y_in=0.05,
            py_in=0,
            v0_in=0,
            compile=True)
    A = pa.mat()
    #A.load("data/test.dat")
    #print(A)
    #plt.contourf(np.array(A).T)
    #plt.show()
if __name__ == "__main__":
    main()