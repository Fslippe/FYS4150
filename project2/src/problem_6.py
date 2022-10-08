import os
import numpy as np
import pyarma as pa
import matplotlib.pyplot as plt

def run(n, compile=True):
    """
    Input n number of discretization steps
    Compiling c++ file if compile = true
    and running jacobi rotation algorithm matrix size n-1
    Saving x, analytic eigenvectors and jacobi eigenvectors and eigenvalues
    """
    if compile == True:
        print("Compiling...\n")
        compile_cpp = os.system("g++ problem_6.cpp functions.cpp -o problem6 -larmadillo"),
    run_cpp = os.system("./problem6 " + "%i" %(n))

def plot_compare(n):
    """
    input n number of discretization steps
    Loading datafiles and plotting the solutions for the 3 smallest eigenvalues
    Saving plots in 3 different files
    """
    v = pa.mat()
    v.load("data/v.dat")
    u = pa.mat()
    u.load("data/u.dat")
    x = pa.mat()
    x.load("data/x.dat")
    lmb = pa.mat()
    lmb.load("data/eigval.dat")

    for i in range(3):
        plt.figure()
        plt.title("$n=$%i     $\lambda=$%.2f" %(n, lmb[i]))
        plt.plot(x, u[:,i], label="analytic")
        plt.plot(x, v[:,i], "--", label="jacobi")
        plt.legend(loc="upper right")
        plt.xlabel("$\hat{x}$")
        plt.ylabel("eigenvector")
        plt.savefig("../figures/eigvec_%i_%i.pdf" %(n, i))
    plt.show()

def main():
    n = 10
    run(n, compile=True)
    plot_compare(n)

    n = 100
    run(n, compile=False)
    plot_compare(n)

if __name__ == '__main__':
    main()
