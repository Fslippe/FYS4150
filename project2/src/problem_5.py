import os
import numpy as np
import pyarma as pa
import matplotlib.pyplot as plt

def run(list, matrix="dense", compile=True):
    """
    Compiling c++ file if compile = true
    and running jacobi rotation algorithm for different matrixsizes.
    Saving files as "N.dat" and "iterations.dat" in data folder
    """
    if compile == True:
        print("Compiling...\n")
        compile_cpp = os.system("g++ problem_5.cpp functions.cpp -o problem5 -larmadillo"),
    run_cpp = os.system("./problem5 " + " ".join(map(str, list)) + " %s" %(matrix))

def plot(label, line = "-"):
    """Importing data and plotting with given label and line"""
    N = pa.mat()
    N.load("data/N.dat")
    iterations = pa.mat()
    iterations.load("data/iterations.dat")
    plt.plot(N, iterations, "%s" %(line), label="%s" %(label))
    plt.xlabel("Matrix size (N)")
    plt.ylabel("Iterations")
    plt.xscale("log")
    plt.yscale("log")

def main():
    """setting up array with different matrix sizes"""
    min_matrix_size = 5
    max_matrix_size = 100
    stepsize = 5
    list = np.arange(min_matrix_size, max_matrix_size, stepsize)

    run(list=list, matrix="tri", compile=True)
    plot("Tridiagonal")
    #plt.legend()
    #plt.savefig("../figures/N_iter_log.pdf")
    #plt.show()

    run(list=list, compile=False)
    plot("Dense", line="--")
    plt.legend()
    #plt.savefig("../figures/N_iter_log_dense.pdf")
    plt.savefig("../figures/N_iter_log_both.pdf")
    plt.show()
if __name__ == '__main__':
    main()
