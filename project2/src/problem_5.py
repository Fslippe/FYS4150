import os
import numpy as np
import pyarma as pa
import matplotlib.pyplot as plt

def run(list, compile=True):
    """
    Compiling c++ file if compile = true
    and running jacobi rotation algorithm for different matrixsizes.
    Saving files as "N.dat" and "iterations.dat" in data folder
    """
    if compile == True:
        print("Compiling...\n")
        compile_cpp = os.system("g++ problem_5.cpp functions.cpp -o problem5 -larmadillo"),
    run_cpp = os.system("./problem5 " + " ".join(map(str, list)))

def main():
    """setting up array with different matrix sizes"""
    min_matrix_size = 5
    max_matrix_size = 100
    stepsize = 5
    list = np.arange(min_matrix_size, max_matrix_size, stepsize)
    run(list=list, compile=False)

    """Importing data and plotting"""
    N = pa.mat()
    N.load("data/N.dat")
    iterations = pa.mat()
    iterations.load("data/iterations.dat")
    plt.plot(N, iterations)
    plt.xlabel("Matrix size (N)")
    plt.ylabel("Iterations")
    plt.xscale("log")
    plt.yscale("log")
    #plt.savefig("../figures/N_iter_log.pdf")
    #plt.savefig("../figures/N_iter_log_dense.pdf")
    plt.show()

if __name__ == '__main__':
    main()
