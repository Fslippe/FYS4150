#include "functions.hpp"
// g++ problem_5.cpp functions.cpp -o problem5 -larmadillo
// This file requires commandline arguments to run
// to run for any integer like N=10 and N=20 run ./problem5 10 20

int main(int argc, char** argv)
{
  //lists for matrix size N and iterations
  arma::vec N_list = arma::vec(argc-1);
  arma::vec iter_list = arma::vec(argc-1);

  for (int i = 1; i < argc; i++) {
    //setting up the A matrix for N equal commandline argument i
    int N = atoi(argv[i]);
    int n = N + 1;
    double h = 1./n;
    double a = -1./(h*h);
    double d = 2./(h*h);

    //setting up either dense or triagonal symetric matrix
    //arma::mat A = create_dense_symetric(N);
    arma::mat A = create_symmetric_tridiagonal(N, a, d);

    // solve eigenvaule problem using Jacobi roatation algorithm
    arma::vec eigval = arma::vec(N);
    arma::mat eigvec = arma::mat(N,N);
    double eps = 1e-9;
    int maxiter = N*N*N;
    int iterations = 0;
    bool converged;
    jacobi_eigensolver(A, eps, eigval, eigvec, maxiter, iterations, converged);

    //saving iterations and matrix size to list
    N_list[i-1] = N;
    iter_list[i-1] = iterations;
  }
  //saving iteration and matrix list to files
  N_list.save("data/N.dat");
  iter_list.save("data/iterations.dat");

  return 0;
}
