#include "functions.hpp"
// Compile and linking with:
// g++ problem_5.cpp functions.cpp -o problem5 -larmadillo
// This file requires commandline arguments to run
// to run for any integer and a tridiagonal matrix  like N=10 and N=20:
// run ./problem5 10 20 tri
// Same N for dense matrix:
// run ./problem5 10 20 dense

int main(int argc, char** argv)
{
  //lists for matrix size N and iterations
  arma::vec N_list = arma::vec(argc-2);
  arma::vec iter_list = arma::vec(argc-2);
  std::string matrix;

  for (int i = 1; i < argc-1; i++) {
    //setting up the A matrix for N equal commandline argument i
    int N = atoi(argv[i]);
    int n = N + 1;

    arma::mat A;

    // setting up either dense or triagonal symetric matrix
    // if last commandlien argument equals tri: setting up tridiagonal
    // Else setting up dense matrix
    if (argv[argc-1] == std::string("tri"))
    {
      double h = 1./n;
      double a = -1./(h*h);
      double d = 2./(h*h);
      A = create_symmetric_tridiagonal(N, a, d);
      matrix = "Used Tridiagonal matrix\n";
    }
    else
    {
      A = create_dense_symetric(N);
      matrix = "Used Dense matrix\n";
    }

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
  std::cout << matrix;

  N_list.save("data/N.dat");
  iter_list.save("data/iterations.dat");

  return 0;
}
