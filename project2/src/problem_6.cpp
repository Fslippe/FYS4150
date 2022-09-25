#include "functions.hpp"
// g++ problem_6.cpp functions.cpp -o problem6 -larmadillo

int main(int argc, char** argv)
{
  //Setting up data for commandlineargument n
  int n = atoi(argv[1]);
  int N = n-1;
  int N_eigval = 3; //Eigenvalues to keep for solution
  double h = 1./n;
  double a = -1./(h*h);
  double d = 2./(h*h);
  arma::mat A = create_symmetric_tridiagonal(N, a, d);
  arma::vec eigval = arma::vec(N);
  arma::mat eigvec = arma::mat(N,N);
  arma::vec idx = arma::vec(N_eigval);
  double eps = 1e-9;
  int maxiter = N*N*N;
  int iterations = 0;
  bool converged;

  //jacobi and analytic solver
  jacobi_eigensolver(A, eps, eigval, eigvec, maxiter, iterations, converged);
  arma::mat analytic_eigvec = analytic_eigenvector(N, a, d);

  // Finding index of 3 smallest eigenvalues
  // adding maximum to avoid finding same index several times
  double max = eigval.max();
  for (int i = 0; i < N_eigval; i++)
  {
    idx[i] = eigval.index_min();
    eigval[idx[i]] = eigval[idx[i]] + max;
  }

  // Returning the eigenvalues back to normal by subtracting maximum
  for (int i = 0; i < N_eigval; i++)
  {
    eigval[idx[i]] = eigval[idx[i]] - max;
  }

  // setting up arrays to save eigenvector and eigenvalue solutions
  arma::vec x = arma::linspace(0, 1, N+2);
  arma::mat v = arma::mat(N+2, N_eigval);
  arma::mat u = arma::mat(N+2, N_eigval);
  arma::vec lamda = arma::vec(N_eigval); //to store N_eigval smallest eigenvalues

  for (int i = 0; i < N_eigval; i++)
  {
    lamda[i] = eigval[idx[i]];

    for (int j = 0; j < N; j++)
    {
      v.col(i)[0] = 0;
      v.col(i)[j+1] = analytic_eigvec.col(idx[i])[j];
      u.col(i)[0] = 0;
      u.col(i)[j+1] = analytic_eigvec.col(idx[i])[j];
    }
  }

  // Saving datafiles
  v.save("data/v.dat");
  lamda.save("data/eigval.dat");
  u.save("data/u.dat");
  x.save("data/x.dat");

  return 0;
}
