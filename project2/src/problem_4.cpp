#include "functions.hpp"
// Compile and linking with:
// g++ problem_4.cpp functions.cpp -o problem4 -larmadillo

int main()
{
  // set up tridiagonal A for N=6
  int N = 6;
  int n = N + 1;
  double h = 1./n;
  double a = -1./(h*h);
  double d = 2./(h*h);
  arma::mat A = create_symmetric_tridiagonal(N, a, d);

  // solve eigenvaule problem using Jacobi roatation algorithm
  arma::vec eigval = arma::vec(N);
  arma::mat eigvec = arma::mat(N,N);
  double eps = 1e-8;
  int maxiter = 1e6;
  int iterations = 0;
  bool converged;

  jacobi_eigensolver(A, eps, eigval, eigvec, maxiter, iterations, converged);

  // arma::normalise for comparing
  eigval = arma::normalise(eigval);
  eigvec = arma::normalise(eigvec);
  arma::vec analytic_eigval = arma::normalise(analytic_eigenval(N, a, d));
  arma::mat analytic_eigvec = arma::normalise(analytic_eigenvector(N, a, d));

  //Print eigenvalues and eigenvectors to see of Armadillo agrees with the analytical result
  std::cout << "JACOBI EIGENVALUE\n";
  eigval.print();
  std::cout << arma::endl;
   std::cout << "ANALYTIC EIGENVALUE\n";
  analytic_eigval.print();
  std::cout << arma::endl;

  std::cout << "JACOBI EIGENVECTOR\n";
  eigvec.raw_print();
  std::cout << arma::endl;
   std::cout << "ANALYTIC EIGENVECTOR\n";
  analytic_eigvec.raw_print();
  std::cout << arma::endl;

  std::cout << "Iterations:\n";
  std::cout << iterations <<std::endl;
  //eigvec = abs(eigvec) - abs(analytic_eigvec);

  // bool test
  std::cout << "Checking if armadillo eigenvectors match analytic soulutions\n";
  for (int i = 0; i < N; i++) {
    assert((fabs(eigval(i)) - fabs(analytic_eigval(i)) < 1e-15));
    for (int j = 0; j < N; j++) {
      assert(fabs(eigvec(i,j)) - fabs(analytic_eigvec(i,j)) < 1e-15);
    }
  }
  std::cout << "Check successfull\n";
 return 0;
}
