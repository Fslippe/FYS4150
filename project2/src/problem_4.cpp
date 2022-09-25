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
  double eps = 1e-9;
  int maxiter = N*N*N;
  int iterations = 0;
  bool converged;

  jacobi_eigensolver(A, eps, eigval, eigvec, maxiter, iterations, converged);

  // The returned eigenvalues and eigenvectors are sorted using arma::sort
  sort_normalise(eigval, eigvec);
  
  arma::vec analytic_eigval = arma::normalise(sort(analytic_eigenval(N, a, d)));
  arma::mat analytic_eigvec = arma::normalise(analytic_eigenvector(N, a, d));
  analytic_eigvec = eigvec.each_col( [](arma::vec& vec){vec = arma::conv_to<arma::vec>::from(arma::sort(vec)); } );


  //Print eigenvalues and eigenvectors to see if jacobi agrees with the analytical result
  std::cout << "JACOBI EIGENVALUE\n";
  eigval.print();
  std::cout << arma::endl;
   std::cout << "ANALYTIC EIGENVALUE\n";
  analytic_eigval.print();
  std::cout << arma::endl;

  std::cout << "JACOBI EIGENVECTOR\n";
  eigvec.print();
  std::cout << arma::endl;
   std::cout << "ANALYTIC EIGENVECTOR\n";
  analytic_eigvec.print();
  std::cout << arma::endl;

  std::cout << "Iterations:\n";
  std::cout << iterations <<std::endl;
  //eigvec = abs(eigvec) - abs(analytic_eigvec);

  // bool test
  std::cout << "Checking if eigenvectors and eigenvalues match analytic soulutions\n";
  for (int i = 0; i < N; i++) {
    assert((fabs(eigval(i)) - fabs(analytic_eigval(i)) < 1e-15));
    for (int j = 0; j < N; j++) {
      assert(fabs(eigvec(i,j)) - fabs(analytic_eigvec(i,j)) < 1e-15);
    }
  }
  std::cout << "Check successfull\n";
 return 0;
}
