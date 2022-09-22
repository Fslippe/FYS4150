#include <assert.h>
#include <armadillo>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include "functions.hpp"


arma::mat create_symmetric_tridiagonal(int N, double a, double d);
arma::mat analytic_eigenvector(int N, double a, double d);
arma::vec analytic_eigenval(int N, double a, double d);

// Determine the the max off-diagonal element of a symmetric matrix A
// - Saves the matrix element indicies to k and l 
// - Returns absolute value of A(k,l) as the function return value
double max_offdiag_symmetric(const arma::mat& A, int& k, int& l);

// Performs a single Jacobi rotation, to "rotate away"
// the off-diagonal element at A(k,l).
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l);

// Jacobi method eigensolver:
// - Runs jacobo_rotate until max off-diagonal element < eps
// - Writes the eigenvalues as entries in the vector "eigenvalues"
// - Writes the eigenvectors as columns in the matrix "eigenvectors"


//   (The returned eigenvalues and eigenvectors are sorted using arma::sort_index)
// - Stops if it the number of iterations reaches "maxiter"
// - Writes the number of iterations to the integer "iterations"
// - Sets the bool reference "converged" to true if convergence was reached before hitting maxiter
void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged);


int main()
{
  // set up tridiagonal A for N=6
  int N = 6;
  int n = N + 1;
  double h = 1./n;
  double a = -1./(h*h);
  double d = 2./(h*h);
  arma::mat A = create_symmetric_tridiagonal(N, a, d);

  // solve eigenvaule problem using Armadillo’s arma::eig_sym
  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval, eigvec, A);

  // arma::normalise for comparing
  eigval = arma::normalise(eigval);
  eigvec = arma::normalise(eigvec);
  arma::vec analytic_eigval = arma::normalise(analytic_eigenval(N, a, d));
  arma::mat analytic_eigvec = arma::normalise(analytic_eigenvector(N, a, d));

  //Print eigenvalues and eigenvectors to see of Armadillo agrees with the analytical result
  eigval.print();
  std::cout << arma::endl;
  analytic_eigval.print();
  std::cout << arma::endl;
  eigvec.raw_print();
  std::cout << arma::endl;
  analytic_eigvec.raw_print();
  std::cout << arma::endl;
  //eigvec = abs(eigvec) - abs(analytic_eigvec);
  eigvec.raw_print();

  // bool test
  std::cout << "Checking if armadillo eigenvectors match analytic soulutions\n";
  for (int i = 0; i < N; i++) {
    assert((abs(eigval(i)) - abs(analytic_eigval(i)) < 1e-15));
    for (int j = 0; j < N; j++) {
      assert(abs(eigvec(i,j)) - abs(analytic_eigvec(i,j)) < 1e-15);
    }
  }
  std::cout << "Check successfull\n";
 return 0;
}




