#include <assert.h>
#include <armadillo>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

//Create a dense symetric matrix with random elements
arma::mat create_dense_symetric(int N);

// Create a symetric tridiagonal matrix with size N x N
// diagonal d
// sub- and super-diagonal a
arma::mat create_symmetric_tridiagonal(int N, double a, double d);

// calculates eigenvectors for a symmetric tridiagonal matrix with size N x N
// diagonal d
// sub- and super-diagonal a
arma::mat analytic_eigenvector(int N, double a, double d);

// Calculates eigenvalues for a symetric tridiagonal matrix with size N x N
// diagonal d
// sub- and super-diagonal a
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
void jacobi_eigensolver(arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors,
                        const int maxiter, int& iterations, bool& converged);


// Sorting and normalising eigenvalues and eigenvectors
// using arma::sort and arma::normalise
void sort_normalise(arma::vec& eigval, arma::mat& eigvec);
