#include <assert.h>
#include <armadillo>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

arma::mat create_tridiagonal(int N, double a, double d, double e);
arma::mat create_symmetric_tridiagonal(int N, double a, double d);
arma::vec analytic_eigenval(int N, double a, double d);
arma::mat analytic_eigenvector(int N, double a, double d);

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




            
                // functions bellow this point 

// Creates a tridiagonal matrix
arma::mat create_tridiagonal(int N, double a, double d, double e)
{
  arma::mat A = arma::mat(N, N, arma::fill::eye) * d;
  A(0,1) = e; // first superdiagonal element

  for (int i = 1; i < N-1; i++)
  {
    A(i,i+1) = e; // filling superdiagonal
    A(i,i-1) = a; // filling subdiagonal
  }

  A(N-1,N-2) = a; //last subdiagonal element
  return A;
}

// Creates a symmetric tridiagonal matrix
arma::mat create_symmetric_tridiagonal(int N, double a, double d)
{
  return create_tridiagonal(N, a, d, a);
}

// Finding the max off-diagonal element of a matrix
double max_offdiag_symmetric(const arma::mat& A, int& k, int &l)
{
  int N = A.n_rows;
  assert(N > 1);
  assert(A.is_square());

  double maxval = std::abs(A(k,l)); //(i,j) (row, col)
  std::cout << maxval << "\n";

  int col_n = 1;

  for (int i = 1; i <= N-1; ++i) //Loops over rows
  {
    for (int j = 0; j < i; ++j) //loops over columns until the subdiagonal
    {
      if(std::abs(A(i,j)) > maxval)
      {
        k = i; l = j;
        maxval = std::abs(A(i,j));
      }
   }
  }
  return maxval;
 }

 // Calculates analytic eigenvalues
arma::vec analytic_eigenval(int N, double a, double d)
{
  const double pi = 2*acos(0.0);
  arma::vec eig_val = arma::vec(N);
  for(int i = 1; i <= N; ++i)
  {
    eig_val(i-1) = d + 2*a*cos(i*pi/(N+1));
  }
  return eig_val;
}

//Calculates analytic eigenvectors
arma::mat analytic_eigenvector(int N, double a, double d)
{
  const double pi = 2*acos(0.0);
  arma::mat eig_mat = arma::mat(N,N);
  for(int i = 1; i <= N; ++i)
  {
    arma::vec col_v = arma::vec(N);
    for(int j = 1; j <= N; ++j)
    {
      col_v(j-1) = sin((i*j*pi) / (N+1));
    }
    eig_mat.col(i-1) = col_v;
  }
  return eig_mat;
}
