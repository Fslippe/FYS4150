// One large file with all code neccesary for the jacobi rotation algorithm.
// debug effort since split file version was not working.
// All function names are capitalized to avoid linhking conflicts with other files (seemed to be an issue for some reason)

// Build: g++ -std=c++11 jacobi_rotation_algo_all_in_one.cpp functions.cpp -o jacobi -larmadillo
// Run: ./jacobi

#include <assert.h>
#include <armadillo>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

// Create a symetric tridiagonal matrix with size N x N
// diagonal d
// sub- and super-diagonal a
arma::mat Create_symmetric_tridiagonal(int N, double a, double d);

// calculates eigenvectors for a symmetric tridiagonal matrix with size N x N
// diagonal d
// sub- and super-diagonal a
arma::mat Analytic_eigenvector(int N, double a, double d);

// Calculates eigenvalues for a symetric tridiagonal matrix with size N x N
// diagonal d
// sub- and super-diagonal a
arma::vec Analytic_eigenval(int N, double a, double d);

// Determine the the max off-diagonal element of a symmetric matrix A
// - Saves the matrix element indicies to k and l
// - Returns absolute value of A(k,l) as the function return value
double Max_offdiag_symmetric(const arma::mat& A, int& k, int& l);

// Performs a single Jacobi rotation, to "rotate away"
// the off-diagonal element at A(k,l).
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R
void Jacobi_rotate(arma::mat& A, arma::mat& R, int& k, int& l);

// Jacobi method eigensolver:
// - Runs jacobo_rotate until max off-diagonal element < eps
// - Writes the eigenvalues as entries in the vector "eigenvalues"
// - Writes the eigenvectors as columns in the matrix "eigenvectors"


//   (The returned eigenvalues and eigenvectors are sorted using arma::sort_index)
// - Stops if it the number of iterations reaches "maxiter"
// - Writes the number of iterations to the integer "iterations"
// - Sets the bool reference "converged" to true if convergence was reached before hitting maxiter
void Jacobi_eigensolver(arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors,
                        const int maxiter, int& iterations, bool& converged);

// ----------------------------------------------------------------------------------------------------

int main()
{
  // set up tridiagonal A for N=6
  int N = 6;
  int n = N + 1;
  double h = 1./n;
  double a = -1./(h*h);
  double d = 2./(h*h);
  arma::mat A = Create_symmetric_tridiagonal(N, a, d);
  // solve eigenvaule problem using Jacobi roatation algorithm
  arma::vec eigval = arma::vec(N);
  arma::mat eigvec = arma::mat(N,N);
  double eps = 1e-8;
  int maxiter = N*N*N;
  int iterations = 0;
  bool converged;

  Jacobi_eigensolver(A, eps, eigval, eigvec, maxiter, iterations, converged);

  // The returned eigenvalues and eigenvectors are sorted using arma::sort_index
  eigval = arma::conv_to<arma::vec>::from(arma::sort_index(eigval));
  eigvec = eigvec.each_col( [](arma::vec& vec){vec = arma::conv_to<arma::vec>::from(arma::sort_index(vec)); } );

  // arma::normalise for comparing
  eigval = sort(arma::normalise(eigval)); //needs to be resorted for some reason
  eigvec = arma::normalise(eigvec);

  arma::vec analytic_eigval = arma::normalise(Analytic_eigenval(N, a, d));
  arma::mat analytic_eigvec = arma::normalise(Analytic_eigenvector(N, a, d));
  //analytic_eigvec = arma::normalise(analytic_eigvec.each_col( [](arma::vec& vec){vec = arma::conv_to<arma::vec>::from(arma::sort_index(vec)); } ));
  

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

// ----------------------------------------------------------------------------------------------------


// Creates a tridiagonal matrix
arma::mat Create_tridiagonal(int N, double a, double d, double e)
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
arma::mat Create_symmetric_tridiagonal(int N, double a, double d)
{
  return Create_tridiagonal(N, a, d, a);
}

// Calculates analytic eigenvalues
arma::vec Analytic_eigenval(int N, double a, double d)
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
arma::mat Analytic_eigenvector(int N, double a, double d)
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

// Finding the max off-diagonal element of a matrix
double Max_offdiag_symmetric(const arma::mat& A, int &k, int &l)
{
  int N = A.n_rows;
  assert(N > 1);
  assert(A.is_square());

  double maxval = 0.; //std::abs(A(k,l)); //(i,j) (row, col)

  int col_n = 1;

  for (int i = 1; i < N; ++i) //Loops over rows
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


 void Jacobi_rotate(arma::mat& A, arma::mat& R, int& k, int& l)
 {
    double tau = (A(l,l) - A(k,k)) / (2. * A(k,l));
    int N = A.n_rows;
    double s, c,t;
    if (A(k,l) != 0.)
    {
           if (tau > 0)
        {
        t = - tau + sqrt(1. + tau*tau);
        //t = 1. / (tau + sqrt(1 + tau*tau));
        }
        else
        {
            t = - tau - sqrt(1. + tau*tau);
            //t = -1. / (-tau + sqrt(1 + tau*tau));
        }
        c = 1 / sqrt(1.+t*t);
        s = c * t;
    }
    else
    {
        c = 1.;
        s = 0.;
    }

    double A_kk, A_ll, A_ik, A_il, R_ik, R_il;
    A_kk = A(k,k);
    A(k,k) = A_kk*c*c - 2*A(k,l)*c*s + A_ll*s*s;
    A(l,l) = A_ll*c*c + 2*A(k,l)*c*s + A_kk*s*s;
    A(k,l) = 0;
    A(l,k) = 0;
    
    //A.print();
    //std::cout << "k: " << k << ", l: " << l << std::endl;
    for (int i = 0; i < N; i++)
    {
        if(i != k && i != l)
        {
            //std::cout << "looping if i = " << i << std::endl;
            
            A_ik = A(i,k);
            A_il = A(i,l);
            A(i,k) = A_ik*c - A_il*s;
            A(k,i) = A(i,k);
            A(i,l) = A_il*c + A_ik*s;
            A(l,i) = A(i,l);
        }
    
        R_ik = R(i,k);
        R_il = R(i,l);
        R(i,k) = R_ik*c - R_il*s;
        R(i,l) = R_il*c + R_ik*s;
    }
    return;    
 }


 void Jacobi_eigensolver(arma::mat& A, double eps, arma::vec& eigval, arma::mat& eigvec, const int maxiter, int& iterations, bool& converged)
 {
    int n = A.n_rows;
    int k, l;

    arma::mat R = arma::mat(n, n, arma::fill::eye);
    double max = Max_offdiag_symmetric(A, k, l);
    
    while(std::abs(max)>= eps && iterations < maxiter)
    {
      max = Max_offdiag_symmetric(A, k, l);
      Jacobi_rotate(A, R, k, l);
      //std::cout << "it's loopy" << std::endl;
      iterations += 1.;
      //std::cout << iterations << std::endl;
    }
    eigval = A.diag();
    eigvec = R;

    if (iterations < maxiter)
    {
      bool converged = true;
      std::cout << "converged\n"; 
    }
    

    }