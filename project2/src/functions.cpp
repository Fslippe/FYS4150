#include "functions.hpp"

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

// Finding the max off-diagonal element of a matrix
double max_offdiag_symmetric(const arma::mat& A, int& k, int &l)
{
  int N = A.n_rows;
  assert(N > 1);
  assert(A.is_square());

  double maxval = std::abs(A(k,l)); //(i,j) (row, col)

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

 void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l)
 {
    double tau = (A(l,l) - A(k,k)) / (2. * A(k,l));
    double A_ik;
    int N = A.n_rows;
    double t;
    if (tau > 0)
    {
       //t = - tau + sqrt(1. + tau*tau);
       t = 1. / (tau + sqrt(1 + tau*tau));
    }
    else
    {
        //t = - tau - sqrt(1. + tau*tau);
        t = -1. / (-tau + sqrt(1 + tau*tau));
    }
    double c = 1 / sqrt(1.+t*t);
    double s = c * t;

    double A_kk = A(k,k);
    A(k,k) = A(k,k)*c*c - 2*A(k,l)*c*s + A(l,l)*s*s;
    A(l,l) = A(l,l)*c*c + 2*A(k,l)*c*s + A_kk*s*s;
    A(k,l) = 0;
    A(l,k) = 0;

    for (int i = 0; i < N; i++)
    {
        if (i != k || i != l)
        {
            A_ik = A(i,k);
            A(i,k) = A(i,k)*c -A(i,l)*s;
            A(k,i) = A(i,k);
            A(i,l) = A(i,l)*c + A_ik*s;
            A(l,i) = A(i,l);
        }
        double R_ik = R(i,k);
        R(i,k) = R(i,k)*c - R(i,l)*s;
        R(i,l) = R(i,l)*c + R_ik*s;
    }
        
 }



 void jacobi_eigensolver(arma::mat& A, double eps, arma::vec eigval, arma::mat eigvec, const int maxiter, int iterations, bool converged)
 {
    int n = A.n_rows;
    int k=0; int l=1;

    arma::mat R = arma::mat(n, n, arma::fill::eye);
    double maxval = max_offdiag_symmetric(A, k, l);
    
    // test prints 
    R.print();
    std::cout << maxval << std::endl;

    // jacobi_rotate(A, R, k, l);
    // std::cout << "1 rot" << std::endl;
    // R.print();
    // A.print();
    // std::cout << maxval << std::endl;
    // jacobi_rotate(A, R, k, l);
    // std::cout << "2 rot" << std::endl;
    // R.print();
    //  A.print();
    // std::cout << maxval << std::endl;
    // jacobi_rotate(A, R, k, l);
    // std::cout << "3 rot" << std::endl;
    // R.print();
    // jacobi_rotate(A, R, k, l);
    // std::cout << "4 rot" << std::endl;
    // R.print();

    // actual loop
    while(std::abs(A(k,l))>= eps && iterations < maxiter)
    {
      jacobi_rotate(A, R, k, l);
      maxval = max_offdiag_symmetric(A, k, l);
      std::cout << "it's loopy" << std::endl;
      iterations += 1.;
      std::cout << iterations << std::endl;
    }
    R.print();
    eigvec = A;
    //eigvec.each_col( [](arma::vec& vec){vec = arma::conv_to<arma::vec>::from(arma::sort_index(vec)); } );
    eigval = arma::conv_to<arma::vec>::from(arma::sort_index(A.diag()));

    if (iterations < maxiter)
    {
      bool converged = true;
      std::cout << "converged\n"; 
    }
    

    }