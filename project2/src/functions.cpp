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

//creates a dense symetric matrix with random numbers
arma::mat create_dense_symetric(int N)
{
  // Generate random N*N matrix
  arma::mat A = arma::mat(N, N).randn();
  // Symmetrize the matrix by reflecting the upper triangle to lower triangle
  A = arma::symmatu(A);

  return A;
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

  double maxval = 0.; //std::abs(A(k,l)); //(i,j) (row, col)

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

 // Performs a single Jacobi rotation, to "rotate away"
 void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l)
 {
    double tau = (A(l,l) - A(k,k)) / (2. * A(k,l));
    int N = A.n_rows;
    double s, c, t;
    double tmp_element;

    if (A(k,l) != 0.)
    {
        if (tau > 0)
        {
          t = - tau + sqrt(1. + tau*tau);
        }
        else
        {
          t = - tau - sqrt(1. + tau*tau);
        }
        c = 1 / sqrt(1. + t*t);
        s = c * t;
    }
    else
    {
        c = 1.;
        s = 0.;
    }

    tmp_element = A(k,k);
    A(k,k) = A(k,k)*c*c - 2*A(k,l)*c*s + A(l,l)*s*s;
    A(l,l) = A(l,l)*c*c + 2*A(k,l)*c*s + tmp_element*s*s;
    A(k,l) = 0;
    A(l,k) = 0;

    //A.print();
    //std::cout << "k: " << k << ", l: " << l << std::endl;
    for (int i = 0; i < N; i++)
    {
        if(i != k && i != l)
        {
            tmp_element = A(i,k);
            A(i,k) = tmp_element*c - A(i,l)*s;
            A(k,i) = A(i,k);
            A(i,l) = A(i,l)*c + tmp_element*s;
            A(l,i) = A(i,l);
        }

        tmp_element = R(i,k);
        R(i,k) = tmp_element*c - R(i,l)*s;
        R(i,l) = R(i,l)*c + tmp_element*s;
    }
    return;
 }

 // Jacobi method eigensolver runs jacobi_rotate until max in A is less than eps:
 void jacobi_eigensolver(arma::mat& A, double eps, arma::vec& eigval, arma::mat& eigvec, const int maxiter, int& iterations, bool& converged)
 {
    int n = A.n_rows;
    int k, l;

    arma::mat R = arma::mat(n, n, arma::fill::eye);
    double max = max_offdiag_symmetric(A, k, l);

    while(std::abs(max)>= eps && iterations < maxiter)
    {
      jacobi_rotate(A, R, k, l);
      iterations += 1.;
      max = max_offdiag_symmetric(A, k, l);
    }


    eigval = A.diag();
    eigvec = R;

    if (iterations < maxiter)
    {
      bool converged = true;
      std::cout << "converged\n";
    }
}


void sort_normalise(arma::vec& eigval, arma::mat& eigvec)
{
  eigval = sort(arma::normalise(eigval));
  eigvec = arma::normalise(eigvec);
  eigvec = eigvec.each_col( [](arma::vec& vec){vec = arma::conv_to<arma::vec>::from(arma::sort(vec)); } );
}
