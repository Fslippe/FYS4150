#include <armadillo>

arma::mat create_tridiagonal(int N, double a, double d, double e);
arma::mat create_symmetric_tridiagonal(int N, double a, double d);

int main()
{
// set up tridiagonal A for N=6
int N = 6
int n = N + 1;
double h = 1./n
double a = -1./(h*h)
double d = 2./(h*h)
arma::mat A = create_symmetric_tridiagonal(N,a,d);

// solve eigenvaule problem using Armadillo’s arma::eig_sym

//checks that the eigenvalues and eigenvectors from Armadillo agrees with the analytical result for N=6
// arma::normalise when comparing
  return 0;
}

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
// not done
arma::mat analytic_eigenvec(int N, double a, double d)
{
  const double pi = 2*acos(0.0);
  arma::mat eig_mat = arma::mat(N,N); 
  for(int i = 1; i <= N; ++i)
  {
        eig_vec.col(i-1) = ;
  }
  return eig_vec;
}
