#include <armadillo>


arma::mat create_tridiagonal(int n, double a, double d, double e)
{
  arma::mat A = arma::mat(n, n, arma::fill::eye) * d;

  A(0,1) = e; // first superdiagonal element

  for (int i = 1; i < n-1; i++)
  {
    A(i,i+1) = e; // filling superdiagonal
    A(i,i-1) = a; // filling subdiagonal
  }

  A(n-1,n-2) = a; //last subdiagonal element 
  
  return A;
}
