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

// print to test
int main()
{
int n = 30;
arma::mat A = create_tridiagonal(n,6,1,9);
for (int i = 0; i < n; i++)
{
        for (int j = 0; j < n; j++) {
            std::cout << A(i,j) << ' ';
        }
        std::cout << std::endl;
}
  return 0;
}