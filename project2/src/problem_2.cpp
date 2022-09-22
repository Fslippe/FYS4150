#include <armadillo>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <assert.h>


arma::mat create_tridiagonal(int N, double a, double d, double e);
arma::mat create_symmetric_tridiagonal(int N, double a, double d);
arma::vec analytic_eigenval(int N, double a, double d);
arma::mat analytic_eigenvector(int N, double a, double d);


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

  // this bool test could work? not working
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
