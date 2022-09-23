#include "functions.hpp"
// Compile and linking with:
// g++ problem_4.cpp functions.cpp -o problem4 -larmadillo

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
