#include <iostream>
#include <armadillo>
#include <string>
#include <fstream>
#include <iomanip>

int main() {
  int n = 1000;
  // Defining a vector of x values.
  arma::vec x = arma::linspace(0,1,n+1); //Declare and will with random uniform values.
  arma::vec u = arma::vec(n+1); //Initialize vector but don't fill.

  double u_func(double x); // Declaration of u(x).

  // Filling u with u(x) values for x in x.
  for (int i = 1; i <= n-1; i++)
  {
    u[i] = u_func(x[i]);

  }
  u[0] = 0;
  u[n] = 0;
  arma::mat xu = arma::mat(n+1, 2);
  xu.col(0) = x;
  xu.col(1) = u;
  xu.save("data/x_u1000.dat");
  std::cout << xu;

  return 0;
}

// u(x)
double u_func(double x){
    return 1 - (1 - exp(-10))*x - exp(-10*x);
}
