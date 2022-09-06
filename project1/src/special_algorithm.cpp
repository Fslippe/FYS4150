#include <iostream>
#include <armadillo>
#include <string>
#include <fstream>
#include <iomanip>

arma::vec general_algorithm(arma::vec v, arma::vec g); // Declaration of u(x).


int main() {
  int n = 10; // Number of steps
  double v0 = 0.;
  double vn = 0.;

  arma::vec x = arma::linspace(0, 1, n+1); //Declare and will with random uniform values.
  arma::vec g = arma::vec(n-1);
  arma::vec v = arma::vec(n+1);
  arma::vec u = arma::vec(n+1);
  double h = x[1] - x[0];

  for (int i = 0; i <= n-2; i++)
  {
    g[i] = 100*exp(-10*x[i+1])*h*h;
  }

  v[0] = 0;
  v[n-1] = 0;

  v = general_algorithm(v, g);
//            -- Output to a datafile --
// filename
  arma::mat data = arma::mat(n+1, 2);
  data.col(0) = x;
  data.col(1) = v;

  //data.save("n10.dat");
  return 0;
}

arma::vec general_algorithm(arma::vec v, arma::vec g)
{
  int a = -1;
  int b = 2;
  int c = -1;
  double b_tilde = b - a/b*c;
  double ab_tilde = a / b_tilde;

  int n_matrix = g.size();

  arma::vec g_tilde = arma::vec(n_matrix);
  b_tilde[0] = b[0];
  g_tilde[0] = g[0];

  for (int i = 1; i <= n_matrix-1; i++)
  {
    g_tilde[i] = g[i] - ab_tilde * g_tilde[i-1];
  }

  // Different index for v becuase of different length
  v[n_matrix] = g_tilde[n_matrix-1] / b[n_matrix-1];

  for (int i = n_matrix-1; i >= 1; i--)
  {
    v[i] = (g_tilde[i-1] - c * v[i+1]) / b_tilde;
  }
return v;
}
