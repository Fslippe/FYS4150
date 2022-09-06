#include <iostream>
#include <armadillo>
#include <string>
#include <fstream>
#include <iomanip>

arma::vec general_algorithm(arma::vec v, arma::vec g, arma::vec a, arma::vec b, arma::vec c); // Declaration of u(x).

int main() {
  int n = 10000; // Number of steps

  arma::vec x = arma::linspace(0, 1, n+1); //Declare and will with random uniform values.
  arma::vec g = arma::vec(n-1);
  arma::vec v = arma::vec(n-1);
  arma::vec v_full = arma::vec(n+1);
  arma::vec a = arma::vec(n-1).fill(-1);
  arma::vec b = arma::vec(n-1).fill(2);
  arma::vec c = arma::vec(n-1).fill(-1);
  double h = x[1] - x[0];

  for (int i = 0; i <= n-2; i++)
  {
    g[i] = 100*exp(-10*x[i+1])*h*h;
  }

  v = general_algorithm(v, g, a, b, c);
  for (int i = 1; i <= n; i++)
  {
    v_full[i] = v[i-1];
  }
  v_full[0] = 0;
  v_full[n] = 0;

  arma::mat data = arma::mat(n+1, 2);
  data.col(0) = x;
  data.col(1) = v_full;
  data.save("n10000.dat");

  return 0;
}

arma::vec general_algorithm(arma::vec v, arma::vec g, arma::vec a, arma::vec b, arma::vec c)
{
  int n_matrix = a.size();
  arma::vec c_tilde = arma::vec(n_matrix);
  arma::vec g_tilde = arma::vec(n_matrix);
  double w;
  c_tilde[0] = c[0]/b[0];
  g_tilde[0] = g[0]/b[0];

  for (int i = 1; i <= n_matrix-1; i++)
  {
    w = b[i] - a[i]*c_tilde[i-1];
    c_tilde[i] = c[i]/w;
    g_tilde[i] = (g[i] - a[i]*g_tilde[i-1]) / w;
  }


  v[n_matrix-1] = g_tilde[n_matrix-1];

  for (int i = n_matrix-2; i >= 0; i--)
  {
    v[i] = g_tilde[i]- c_tilde[i]*v[i+1];
  }

return v;
}
