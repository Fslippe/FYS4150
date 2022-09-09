#include <iostream>
#include <armadillo>
#include <string>
#include <fstream>
#include <iomanip>
#include <chrono>

arma::vec special_algorithm(arma::vec v, arma::vec g); // Declaration of u(x).

int main() {
  int n = 10; // Number of steps
  int ti = 1000;

  arma::vec x = arma::linspace(0, 1, n+1); //Declare and will with random uniform values.
  arma::vec g = arma::vec(n-1);
  arma::vec v = arma::vec(n-1);
  arma::vec v_full = arma::vec(n+1);
  double h = x[1] - x[0];

  for (int i = 0; i <= n-2; i++)
  {
    g[i] = 100*exp(-10*x[i+1])*h*h;
  }

  arma::vec time = arma::vec(ti);
  for (int i = 0; i <= ti-1; i++)
  {
    clock_t t1 = clock();
    v = special_algorithm(v, g);
    clock_t t2 = clock();
    double duration_seconds = ((double) (t2 - t1)) / CLOCKS_PER_SEC;
    time[i] = duration_seconds;
  }
  time.save("special_algorithm_time10.dat");

  for (int i = 1; i <= n; i++)
  {
    v_full[i] = v[i-1];
  }
  v_full[0] = 0;
  v_full[n] = 0;

  arma::mat data = arma::mat(n+1, 2);
  data.col(0) = x;
  data.col(1) = v_full;
  data.save("special_n100.dat");

  return 0;
}

//General algorithm to solve the matrix equation Au = g for a tridiagonal n x n matrix A and known g
arma::vec special_algorithm(arma::vec v, arma::vec g)
{
  int n_matrix = g.size();
  arma::vec b = arma::vec(n_matrix);
  double w;
  b[0] = 2;
  for (int i = 1; i <= n_matrix-1; i++)
  {
    w = -1/b[i-1];
    b[i] = 2 + w;
    g[i] = g[i] - w * g[i-1];
  }

  v[n_matrix-1] = g[n_matrix-1] / b[n_matrix-1];

  for (int i = n_matrix-2; i >= 0; i--)
  {
    v[i] = (g[i] + v[i+1]) / b[i];
  }
return v;
}
