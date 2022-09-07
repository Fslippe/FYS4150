#include <iostream>
#include <armadillo>
#include <string>
#include <fstream>
#include <iomanip>

arma::vec special_algorithm(arma::vec v, arma::vec g); // Declaration of u(x).

int main() {
  int n = 100; // Number of steps

  arma::vec x = arma::linspace(0, 1, n+1); //Declare and will with random uniform values.
  arma::vec g = arma::vec(n-1);
  arma::vec v = arma::vec(n-1);
  arma::vec v_full = arma::vec(n+1);
  double h = x[1] - x[0];

  for (int i = 0; i <= n-2; i++)
  {
    g[i] = 100*exp(-10*x[i+1])*h*h;
  }

  v = special_algorithm(v, g);

  for (int i = 1; i <= n; i++)
  {
    v_full[i] = v[i-1];
  }
  v_full[0] = 0;
  v_full[n] = 0;

  arma::mat data = arma::mat(n+1, 2);
  data.col(0) = x;
  data.col(1) = v_full;
  data.save("n100.dat");

  return 0;
}

//General algorithm to solve the matrix equation Au = g for a tridiagonal n x n matrix A and known g
arma::vec special_algorithm(arma::vec v, arma::vec g)
{
  int n_matrix = g.size();
  double a = -1; double b = 2; double c = -1;

  double w;

  arma::vec c_tilde = arma::vec(n_matrix);
  arma::vec g_tilde = arma::vec(n_matrix);


  w = a/b;
  b_tilde = b - w*c ;
  for (int i = 1; i <= n_matrix-1; i++)
  {

    g_tilde[i] = g[i] - w*g_tilde[i-1];
  }

  v[n_matrix-1] = g_tilde[n_matrix-1]/b_tilde;

  for (int i = n_matrix-2; i >= 0; i--)
  {
    v[i] = g_tilde/b - c/b*v[i+1];
  }
  std::string filename = "x_v.txt";

  std::ofstream ofile;
  ofile.open(filename);
  int width = 12;
  int prec  = 4;
  // Loop
  for (int i = 0; i <= n_matrix-1; i++)
  {
    ofile << std::setw(width) << std::setprecision(prec) << std::scientific << v(i)
          << std::setw(width) << std::setprecision(prec) << std::scientific << v(i)
          << std::endl;
  }
  ofile.close();
return v;
}
