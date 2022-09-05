#include <iostream>
#include <armadillo>
#include <string>
#include <fstream>
#include <iomanip>


arma::vec general_algorithm(arma::vec x); // Declaration of u(x).

int n = 101;
arma::vec x_vec = arma::linspace(0,1,n); //Declare and will with random uniform values.
arma::vec v_vec = arma::vec(n);

int main() {
  arma::vec v_vec = general_algorithm(x_vec);
//            -- Output to a datafile --
// filename
  std::string filename = "x_v.txt";

  std::ofstream ofile;
  ofile.open(filename);
  int width = 12;
  int prec  = 4;
  // Loop
  for (int i = 0; i <= n+1; i++)
  {
    ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x_vec(i)
          << std::setw(width) << std::setprecision(prec) << std::scientific << v_vec(i)
          << std::endl;
  }
  ofile.close();

  return 0;
}

arma::vec general_algorithm(arma::vec x)
{
  // Defining a vector a ,b and c
  double h = x[1] - x[0];
  arma::vec a_vec = arma::vec(n-2).fill(-1.);
  arma::vec b_vec = arma::vec(n-2).fill(2.);
  arma::vec c_vec = arma::vec(n-2).fill(-1.);
  arma::vec g_vec = arma::vec(n-2);
  arma::vec b_tilde_vec = arma::vec(n);
  arma::vec g_tilde_vec = arma::vec(n);

  b_tilde_vec[0] = b_vec[0];
  g_vec[0] = 100.;
  g_tilde_vec[0] = g_vec[0];
  for (int i = 1; i <= n-1; i++)
  {
    b_tilde_vec[i] = b_vec[i] - a_vec[i] / b_vec[i-1] * c_vec[i-1];
    g_vec[i] = 100 * exp(-10*x[i]) * h*h;
    g_tilde_vec[i] = g_vec[i] - a_vec[i] / b_vec[i-1] * g_tilde_vec[i-1];
  }
  v_vec[n+1] = 0.; //final boudary
  v_vec[0] = 0.; // initial boundary
  v_vec[n] = g_tilde_vec[n-1] / b_vec[n-1];
  for (int i = n-1; i >= 1; i--)
  {
    v_vec[i] = (g_tilde_vec[i] - c_vec[i] * v_vec[i+1]) / b_tilde_vec[i];
  }
return v_vec;
}
