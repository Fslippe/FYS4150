#include <iostream>
#include <armadillo>
#include <string>
#include <fstream>
#include <iomanip>

int n = 101;
// Defining a vector of x values.
arma::vec x_vec = arma::linspace(0,1,n); //Declare and will with random uniform values.
arma::vec u_vec = arma::vec(n); //Initialize vector but don't fill.

double u_func(double x); // Declaration of u(x).

// Filling u_vec with u(x) values for x in x_vec.
int main() {
  for (int i = 0; i < n-1; i++)
  {
    u_vec(i) = u_func(x_vec(i));
  }
//            -- Output to a datafile --
// Set a filename
  std::string filename = "x_u.txt";
  // Create and open the output file. Or, technically, create
  // an "output file stream" (type std::ofstream) and connect it to our filename.
  std::ofstream ofile;
  ofile.open(filename);
  // Some width and precision parameters we will use to format the output
  int width = 12;
  int prec  = 4;
  // Loop over steps
  for (int i = 0; i <= n-1; i++)
  {
    // Write a line with the current x and u values (nicely formatted) to file
    ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x_vec(i)
          << std::setw(width) << std::setprecision(prec) << std::scientific << u_vec(i)
          << std::endl;
  }
  // Close the output file
  ofile.close();

  return 0;
}

// u(x)
double u_func(double x){
    return 1 - (1 - exp(-10))*x - exp(-10*x);
}
