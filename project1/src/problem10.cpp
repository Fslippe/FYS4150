#include <time.h>
#include <iostream>
#include <armadillo>
#include <string>
#include <fstream>
#include <iomanip>
// include special and general alrgorithm (create hpp file with algorithm class?)


arma::vec general_algorithm(arma::vec v, arma::vec g, arma::vec a, arma::vec b, arma::vec c); // Declaration
arma::vec special_algorithm(arma::vec v, arma::vec g); // Declaration

int main ()
{

  int n = 100 // Numer of steps. up to 10^6
  int runs = 1000 // Number of runs

// loop for repeated runs of each choice of n
double g_total = 0
double s_total = 0
for (int j = 1; j <= runs; j++)
  {

  clock_t t1_g = clock();  // Start

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

  clock_t t2_g = clock();  // Stop

  clock_t t1_s = clock(); // Start

  
  // special_algorithm here


  clock_t t2_s = clock();   // Stop

  double duration_seconds_g = ((double) (t2_g - t1_g)) / CLOCKS_PER_SEC; // time for general algorithm
  double duration_seconds_s = ((double) (t2_s - t1_s)) / CLOCKS_PER_SEC; // time for speciel algorithm

  g_total += duration_seconds_g
  s_total += duration_seconds_s
  }
// average out time mesurments
double g_final = g_total / runs
double s_final = s_total / runs
  // write to file
  arma::mat data = arma::mat(n+1, 2);
  data.col(0) = g_final;
  data.col(1) = s_final;
  data.save("n100_time.dat");
}
