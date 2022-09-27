#include "functions.hpp"
// Compile and linking with:
// g++ problem_3.cpp functions.cpp -o problem3 -larmadillo

// test func using A given in prob 3 b).
int main()
{
  int k = 0; int l = 1;
  arma::mat A = {{1, 0, 0, 0.5},
                {0, 1, -0.7, 0},
                {0, -0.7, 1, 0},
                {0.5, 0, 0, 1}};


  double max = max_offdiag_symmetric(A, k, l);

  A.print();
  std::cout << arma::endl;
  std::cout << "maxval: "<< max << arma::endl;
  return 0;
}
