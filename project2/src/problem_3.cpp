#include <assert.h>
#include <armadillo>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

double max_offdiag_symmetric(const arma::mat& A, int& k, int &l);

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

// Finding the max off-diagonal element of a matrix
double max_offdiag_symmetric(const arma::mat& A, int& k, int &l)
{
  int N = A.n_rows;
  assert(N > 1);
  assert(A.is_square());

  double maxval = std::abs(A(k,l)); //(i,j) (row, col)
  //std::cout << maxval << "\n";

  int col_n = 1;

  for (int i = 1; i <= N-1; ++i) //Loops over rows
  {
    for (int j = 0; j < i; ++j) //loops over columns until the subdiagonal
    {
      std::cout << A(i,j) << std::endl;
      if(std::abs(A(i,j)) > maxval)
      {
        k = i; l = j;
        maxval = std::abs(A(i,j));
      }
   }
  }
  return maxval;
 }
