#include <assert.h>
#include <armadillo>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

double max_offdiag_symmetric(const arma::mat A, int k, int l);

// test func using A given in prob 3 b).
int main()
{
 int k = 0; int l = 1;
 arma::mat A = {{1, 0, 0, 0.5},
                {0, 1, -0.7, 0},
                {0, -0.7, 1, 0},
                {0.5, 0, 0, 1}};

 int max = max_offdiag_symmetric(A, k, l);

  A.print();
 std::cout << arma::endl;
 std::cout << "maxval: "<< max << arma::endl;
    return 0;
}

// Finding the max off-diagonal element of a matrix
double max_offdiag_symmetric(const arma::mat A, int k, int l)
{
 int N = A.n_rows;
 assert(N > 1); 
 assert(A.is_square());

 double maxval = A(k,l); //(i,j) (row, col)
 int col_n = 1;

 for (int i = 0; i <= N-2; ++i){
  for (int j = col_n; j <=N-1; ++j){
    std::cout << A(i,j) << arma::endl;
    if(std::abs(A(i,j)) > std::abs(maxval)){
        k = i; l = j;
    }
    col_n += 1;
   }
  }
  return maxval;
 }
