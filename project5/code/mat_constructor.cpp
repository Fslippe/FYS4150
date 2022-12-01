#include <armadillo>
#include <cmath>
#include <complex>


int main()
{
    // input M, h, dt, V
    int M = 5;
    double dt = 1.5;
    double h = 1.6;
   //arma::cx_mat V;
//--------------------------------------------
    int n = (M-2);
    int N = n*n; // matix size NxN
    //int r = 3;

    // arma::vec P(10, arma::fill::randu);
    // arma::cx_vec a = roots(P);
    // arma::cx_vec b = roots(P);
//-------------------------------------------- VECTOR FILLER
    // complex values
    std::complex<double> r(0, (dt / (2.*h*h)));
    std::complex<double> a_term = 1. + 4.*r;
    std::cout << "\n \n";
    std::cout << a_term;
    std::complex<double> b_term = 1. - 4.*r;
    std::complex<double> ab_term(0, dt / 2.);
    
    // create a and b vectors
    arma::cx_vec a(N);
    arma::cx_vec b(N);

    // // fill a and b vectors
    // for (int i = 1; i <= n; i++ )
    // {
    //      for (int j = 1; i <= n; i++ )
    //      {
    //         a( translate_indices( i, j) ) = a_term + ab_term * V(i,j);
    //         b( translate_indices( i, j) ) = b_term - ab_term * V(i,j);
    //      }
    // }

//--------------------------------------------- MATRIX FILLER
    // Create A and B matrix
    arma::sp_cx_mat A(N,N); // change to cx_mat for more readable print()
    arma::sp_cx_mat B(N,N);

    //fill diagonal and +-(M-2) super- and subdiagonal
    A.diag() = a;
    A.diag(n).fill(-r);
    A.diag(-n).fill(-r);

    B.diag() = b;
    B.diag(n).fill(r);
    B.diag(-n).fill(r);

    //fill first super and sub diagonals
    A(0,1) = -r;
    A(1,0) = -r;
    B(0,1) = r;
    B(1,0) = r;

    for (int i = 0; i < N; i++)
    {
        if ( (i+1)%n !=0 )
        {
            A.diag(1)(i) = -r;
            A.diag(-1)(i) = -r;
            B.diag(1)(i) = r;
            B.diag(-1)(i) = r;
        }

    }

//---------------------------------------------
  

    // A.print();
    // std::cout << "\n \n";
    // B.print();
    
    return 0;
}