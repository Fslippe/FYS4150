#include "double_slit_box.hpp"

// Constructor
DoubleSlitBox::DoubleSlitBox()
{
    
}

// Tranlates a pair of indices (i,j) from an u^n (M-2)x(M-2) matrix to a corresponding single index k in a u^n vector.
int DoubleSlitBox::translate_indices(int i, int j) 
{
    return (i-1) + (M-2)*(j-1);
}

// Fills mnatrices A and B for CN scheme 
void DoubleSlitBox::fill_A_B( int M, double h, double dt, arma::cx_mat V)
{
    int n = (M-2);
    int N = n*n; // matix size NxN

//-------------------------------------------- VECTOR FILLER
    // complex values
    std::complex<double> r(0, (dt / (2.*h*h)));
    std::complex<double> a_term = 1. + 4.*r;
    std::complex<double> b_term = 1. - 4.*r;
    std::complex<double> ab_term(0, dt / 2.);
    
    // create a and b vectors
    arma::cx_vec a(N);
    arma::cx_vec b(N);

    // fill a and b vectors
    for (int i = 1; i <= n; i++ )
    {
         for (int j = 1; i <= n; i++ )
         {
            a( translate_indices( i, j) ) = a_term + ab_term * V(i,j);
            b( translate_indices( i, j) ) = b_term - ab_term * V(i,j);
         }
    }

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

}

// Evolves the system one time step (dt) using the Crank-Nicholson scheme
void DoubleSlitBox::evolve_CN(double dt)
{
    b = B * u;
    u = arma::spsolve(A,b);
}

// Sets up initial state based on an unnormalised Gaussian wave packet epression
void DoubleSlitBox::init_box()
{

}

// initializes the potential V
void DoubleSlitBox::init_V()
{

}