#include "double_slit_box.hpp"

// Constructor
DoubleSlitBox::DoubleSlitBox(double h_in, double dt_in, double T_in, double xc_in, double sigmax_in, double px_in, double yc_in, double sigmay_in, double py_in, int v0_in)
{
    //Member variables
    h = h_in;
    dt = dt_in;
    T = T_in;
    xc = xc_in;
    sigmax = sigmax_in;
    px = px_in;
    yc = yc_in;
    sigmay = sigmay_in;
    py = py_in;
    v0 = v0_in;

    //Sizes
    M = 1/h; // total matrix size MxM
    n = (M-2); // inner box nxn
    N = n*n; // CN matix size NxN, vector size N

    //Double slix box matrices (total matrix size)
    //arma::cx_mat U(M,M);
    arma::cx_mat V(M,M);

    //Crank Nicholson matrices and vectors (internal matrix size)
    arma::cx_vec b(N);
    arma::cx_vec u(N);
    arma::sp_cx_mat A(N,N);
    arma::sp_cx_mat B(N,N);

}

// Tranlates a pair of indices (i,j) from an u^n (M-2)x(M-2) matrix to a corresponding single index k in a u^n vector.
int DoubleSlitBox::translate_indices(int i, int j) 
{
    return (i-1) + (M-2)*(j-1);
}

// Fills mnatrices A and B for CN scheme 
// NB! needs to return A B or A B are created in constructor
void DoubleSlitBox::fill_A_B( int M, double h, double dt, arma::cx_mat V)
{
    // int n = (M-2);
    // int N = n*n; // matix size NxN

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
            a( translate_indices( i, j) ) = a_term + ab_term * V( i, j);
            b( translate_indices( i, j) ) = b_term - ab_term * V( i, j);
         }
    }

    //--------------------------------------------- MATRIX FILLER
    // // Create A and B matrix
    // arma::sp_cx_mat A(N,N); // change to cx_mat for more readable print()
    // arma::sp_cx_mat B(N,N);

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
    u = arma::spsolve(A,b); // spslover is well suited since A is a sparse matrix
}

// Sets up initial state vector u0 based on an unnormalised Gaussian wave packet epression
void DoubleSlitBox::init_wave()
{
    // real and imaginary components
    double wp_Re;
    std::complex<double> wp_Im;
    // wave packet value at point (x,y)
    std::complex<double> wp_val;
    double x;
    double y;

    // loop over all (x,y) points
    for (int i = 1; i <= n; i++ )
    {
         for (int j = 1; i <= n; i++ )
         {
            x = i*h;
            y = i*h;
            
            // not sure if the below compelx nr esp works
            wp_Re = - (std::pow(x-xc,2)/ (2*std::pow(sigmax,2))) - (std::pow(y-yc,2)/ (2*std::pow(sigmay,2)));
            wp_Im = (0,px*(x-xc) + py*(y-yc));
            wp_val = std::exp(wp_Re + wp_Im);
            u( translate_indices( i, j) ) = wp_val;
         }
    }
    // NB! u needs to be normalized !!!!
   
}

// initializes the potential V matrix
void DoubleSlitBox::init_V()
{

}