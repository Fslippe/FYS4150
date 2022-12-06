#ifndef __double_slit_box__
#define __double_slit_box__

#include <string>
#include <armadillo>
//#include <time.h>
#include <cmath>
#include <complex>

class DoubleSlitBox
{
public:
    // Member variables
    double h;
    double dt;
    double T;
    double xc;
    double sigma_x;
    double px;
    double yc;
    double sigma_y;
    double py;
    int v0;

    // Sizes
    int M; // total matrix size MxM
    int n; // inner box nxn
    int N; // CN matix size NxN, , vector size N

    // Double slix box matrices (total matrix size)
    // arma::cx_mat U;
    arma::mat V;

    // Crank Nicholson matrices and vectors (internal matrix size)
    arma::cx_vec b;
    arma::cx_vec u;
    arma::sp_cx_mat A;
    arma::sp_cx_mat B;

    // Consturctor
    DoubleSlitBox(double h_in, double dt_in, double T_in, double xc_in, double simgax_in, double px_in, double yx_in, double sigmay_in, double py_in, int v0_in);

    // Tranlates a pair of indices (i,j) into a corresponding single index k
    int translate_indices(int i, int j);

    // Fills mnatrices A and B fir CN scheme
    // Input: enire state size M, spatial step size h, time step size dt and potential matrix V
    void fill_A_B();

    // Evolves the system one time step (dt) using the Crank-Nicholson scheme
    void evolve_CN();

    // Sets up initial state vector u0 based on an unnormalised Gaussian wave packet epression
    void init_wave();

    // initializes the potential V matrix
    void init_V();
    void evolve_CN_to_time();
};

#endif