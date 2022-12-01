#include "double_slit_box.hpp"

    // Constructor
    DoubleSlitBox::DoubleSlitBox()
    {

    }

    // Tranlates a pair of indices (i,j) into a corresponding single index k
    int DoubleSlitBox::translate_indices(int i, int j) 
    {
        return (i-1) + (M-2)*(j-1);
    }

    // Fills mnatrices A and B
    DoubleSlitBox::fill_A_B( int M, double h, double dt, arma::mat V)
    {

    }

    // Evolves the system one time step (dt) using the Crank-Nicholson scheme
    void DoubleSlitBox::evolve_CN(double dt)
    {

    }

    // Sets up initial state based on an unnormalised Gaussian wave packet epression
    void DoubleSlitBox::init_box()
    {

    }

    // initializes the potential V
    void DoubleSlitBox::init_V()
    {

    }