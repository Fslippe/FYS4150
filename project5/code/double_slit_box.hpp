#ifndef __double_slit_box__
#define __double_slit_box__

#include <string>
#include <armadillo>
//#include <time.h>
#include <cmath>

class DoubleSlitBox
{
    private:



    public:
    // Consturctor
    DoubleSlitBox();

    // Tranlates a pair of indices (i,j) into a corresponding single index k
    int translate_indices(int, int); 

    // Fills mnatrices A and B
    fill_A_B( int M, double h, double dt, arma::mat V);
    
    // Evolves the system one time step (dt) using the Crank-Nicholson scheme
    void evolve_CN(double dt);
    
    // Sets up initial state based on an unnormalised Gaussian wave packet epression
    void init_box();

    // initializes the potential V
    void init_V();

};

#endif