#ifndef __ising_model_hpp__
#define __ising_model_hpp__

#include <string>
#include <armadillo>
#include <time.h>
#include <cmath>
class IsingModel
{
public:
    double up;                  // spin of index above
    double down;                // spin of index under
    double left;                // spin of index to the left
    double right;               // spin of index to the right
    int ix;                     // random x index generated
    int iy;                     // random y index generated
    double rand_n;              // random number generated
    int dE;                     // Energy of random index ix iy
    double T;                   // Temperature of system
    int n_spins;                // total number of spins in lattice
    arma::vec sum;              // sum for expectation values in MC sample
    int dim;                    // dimension of lattice
    arma::vec average_norm;     // normalized expectation values
    arma::vec histogram_values; // histogram values for each MC cycle

    arma::mat lattice;                          // lattice
    arma::vec w;                                // containing possible energies of system
    int M;                                      // magnetization
    int E;                                      // Energy
    double e_avg;                               // expected e
    double m_abs;                               // expected |m|
    double Cv_avg;                              // expected C_v
    double chi_avg;                             // expected chi
    arma::vec val_vec;                          // vector to contain 4 values above
    std::mt19937 generator;                     // random number generator
    std::uniform_real_distribution<double> rnd; // uniform distributuion

    IsingModel(int dim_in);                                 // constructor
    void init_lattice(double T, int seed, bool spin_order); // initialize lattice and values
    void MC_sample(int cycles, bool histogram);             // perform MC samples for any nymber of cycles
    void output(int cycles);                                // output run at end of MC_sample
    int energy(int ix, int iy);                             // Energy for any index ix, iy
    int periodic(int idx, int lim, int offset);             // periodic boundary condtions
    void metropolis();                                      // metropolis algorithm contained in MC_sample method
};

#endif