#ifndef __ising_model_hpp__
#define __ising_model_hpp__

#include <string>
#include <armadillo>
#include <time.h>
#include <cmath>
class IsingModel
{
private:
    double up;
    double down;
    double left;
    double right;
    int ix;
    int iy;
    double rand_n;
    int dE;
    double T;
    int n_spins;
    int temp_idx;
    arma::vec average;

public:
    int N;
    int dim;
    arma::vec average_norm;
    arma::mat histogram_values;

    // arma::vec dE;
    arma::mat lattice;
    arma::vec E_diff;
    int M;
    int E;
    double e_avg;
    double m_abs;
    double Cv_avg;
    double chi_avg;
    arma::vec val_vec;
    arma::vec s;
    std::mt19937 generator;
    std::uniform_real_distribution<double> rnd;

    IsingModel(int dim_in);
    void init_lattice(double T, int seed, bool spin_order);
    void MC_sample(int cycles, bool histogram);
    void output(int cycles);
    int energy(int ix, int iy);
    void total_energy(int i);
    void total_magnetization(int i);
    int periodic(int idx, int lim, int offset);
    void metropolis();
};

#endif