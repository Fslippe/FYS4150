
#ifndef __ising_model_hpp__
#define __ising_model_hpp__

#include <string>
#include <armadillo>
#include <time.h>


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

    public:
        int N;
        int dim;
        arma::vec M; 
        arma::vec E; 
        arma::vec dE;
        arma::mat lattice; 
        arma::vec E_diff;
        int M;

        arma::vec s; 
        std::mt19937 generator;
        std::uniform_real_distribution<double> rnd;

    IsingModel(int dim);
    int energy(int ix, int iy);
    void total_energy(int i);
    void total_magnetization(int i);
    void init_lattice(double T);
    int periodic(int idx, int lim, int offset);
    void metropolis();
};

#endif