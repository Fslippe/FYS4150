// g++ main.cpp ising_model.cpp -o -larmadillo

#include "ising_model.hpp"

int main()
{
    int seed = 1032;
    double T = 2.4;
    int cycles = 50;
    int lattice_dim = 2;
    IsingModel IM = IsingModel(lattice_dim);
    IM.init_lattice(T, seed, true);
    IM.MC_sample(cycles);
    return 0;
}