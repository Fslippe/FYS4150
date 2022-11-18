// g++ main.cpp ising_model.cpp -o -larmadillo

#include "ising_model.hpp"
#include "omp.h"
#include "time.h"

int main(int argc, char *argv[])
{
    int threads = atoi(argv[1]);
    double T = atof(argv[2]);
    int lattice_dim = atoi(argv[3]);
    std::string order(argv[4]);
    if (order == "false")
    {
        bool order = false;
    }
    else
    {
        bool order = true;
    }
    std::string save1(argv[5]);
    std::string save2(argv[6]);

    int seed = 2;
    int cycles = 1000000;

    int n_T = 40;
    int n_cycles = 40;
    arma::vec cycle_array = arma::logspace(1, 6, n_cycles);

    arma::mat cycle_val = arma::mat(5, n_cycles);
    cycle_array.print();
    clock_t t1 = clock();

#pragma omp parallel for num_threads(threads)
    for (int j = 0; j < n_cycles; j++)
    {
        IsingModel IM = IsingModel(lattice_dim);
        cycle_val(0, j) = cycle_array[j];
        IM.init_lattice(T, seed, false);
        IM.MC_sample(cycle_array[j]);
        for (int i = 1; i < 5; i++)
        {
            cycle_val(i, j) = IM.val_vec(i - 1);
        }
    }
    clock_t t2 = clock();
    double duration_seconds = ((double)(t2 - t1)) / CLOCKS_PER_SEC / threads;
    std::cout << "Cycles loop time used: " << duration_seconds << "s\n";
    std::cout << "\nTEMPERATURE DIFFERENCE\n";
    cycle_val.save(save1);

    arma::vec temp = arma::linspace(0.5, 4, n_T);
    arma::mat T_val = arma::mat(5, n_T);

    t1 = clock();
#pragma omp parallel for num_threads(threads)
    for (int j = 0; j < n_T; j++)
    {
        IsingModel IM = IsingModel(lattice_dim);
        T_val(0, j) = temp[j];
        IM.init_lattice(temp[j], seed, false);
        IM.MC_sample(cycles);
        for (int i = 1; i < 5; i++)
        {
            T_val(i, j) = IM.val_vec(i - 1);
        }
    }

    t2 = clock();
    T_val.save(save2);
    duration_seconds = ((double)(t2 - t1)) / CLOCKS_PER_SEC / threads;
    std::cout << "Temperature loop time used: " << duration_seconds << "s\n";
    return 0;
}