// g++ main.cpp ising_model.cpp -o -larmadillo

#include "ising_model.hpp"
#include "omp.h"
#include "time.h"

void temp_loop(int n_T, int cycles, int seed, double T, int lattice_dim, bool order, std::string save1, int threads, double T_min, double T_max);
void cycle_loop(int n_cycles, int seed, double T, int lattice_dim, bool order, std::string save2, int threads);

int main(int argc, char *argv[])
{
    // Commandline arguments
    int threads = atoi(argv[1]);
    double T = atof(argv[2]);
    int lattice_dim = atoi(argv[3]);
    std::string order_in(argv[4]);
    std::string save1(argv[5]);
    std::string save2(argv[6]);
    std::string save3(argv[7]);
    double T_min = atof(argv[8]);
    double T_max = atof(argv[9]);

    bool order;
    if (order_in == "false")
    {
        order = false;
    }
    else
    {
        order = true;
    }

    // Parameters
    int seed = 2;
    int cycles = 1000000;
    int n_T = 40;
    int n_cycles = 40;
    double duration_seconds;

    // Values depending on number of cycles
    clock_t t1;
    clock_t t2;
    if (save1 != "none")
    {
        std::cout << "\nCYCLE DIFFERENCE\n";
        t1 = clock();
        cycle_loop(n_cycles, seed, T, lattice_dim, order, save1, threads);
        t2 = clock();
        duration_seconds = ((double)(t2 - t1)) / CLOCKS_PER_SEC / threads;
        std::cout << "Cycles loop time used: " << duration_seconds << "s\n";
    }

    // Values depending on temperature
    if (save2 != "none")
    {
        std::cout << "\nTEMPERATURE DIFFERENCE\n";
        t1 = clock();
        temp_loop(n_T, cycles, seed, T, lattice_dim, order, save2, threads, T_min, T_max);
        t2 = clock();
        duration_seconds = ((double)(t2 - t1)) / CLOCKS_PER_SEC / threads;
        std::cout << "Temperature loop time used: " << duration_seconds << "s\n";
    }

    // Histogram
    if (save3 != "none")
    {
        int cycles = 10000000;

        arma::mat hist = arma::mat(cycles, 1);
        std::cout << "\nHISTOGRAM\n";
        t1 = clock();
        IsingModel IM = IsingModel(lattice_dim);
        IM.init_lattice(T, seed, order);
        IM.MC_sample(cycles, true);
        hist.col(0) = IM.histogram_values;
        hist.save(save3);
        t2 = clock();
        duration_seconds = ((double)(t2 - t1)) / CLOCKS_PER_SEC;
        std::cout << "Temperature loop time used: " << duration_seconds << "s\n";
    }

    return 0;
}

void temp_loop(int n_T, int cycles, int seed, double T, int lattice_dim, bool order, std::string save1, int threads, double T_min, double T_max)
{
    arma::mat T_val;
    arma::vec temp;
    if (T_min == 2.1)
    {
        arma::vec a = arma::linspace(T_min, 2.23, 15);
        arma::vec b = arma::linspace(2.235, 2.33, 30);
        arma::vec c = arma::linspace(2.335, T_max, 10);
        temp = arma::join_cols(a, b, c);
        // change
        n_T = 55;
        T_val = arma::mat(5, n_T);
    }
    else
    {
        temp = arma::linspace(T_min, T_max, n_T); // change
        T_val = arma::mat(5, n_T);
    }

#pragma omp parallel for num_threads(threads)
    for (int j = 0; j < n_T; j++)
    {
        IsingModel IM = IsingModel(lattice_dim);
        T_val(0, j) = temp[j];
        IM.init_lattice(temp[j], seed, order);
        IM.MC_sample(cycles, false);
        for (int i = 1; i < 5; i++)
        {
            T_val(i, j) = IM.val_vec(i - 1);
        }
    }
    T_val.save(save1);
}
void cycle_loop(int n_cycles, int seed, double T, int lattice_dim, bool order, std::string save2, int threads)
{
    arma::vec cycle_array;
    if (lattice_dim == 2)
    {
        cycle_array = arma::logspace(1, 7, n_cycles);
    }
    else
    {
        cycle_array = arma::logspace(1, 6, n_cycles);
    }

    arma::mat cycle_val = arma::mat(5, n_cycles);

#pragma omp parallel for num_threads(threads)
    for (int j = 0; j < n_cycles; j++)
    {
        IsingModel IM = IsingModel(lattice_dim);
        cycle_val(0, j) = cycle_array[j];
        IM.init_lattice(T, seed, order);
        IM.MC_sample(cycle_array[j], false);
        for (int i = 1; i < 5; i++)
        {
            cycle_val(i, j) = IM.val_vec(i - 1);
        }
    }
    cycle_val.save(save2);
}
