// g++ main.cpp ising_model.cpp -o -larmadillo

#include "ising_model.hpp"

int main()
{
    int seed = 2;
    double T = 1;
    int cycles = 1000000;
    int lattice_dim = 2;

    IsingModel IM = IsingModel(lattice_dim);
    int n_T = 20;
    int n_cycles = 40;
    arma::vec cycle_array = arma::logspace(1, 6, n_cycles);

    arma::mat cycle_val = arma::mat(5, n_cycles);
    cycle_array.print();
    for (int j = 0; j < n_cycles; j++)
    {
        cycle_val(0, j) = cycle_array[j];
        IM.init_lattice(T, seed, false);
        IM.MC_sample(cycle_array[j]);
        for (int i = 1; i < 5; i++)
        {
            cycle_val(i, j) = IM.val_vec(i - 1);
        }
    }
    std::cout << "FERDIGE\n";
    // cycle_val.print();
    cycle_val.save("data/cycle_val.dat");
    std::cout << "NÅ\n";
    std::cout << "\nTEMPERATURE DIFFERENCE";
    arma::vec temp = arma::linspace(0.5, 4, n_T);
    arma::mat T_val = arma::mat(5, n_T);
    for (int j = 0; j < n_T; j++)
    {
        T_val(0, j) = temp[j];
        IM.init_lattice(temp[j], seed, false);
        IM.MC_sample(cycles);
        for (int i = 1; i < 5; i++)
        {

            T_val(i, j) = IM.val_vec(i - 1);
        }
    }
    T_val.save("data/T_val_1mill.dat");

    return 0;
}