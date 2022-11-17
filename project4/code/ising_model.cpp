#include "ising_model.hpp"

const double k_b = 1; // 1.38064903 * std::pow(10, -23);

IsingModel::IsingModel(int dim_in)
{

    dim = dim_in;
    n_spins = dim * dim;
    temp_idx = 0;
    dE = 1;
    E = 0;
    M = 0;

    lattice = arma::mat(dim, dim);
    E_diff = arma::zeros(17);
    average = arma::zeros(5);
}

void IsingModel::init_lattice(double T_in, int seed, bool spin_order)
{
    T = T_in;
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> rnd(0.0, 1.0);

    //  Initialize random spin and total magnetization
    if (spin_order)
    {
        lattice.ones();
        M = arma::accu(lattice);
    }
    else
    {
        for (int j = 0; j < dim; j++)
        {
            for (int i = 0; i < dim; i++)
            {
                rand_n = rnd(generator);
                lattice(i, j) = (rand_n < 0.5) ? 1 : -1;
                M += lattice(i, j);
            }
        }
    }
    std::cout << lattice << "\n";
    // Initialize Energy
    for (int j = 0; j < dim; j++)
    {
        for (int i = 0; i < dim; i++)
        {
            E -= lattice(i, j) * lattice(periodic(i, dim, -1), j) + lattice(i, periodic(j, dim, -1));
            // std::cout << "i: " << i << " ";
            // std::cout << "j: " << j << " ";
            // std::cout << "E: " << E << "\n";
        }
    }

    for (int i = 0; i < E_diff.size(); i += 4)
    {
        E_diff[i] = std::exp(-(i - 8) / T);
        // std::cout << "E after: " << E_diff << "\n";
    }
    E_diff.print();
}

int IsingModel::periodic(int idx, int lim, int offset)
{
    return (idx + lim + offset) % offset;
}

int IsingModel::energy(int ix, int iy)
{
    up = lattice(ix, periodic(iy, dim, 1));
    down = lattice(ix, periodic(iy, dim, -1));
    left = lattice(periodic(ix, dim, -1), iy);
    right = lattice(periodic(ix, dim, 1), iy);
    return 2 * lattice(ix, iy) * (up + down + left + right);
}

void IsingModel::metropolis()
{
    for (int y = 0; y < dim; y++)
    {
        for (int x = 0; x < dim; x++)
        {
            ix = rnd(generator) * dim;
            iy = rnd(generator) * dim;
            // std::cout << "ix iy " << ix << " " << iy << "\n";
            dE = energy(ix, iy);
            std::cout << "dE " << dE << "\n";
            if (rnd(generator) <= E_diff[dE + 8])
            {
                // std::cout << "FLIP\n";
                lattice(ix, iy) *= -1;
                M += 2 * lattice(ix, iy);
                E += dE;
            }
            std::cout << "E:  " << E << "\n";
        }
    }
}

void IsingModel::MC_sample(int cycles)
{

    for (int i = 0; i < cycles; i++)
    {
        metropolis();
        average(0) += E;
        average(1) += E * E;
        average(2) += M;
        average(3) += M * M;
        average(4) += std::fabs(M);
    }
    output(cycles);
}

void IsingModel::output(int cycles)
{
    temp_idx += 1;
    double norm = 1 / ((double)cycles);
    arma::vec average_norm = average * norm;

    double E_avg = average_norm(0) / n_spins;
    double E2_avg = average_norm(1) / n_spins;
    double M_avg = average_norm(2) / n_spins;
    double M2_avg = average_norm(3) / n_spins;
    double M_abs = average_norm(4) / n_spins;
    double Cv_avg = (E2_avg - E_avg * E_avg) / (n_spins * T * T);
    double chi_avg = (M2_avg - M_abs * M_abs) / (n_spins * T);
    arma::vec average = arma::vec(7);
    int J = 1;
    double beta = 1 / T; // k_b / T;
    double Z = 4 * (std::cosh(8 * beta * J) + 3);
    double E_an = 16 / Z * (std::exp(-beta * 8) - std::exp(beta * 8)) / 4;
    std::cout << "E_ANALYTIC" << E_an << "\n";

    average[0] = E_avg;
    average[1] = M_avg;
    average[2] = E2_avg;
    average[3] = M2_avg;
    average[4] = M_abs;
    average[5] = Cv_avg;
    average[6] = chi_avg;
    average.print();
    average.save("data/average.dat");
}
