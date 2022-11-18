#include "ising_model.hpp"
#include "omp.h"
IsingModel::IsingModel(int dim_in)
{

    dim = dim_in;
    n_spins = dim * dim;

    lattice = arma::mat(dim, dim);
    E_diff = arma::zeros(17);
    val_vec = arma::vec(4);
}

void IsingModel::init_lattice(double T_in, int seed, bool spin_order)
{
    average = arma::zeros(5);
    dE = 0;
    E = 0;
    M = 0;
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
    // std::cout << lattice << "\n";
    //  Initialize Energy
    for (int j = 0; j < dim; j++)
    {
        for (int i = 0; i < dim; i++)
        {
            E -= lattice(i, j) * (lattice(periodic(i, dim, -1), j) + lattice(i, periodic(j, dim, -1)));
            // std::cout << "i: " << i << " ";
            // std::cout << "j: " << j << " ";
            // std::cout << "E: " << E << "\n";
        }
    }

    for (int i = 0; i < E_diff.size(); i += 4)
    {
        E_diff[i] = std::exp(-(i - 8) / T);
    }
    // std::cout << "E after: " << E_diff << "\n";
}

int IsingModel::periodic(int idx, int lim, int offset)
{
    return (idx + lim + offset) % lim;
}

int IsingModel::energy(int ix, int iy)
{
    up = lattice(ix, periodic(iy, dim, 1));
    down = lattice(ix, periodic(iy, dim, -1));
    left = lattice(periodic(ix, dim, -1), iy);
    right = lattice(periodic(ix, dim, 1), iy);
    return 2 * lattice(ix, iy) * (up + down + left + right);
}

double IsingModel::metropolis(double E_loc, double M_loc)
{
    // E_diff.print();
    // lattice.print();

    for (int y = 0; y < dim; y++)
    {
        for (int x = 0; x < dim; x++)
        {
            ix = rnd(generator) * (double)dim;
            iy = rnd(generator) * (double)dim;
            // std::co  ut << "ix iy " << ix << " " << iy << "\n";
            dE = energy(ix, iy);

            // std::cout << "\ndE " << dE << "\n\n";

            if (rnd(generator) <= E_diff[dE + 8])
            {
                lattice(ix, iy) *= -1;
                M_loc += 2 * lattice(ix, iy);
                E_loc += dE;
            }
            // lattice.print();
            //   std::cout << "E:  " << E << "\n";
        }
    }
    return E_loc, M_loc;
}

void IsingModel::MC_sample(int cycles)
{

    average_global = arma::zeros(5);
    arma::vec average_loc;
#pragma omp parallel private(average_loc) shared(average_global)
    {
        double E_loc = E, M_loc = M;
        double E_test = 0;
        average_loc = arma::zeros(5);
        const int n_threads = omp_get_num_threads();
        const double weight = 1.0 / n_threads;
        const int cycle_threads = cycles * weight;

#pragma omp for
        for (int i = 0; i < cycle_threads; i++)
        {
            E_loc, M_loc = metropolis(E_loc, M_loc);
            // std::cout << E << "\n";
            average_loc(0) += E_loc;
            average_loc(1) += E_loc * E_loc;
            average_loc(2) += M_loc;
            average_loc(3) += M_loc * M_loc;
            average_loc(4) += std::fabs(M_loc);
        }
#pragma omp critical
        average_global += weight * average_loc;
    }
    output(cycles);
}

void IsingModel::output(int cycles)
{

    double norm = 1 / ((double)cycles);
    arma::vec average_norm = average_global * norm;

    double E_avg = average_norm(0);
    double E2_avg = average_norm(1);
    double M_avg = average_norm(2);
    double M2_avg = average_norm(3);
    double M_abs = average_norm(4);
    e_avg = E_avg / n_spins;
    m_abs = M_abs / n_spins;
    Cv_avg = (E2_avg - E_avg * E_avg) / (n_spins * T * T);
    chi_avg = (M2_avg - M_abs * M_abs) / (n_spins * T);

    val_vec[0] = e_avg;
    val_vec[1] = m_abs;
    val_vec[2] = Cv_avg;
    val_vec[3] = chi_avg;
    // val_vec.print();
    int J = 1;
    double beta = 1 / T;
    double Z = 4 * (std::cosh(8 * beta * J) + 3);
    double E_an = 16 / Z * (std::exp(-beta * 8) - std::exp(beta * 8)) / 4;
}
