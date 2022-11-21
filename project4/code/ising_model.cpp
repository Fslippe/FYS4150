#include "ising_model.hpp"

IsingModel::IsingModel(int dim_in)
{

    dim = dim_in;
    n_spins = dim * dim;

    lattice = arma::mat(dim, dim);
    w = arma::zeros(17);
    val_vec = arma::vec(4);
}

void IsingModel::init_lattice(double T_in, int seed, bool spin_order)
{
    sum = arma::zeros(5);
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

    for (int j = 0; j < dim; j++)
    {
        for (int i = 0; i < dim; i++)
        {
            E -= lattice(i, j) * (lattice(periodic(i, dim, -1), j) + lattice(i, periodic(j, dim, -1)));
        }
    }

    for (int i = 0; i < w.size(); i += 4)
    {
        w[i] = std::exp(-(i - 8) / T);
    }
}

int IsingModel::periodic(int idx, int lim, int offset)
{
    return (idx + lim + offset) % lim;
}

// Calculate energy at the given index using periodic boundary conditions
//
int IsingModel::energy(int ix, int iy)
{
    up = lattice(ix, periodic(iy, dim, 1));
    down = lattice(ix, periodic(iy, dim, -1));
    left = lattice(periodic(ix, dim, -1), iy);
    right = lattice(periodic(ix, dim, 1), iy);
    return 2 * lattice(ix, iy) * (up + down + left + right);
}

// Metropolis algorithm to sweep through the lattice of size N=dim^2
void IsingModel::metropolis()
{
    for (int y = 0; y < dim; y++)
    {
        for (int x = 0; x < dim; x++)
        {
            ix = rnd(generator) * (double)dim;
            iy = rnd(generator) * (double)dim;
            dE = energy(ix, iy);

            if (rnd(generator) <= w[dE + 8])
            {
                lattice(ix, iy) *= -1;
                M += 2 * lattice(ix, iy);
                E += dE;
            }
        }
    }
}

void IsingModel::MC_sample(int cycles, bool histogram)
{
    if (histogram)
    {
        // Instead of sum - vector containg E at each cycle
        histogram_values = arma::vec(cycles);
        for (int i = 0; i < cycles; i++)
        {
            metropolis();
            histogram_values(i) = E;
        }
        histogram_values /= n_spins;
    }
    else
    {
        // cumulative sum over the cycles
        for (int i = 0; i < cycles; i++)
        {
            metropolis();
            sum(0) += E;
            sum(1) += E * E;
            sum(2) += M;
            sum(3) += M * M;
            sum(4) += std::fabs(M);
        }
        output(cycles);
    }
}

void IsingModel::output(int cycles)
{
    // normalizing the sum
    double norm = 1 / ((double)cycles);
    arma::vec average_norm = sum * norm;

    double E_avg = average_norm(0);
    double E2_avg = average_norm(1);
    double M_avg = average_norm(2);
    double M2_avg = average_norm(3);
    double M_abs = average_norm(4);

    // Calculating values of interest
    e_avg = E_avg / n_spins;
    m_abs = M_abs / n_spins;
    Cv_avg = (E2_avg - E_avg * E_avg) / (n_spins * T * T);
    chi_avg = (M2_avg - M_abs * M_abs) / (n_spins * T);

    // Vector containing values of interest
    val_vec[0] = e_avg;
    val_vec[1] = m_abs;
    val_vec[2] = Cv_avg;
    val_vec[3] = chi_avg;
}
