#include "ising_model.hpp"

const double k_b = 1.38064903*std::pow(10, -23);

IsingModel::IsingModel(int dim_in)
{
    std::mt19937 generator(1583023);
    std::uniform_real_distribution<double> rnd(0.0, 1.0);   
    dim = dim_in;

    dE = 1;
    E = 0;
    M = 0;

    lattice = arma::mat(dim, dim);
    E_diff = arma::vec(16);
    average = arma::zeros(5);
    average = arma::zeros(5);


}
int IsingModel::periodic(int idx, int lim, int offset)
{
    return (idx + lim + offset) % offset;
}
int IsingModel::test()
{
    std::cout << E << "\n";

 return M;   
}
void IsingModel::init_lattice(double T_in)
{   
    T = T_in;
    std::cout << E << "\n";
    // Initialize random spin and total magnetization
    for (int j = 0; j < dim; j++)
    {
        for (int i = 0; i < dim; i++)
        {
            rand_n = rnd(generator);
            lattice(i,j) = (rand_n < 0.5) ? 1 : -1;
            M += lattice(i,j);
        }

    }

    // Initialize Energy
    for (int j = 0; j < dim; j++)
    {
        for (int i = 0; i < dim; i++)
        {

            E -= lattice(i,j)*lattice(periodic(i, dim, -1), j) + lattice(i, periodic(j, dim, -1));
        }
    }   
    std::cout << dim << "\n";

    for (int i = -8; i <= 8; i++)
    {
        E_diff[i+8] = 0;
    }

    for (int i = -8; i <= 8; i+=4)
    {
        E_diff[i+8] = std::exp(-i / T);
    }
    
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
    //std::cout << "METROPOLIS" << dim << "\n";

    for(int y = 0; y < dim; y++) 
    {
        for (int x = 0; x < dim; x++)
        {
            ix = rnd(generator) * dim;
            iy = rnd(generator) * dim;
            dE = energy(ix, iy);
            if (rnd(generator) <= E_diff[dE + 8])
            {
                lattice(ix, iy) *= -1; 
                M +=  2* lattice(ix, iy);
                E +=  dE;

            }

        }
    }
}

void IsingModel::MC_sample(int cycles)
{
    std::cout  << "E  " << E << "\n"; 
    
    for (int i = 0; i < cycles; i++)
    {
        
        metropolis();
        average(0) += E;
        average(1) += E*E;
        average(2) += M;
        average(3) += M*M;
        average(4) += std::fabs(M);
    }
    output(cycles);
}
    

void IsingModel::output(int cycles)
{   
    double norm = 1 / ((double) cycles);
    arma::vec average_norm = average * norm; 
    int dim2 = dim*dim;
    double E_avg =  average_norm(0) / dim2;
    double E2_avg = average_norm(1) / dim2;
    double M_avg = average_norm(2) / dim2;
    double M2_avg = average_norm(3) / dim2;
    double M_abs = average_norm(4) / dim2;
    double Cv_avg = (E2_avg - E_avg*E_avg) /(dim2 *  T*T);
    double chi_avg = (M2_avg - M_abs*M_abs) /(dim2 *  T);
    arma::vec average = arma::vec(7);
    average[0] = E_avg;
    average[1] = M_avg;
    average[2] = E2_avg;
    average[3] = M2_avg;
    average[4] = M_abs;
    average[5] = Cv_avg;
    average[6] = chi_avg;
    //average.print();
    average.save("data/average.dat");
}
