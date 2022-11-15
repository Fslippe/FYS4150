#include "ising_model.hpp"

IsingModel::IsingModel(int dim)
{
    dim = dim;
    dE = 1;
    std::mt19937 generator(1583023);
    std::uniform_real_distribution<double> rnd(0.0, 1.0);   
    lattice = arma::mat(dim, dim);
    E_diff = arma::vec(16);

}
int IsingModel::periodic(int idx, int lim, int offset)
{
    return (idx + lim + offset) % offset;
}

void IsingModel::init_lattice(double T)
{
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

    for (int i = -8; i <= 8; i++)
    {
        E_diff[i+8] = 0;
    }

    for (int i = -8; i <= 8; i+=4)
    {
        E_diff[i+8] = std::exp(-i / T);
    }
    
}

void IsingModel::total_energy(int i)
{
    for (int j = 0; j < N; j++)
    {
        E[i] += s[j] * s[j+1];
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
    for(int y = 0; y < dim; y++) 
    {
        for (int x = 0; x < dim; x++)
        {
            ix = rnd(generator) * dim;
            iy = rnd(generator) * dim;
            int dE = energy(ix, iy);
            if (rnd(generator) <= E_diff[dE + 8])
            {
                lattice(ix, iy) *= -1; 
                M +=  2* lattice(ix, iy);
                E +=  dE;
            }
        }
    }
}


// function to initialise energy, spin matrix and magnetization
void initialize(int n_spins, double temp, int **spin_matrix,double& E, double& M)
{
    // setup spin matrix and intial magnetization
    for(int y =0; y < n_spins; y++)
    {
        for (int x= 0; x < n_spins; x++)
        {
            if (temp < 1.5) spin_matrix[y][x] = 1; // spin orientation for the ground state
            M += (double) spin_matrix[y][x];
        }
    }
        // setup initial energy
    for(int y =0; y < n_spins; y++) 
    {
        for (int x= 0; x < n_spins; x++)
        {
            E -= (double) spin_matrix[y][x]*
            (spin_matrix[periodic(y,n_spins,-1)][x] +
            spin_matrix[y][periodic(x,n_spins,-1)]);
        }
    }
}// end function initialise


// Function to read in data from screen
void read_input(int&, int&, double&, double&, double&);
// Function to initialise energy and magnetization
void initialize(int, double, int **, double&, double&);
// The metropolis algorithm
void Metropolis(int, long&, int **, double&, double&, double *);
// prints to file the results of the calculations
void output(int, int, double, double *);

// main program
int main(int argc, char* argv[])
{
    char *outfilename;
    long idum;
    int **spin_matrix, n_spins, mcs;
    double w[17], average[5], initial_temp, final_temp, E, M, temp_step;
    // Read in output file, abort if there are too few command-line arguments
    if( argc <= 1 )
    {
        cout << "Bad Usage: " << argv[0] <<
        " read also output file on same line" << endl;
        exit(1);
    }
    else
    {
        outfilename=argv[1];
    }
    ofile.open(outfilename);
    // Read in initial values such as size of lattice, temp and cycles
    read_input(n_spins, mcs, initial_temp, final_temp, temp_step);
    spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
    idum = -1; // random starting point
    for ( double temp = initial_temp; temp <= final_temp; temp+=temp_step)
    {
        // initialise energy and magnetization
        E = M = 0.;
        // setup array for possible energy changes
        for( int de =-8; de <= 8; de++) w[de+8] = 0;
        for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temp);

        // initialise array for expectation values
        for( int i = 0; i < 5; i++) average[i] = 0.;
        initialize(n_spins, double temp, spin_matrix, E, M);

        // start Monte Carlo computation
        for (int cycles = 1; cycles <= mcs; cycles++)
        {
            Metropolis(n_spins, idum, spin_matrix, E, M, w);
            // update expectation values
            average[0] += E; average[1] += E*E;
            average[2] += M; average[3] += M*M; average[4] += fabs(M);
        }
        // print results
        output(n_spins, mcs, temp, average);
    }
    free_matrix((void **) spin_matrix); // free memory
    ofile.close(); // close output file
    return 0;
}