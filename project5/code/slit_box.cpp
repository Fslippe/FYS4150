#include "slit_box.hpp"

// Constructor
SlitBox::SlitBox(double h_in, double dt_in, double T_in, double xc_in, double sigma_x_in, double px_in, double yc_in, double sigma_y_in, double py_in, int v0_in)
{
    // Member variables
    h = h_in;
    dt = dt_in;
    T = T_in;
    xc = xc_in;
    sigma_x = sigma_x_in;
    px = px_in;
    yc = yc_in;
    sigma_y = sigma_y_in;
    py = py_in;
    v0 = v0_in;

    // Sizes
    M = 1 / h + 1; // total matrix size MxM
    n = (M - 2);   // inner box nxn
    N = n * n;     // CN matix size NxN, vector size N

    // Double slix box matrices (total matrix size)
    V = arma::zeros(M, M);
    // Crank Nicholson matrices and vectors (internal matrix size)
    b = arma::cx_vec(N);
    u = arma::cx_vec(N);
    A = arma::sp_cx_mat(N, N);
    B = arma::sp_cx_mat(N, N);
}

// Tranlates a pair of indices (i,j) from an u^n (M-2)x(M-2) matrix to a corresponding single index k in a u^n vector.
int SlitBox::translate_indices(int i, int j)
{
    return (i - 1) + (M - 2) * (j - 1);
}

// Fills matrices A and B for CN scheme
void SlitBox::fill_A_B()
{
    //-------------------------------------------- VECTOR FILLER
    // complex values
    std::complex<double> r(0, (dt / (2. * h * h)));
    std::complex<double> a_term = 1. + 4. * r;
    std::complex<double> b_term = 1. - 4. * r;
    std::complex<double> ab_term(0, dt / 2.);

    // create a and b vectors
    arma::cx_vec a(N);
    arma::cx_vec b(N);

    // fill a and b vectors
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            a(translate_indices(i, j)) = a_term + ab_term * V(i, j);
            b(translate_indices(i, j)) = b_term - ab_term * V(i, j);
        }
    }

    //--------------------------------------------- MATRIX FILLER
    //  Create A and B matrix
    // fill diagonal and +-(M-2) super- and subdiagonal
    A.diag() = a;
    A.diag(n).fill(-r);
    A.diag(-n).fill(-r);

    B.diag() = b;
    B.diag(n).fill(r);
    B.diag(-n).fill(r);

    // fill first super and sub diagonals
    A(0, 1) = -r;
    A(1, 0) = -r;
    B(0, 1) = r;
    B(1, 0) = r;

    for (int i = 0; i < N; i++)
    {
        if ((i + 1) % n != 0)
        {
            A.diag(1)(i) = -r;
            A.diag(-1)(i) = -r;
            B.diag(1)(i) = r;
            B.diag(-1)(i) = r;
        }
    }
}

// Evolves the system one time step (dt) using the Crank-Nicholson scheme
void SlitBox::evolve_CN()
{
    b = B * u;
    u = arma::spsolve(A, b); // spslover is well suited since A is a sparse matrix
}

// Use evolve_CN a to evolve system to the chosen time T_in
void SlitBox::evolve_CN_to_time()
{
    // Initialize time vector to save
    int n_T = T / dt;
    arma::vec t = arma::linspace(0, T, n_T + 1);
    arma::mat t_save = t;
    t_save.save("data/t8.dat");

    // Initialize U cube to save
    u_save = arma::cube(n_T + 1, M, M);
    u_real = arma::cube(n_T + 1, M, M);
    u_imag = arma::cube(n_T + 1, M, M);

    arma::vec p_val(u.size());
    std::cout << "number of timesteps: " << n_T << "\n";
    for (int i = 0; i <= n_T; i++)
    {
        std::cout << i << "\n";
        evolve_CN();
        p_val = arma::vec(arma::real(u % arma::conj(u)));

        // Fill U matrices
        for (int k = 1; k < M - 1; k++)
        {
            for (int j = 1; j < M - 1; j++)
            {
                u_save(i, j, k) = p_val(translate_indices(j, k));
                u_real(i, j, k) = arma::vec(arma::real(u))(translate_indices(j, k));
                u_imag(i, j, k) = arma::vec(arma::imag(u))(translate_indices(j, k));
            }
        }
    }
}

// Sets up initial state vector u0 based on an unnormalised Gaussian wave packet epression
void SlitBox::init_wave()
{
    // real and imaginary wave packet components
    double wp_Re;
    double wp_Im;
    // wave packet value at point (x,y)
    std::complex<double> wp_val;
    double x;
    double y;
    double sum = 0;
    // loop over all (x,y) points
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            x = i * h;
            y = j * h;
            wp_Re = -(std::pow(x - xc, 2) / (2 * std::pow(sigma_x, 2))) - (std::pow(y - yc, 2) / (2 * std::pow(sigma_y, 2)));
            wp_Im = px * (x - xc) + py * (y - yc);
            wp_val = std::exp(std::complex<double>(wp_Re, wp_Im));
            sum += std::real(wp_val * std::conj(wp_val));
            u(translate_indices(i, j)) = wp_val;
        }
    }
    u /= std::sqrt(sum); // Normalizing sum of u* x u =1
}

// initializes the potential V matrix for any given number of slits
void SlitBox::init_V(int slits)
{
    double width = 0.02;
    double center_x = 0.5;
    double center_y = 0.5;
    double distance = 0.05;
    double aperture = 0.05;
    int odd = slits % 2;
    std::cout << odd;

    // runs over all indices
    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < M; i++)
        {
            // only fill for middle x-values
            if (i * h >= center_x - width / 2 && i * h <= center_x + width / 2)
            {
                // open middle for odd number of slits
                if (odd == 1 && j * h < center_y + aperture / 2 && j * h > center_y - aperture / 2)
                {
                    V(i, j) = 0;
                }
                else
                {
                    V(i, j) = v0;
                }
                // Fill out rest of slits
                for (int k = 0; k < (slits - odd) / 2; k++)
                {
                    if (
                        // upper slits
                        j * h >= (k + 0.5 * odd) * aperture + (k + 1 + 0.5 * (odd - 1)) * distance + center_y           // lower boundary
                            && j * h < (k + 1 + 0.5 * odd) * aperture + (k + 1 + 0.5 * (odd - 1)) * distance + center_y // upper boundary
                        // Lower slits
                        || (j * h < -(k + 0.5 * odd) * aperture - (k + 1 + 0.5 * (odd - 1)) * distance + center_y           // upper boundary
                            && j * h >= -(k + 1 + 0.5 * odd) * aperture - (k + 1 + 0.5 * (odd - 1)) * distance + center_y)) // lowe boundary
                    {
                        V(i, j) = 0; // open up the wall
                    }
                }
            }
        }
    }
    // save file
    V.save(std::string("data/V_") + std::to_string(slits) + std::string(".dat"));
}