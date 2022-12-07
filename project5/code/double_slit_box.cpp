#include "double_slit_box.hpp"

// Constructor
DoubleSlitBox::DoubleSlitBox(double h_in, double dt_in, double T_in, double xc_in, double sigma_x_in, double px_in, double yc_in, double sigma_y_in, double py_in, int v0_in)
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
    M = 1 / h;   // total matrix size MxM
    n = (M - 2); // inner box nxn
    N = n * n;   // CN matix size NxN, vector size N

    // Double slix box matrices (total matrix size)
    // arma::cx_mat U(M,M);
    V = arma::mat(M, M);
    // Crank Nicholson matrices and vectors (internal matrix size)
    b = arma::cx_vec(N);
    u = arma::cx_vec(N);
    A = arma::sp_cx_mat(N, N);
    B = arma::sp_cx_mat(N, N);
}

// Tranlates a pair of indices (i,j) from an u^n (M-2)x(M-2) matrix to a corresponding single index k in a u^n vector.
int DoubleSlitBox::translate_indices(int i, int j)
{
    return (i - 1) + (M - 2) * (j - 1);
}

// Fills mnatrices A and B for CN scheme
// NB! needs to return A B or A B are created in constructor
void DoubleSlitBox::fill_A_B()
{
    // int n = (M-2);
    // int N = n*n; // matix size NxN

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
            // std::cout << V << "\n";
            // std::cout << i << j << translate_indices(i, j) << "\n";
            // std::cout << a_term << ab_term << "\n";

            a(translate_indices(i, j)) = a_term + ab_term * V(i, j);
            b(translate_indices(i, j)) = b_term - ab_term * V(i, j);
        }
    }

    //--------------------------------------------- MATRIX FILLER
    //  Create A and B matrix
    // arma::sp_cx_mat A(N,N); // change to cx_mat for more readable print()
    // arma::sp_cx_mat B(N,N);
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
void DoubleSlitBox::evolve_CN()
{
    b = B * u;

    u = arma::spsolve(A, b); // spslover is well suited since A is a sparse matrix
}

void DoubleSlitBox::evolve_CN_to_time()
{
    int n_T = T / dt;
    arma::vec t = arma::linspace(0, T, n_T + 1);
    arma::mat t_save = t;
    t_save.save("data/t8.dat");
    arma::cube u_save = arma::zeros(n_T + 1, M, M);
    arma::cube u_real = arma::zeros(n_T + 1, M, M);
    arma::cube u_imag = arma::zeros(n_T + 1, M, M);

    arma::vec p_val(u.size());
    std::cout << "number of timesteps: " << n_T << "\n";
    for (int i = 0; i <= n_T; i++)
    {
        std::cout << i << "\n";
        evolve_CN();
        p_val = arma::vec(arma::real(u % arma::conj(u)));
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
    std::cout << n_T << M << "SIZE " << arma ::size(u_save);
    u_save.save("data/double_slit_sigmay_02.dat");
    u_real.save("data/double_slit_sigmay_02_real.dat");
    u_imag.save("data/double_slit_sigmay_02_imag.dat");
}
// Sets up initial state vector u0 based on an unnormalised Gaussian wave packet epression
void DoubleSlitBox::init_wave()
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
            // not sure if the below compelx nr esp works
            wp_Re = -(std::pow(x - xc, 2) / (2 * std::pow(sigma_x, 2))) - (std::pow(y - yc, 2) / (2 * std::pow(sigma_y, 2)));
            wp_Im = px * (x - xc) + py * (y - yc);
            wp_val = std::exp(std::complex<double>(wp_Re, wp_Im));
            // std::cout << wp_val * std::conj(wp_val) << "\n";
            sum += std::real(wp_val * std::conj(wp_val));
            u(translate_indices(i, j)) = wp_val; // * std::pow(10, 20);
            // std::cout << u(translate_indices(i, j)) << "\n";
        }
    }
    u /= std::sqrt(sum);
    // u.print();
    // std::cout << arma::sum(arma::real(u % arma::conj(u)));
    // SOMETHING WRONG WITH THE TEST FOR THE SUM OVER U LENGTH
    // std::cout << sum << "\n";
    // arma::vec im = arma::imag(u);
    // std::cout << im * im << "test\n";
    // std::cout << arma::real(u) * arma::real(u) + arma::imag(u) * arma::imag(u);

    // std::cout << u(1) << "vsd" << u;

    // arma::cx_vec u_conj = arma::conj(u);
    // std::cout << u << wp_Re << "\n\nTEST\n";

    // std::cout << (arma::conj(u));
    //  NB! u needs to be normalized !!!!
}

// initializes the potential V matrix
void DoubleSlitBox::init_V()
{ // double wall_width, double wall_center_x, wall_center_y
    // std::cout << V.size();

    double wall_width = 0.02;
    double wall_center_x = 0.5;
    double wall_center_y = 0.5;
    double slit_distance = 0.05;
    double slit_aperture = 0.05;
    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < M; i++)
        {
            if ((i * h >= wall_center_x - wall_width / 2) && (i * h <= wall_center_x + wall_width / 2))
            {
                if (j * h <= wall_center_y - slit_distance / 2 - slit_aperture || j * h >= wall_center_y + slit_distance / 2 + slit_aperture || (j * h <= wall_center_y + slit_distance / 2 && j * h >= wall_center_y - slit_distance / 2))
                {
                    V(i, j) = v0;
                }
            }
            else
            {
                // std::cout << i << j << V.size() << "\n";

                V(i, j) = 0;
                // std::cout << i << j << V.size() << "\n";
            }
            // std::cout << i << j << V(i, j) << "\n";
        }
    }

    V.save("data/V.dat");
}