#include "slit_box.hpp"
// g++ main.cpp double_slit_box.cpp -o main -larmadillo -O2

int main(int argc, char *argv[])
{

    double h_in = atof(argv[1]);
    double dt_in = atof(argv[2]);
    double T_in = atof(argv[3]);
    double xc_in = atof(argv[4]);
    double sigma_x_in = atof(argv[5]);
    double px_in = atof(argv[6]);
    double yc_in = atof(argv[7]);
    double sigma_y_in = atof(argv[8]);
    double py_in = atof(argv[9]);
    int v0_in = atoi(argv[10]);
    int slits = atoi(argv[11]);
    std::string savename = argv[12];

    SlitBox DS = SlitBox(h_in, dt_in, T_in, xc_in, sigma_x_in, px_in, yc_in, sigma_y_in, py_in, v0_in);
    if (slits > 0)
    {
        DS.init_V(slits);
    }

    DS.fill_A_B();
    DS.init_wave();
    DS.evolve_CN_to_time();
    arma::cube p_val = DS.u_save;
    arma::cube u_real = DS.u_real;
    arma::cube u_imag = DS.u_imag;
    p_val.save("data/" + savename + ".dat");
    u_real.save("data/" + savename + "_real" + ".dat");
    u_imag.save("data/" + savename + "_imag" + ".dat");

    return 0;
}
