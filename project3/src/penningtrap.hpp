// PenningTrap class
extern const double k_e;

#ifndef __penningtrap_hpp__
#define __tenningtrap_hpp__

#include <string>
#include <armadillo>
#include <time.h>
#include "particle.hpp"

class PenningTrap
{
  private:
    double V_0_time_dep;
    arma::vec E_ext;
    double E_x;
    double E_y;
    double E_z;
    arma::vec r_diff;
    arma::vec E_i;
    arma::vec F_int;
    arma::vec B;
    arma::vec F_ext;
    arma::vec F_p;
    arma::vec F_tot;

    arma::vec tmp_pos;
    arma::vec tmp_vel;
    arma::vec kr_1;
    arma::vec kv_1;
    arma::vec kr_2;
    arma::vec kv_2;
    arma::vec kr_3;
    arma::vec kv_3;
    arma::vec kr_4;
    arma::vec kv_4;
    int counter;

  public:
    // Member variables
    int n; // Number of particles 
    double dt2; //perform 4n less FLOPs by cutting divide in k2 and k3
    double dt6; //perform 2n less FLOPs by cutting divide in evolve
    double B_0; // magnetic field strength
    double V_0; // applied potetial
    double d; // characteristic dimension
    bool interaction;
    bool time_dependency; 
    double V_0_d2; // V_0 / (d*d) for time independency
    double time = 0; // time value to be updated
    double f; // amplitude_in
    double omega_v; // frequency_in

    
    std::vector<Particle> p; // conntains Particle objects in the Penning trap

    // Constructor
    PenningTrap(double B0_in, double V0_in, double d_in, bool interaction_in, bool time_dependency_in);

    // Add a particle to the trap
    void add_particle(Particle p_in);

    // Add n random particles to the trap
    void add_n_random_particles(int n);

    // External electric field at point r=(x,y,z)
    arma::vec external_E_field(arma::vec r);

    //Method to set amplitude f and angular frequency omega_v
    void set_amplitude_and_frquency(double amplitude_in, double frequency_in);

    // External electric field at point r=(x,y,z) with time dependence
    arma::vec external_E_field(arma::vec r, double t);

    // External magnetic field at point r=(x,y,z)
    arma::vec external_B_field(arma::vec r);

    // Force on particle_i from particle_j
    arma::vec force_particle(int i, int j);

    // The total force on particle_i from the external fields
    arma::vec total_force_external(int i);

    // The total force on particle_i from the other particles
    arma::vec total_force_particles(int i);

    // The total force on particle_i from both external fields and other particles
    arma::vec total_force(int i);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order
    void evolve_RK4(double dt);

    // Evolve the system one time step (dt) using Forward Euler
    void evolve_forward_Euler(double dt);

    // Takes in t vector with time steps, returns analytic solution of single particle motion over time as vectros x, y, z. 
    arma::mat analytic(arma::vec t);

    //Returns number of particles still inside trap
    int particles_left_in_trap();

    //Creates a matrix with fraction of particles left in trap for different frequencies and writes it to a .dat file.
    void parcticles_left_for_omega_v(double dt, int N, double omega_min, double omega_max, double omega_step);

};

#endif
